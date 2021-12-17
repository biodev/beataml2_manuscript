clinical.for.genomic.analysis <- function(clin.wkst){
  
  clin.dt <- data.table(read.xlsx(clin.wkst))
  
  #specimenGroup filtered on text contains "Initial". I got 582 out of 913 AML/ALALs for 63.75%.
  
  clin.dt[,isInitial:=grepl("Initial", specimenGroups)]
  clin.dt[,isTherapy:=SpecificDxAtAcquisition=="Therapy-related myeloid neoplasms"]
  
  stopifnot(clin.dt[dxAtSpecimenAcquisition %in% c("ACUTE MYELOID LEUKAEMIA (AML) AND RELATED PRECURSOR NEOPLASMS", "ACUTE LEUKAEMIAS OF AMBIGUOUS LINEAGE") & isInitial==T,.N] == 582)
  
  #limit to one sample per datatype
  
  use.clin <- clin.dt[used_manuscript_analyses=="yes"]
  
  #make a more R-friendly patient name
  use.clin[,ptid:=paste0("pt", dbgap_subject_id)]
  
  #not sure what best to do with duplicate patients?  Send down the line given that the samples are not duplicated?
  
  use.clin[ELN2017 %in% c("FavorableOrIntermediate", "MissingMutations", "NonInitial", "MissingKaryo", "NonAML", "IntermediateOrAdverse"), ELN2017:=NA_character_]
  
  use.clin[,age_cut:=cut(ageAtDiagnosis, labels=c("young", "middle", "older", "oldest"), breaks=c(-1, 45, 60, 75, 100))]
  
  #as the other diagnoses are reported differently
  use.clin[dxAtSpecimenAcquisition %in% c("ACUTE MYELOID LEUKAEMIA (AML) AND RELATED PRECURSOR NEOPLASMS", "ACUTE LEUKAEMIAS OF AMBIGUOUS LINEAGE") == F, `:=`(overallSurvival=NA_real_, vitalStatus=NA_character_, isInitial=NA)]
  
  use.clin[vitalStatus=="Unknown",vitalStatus:=NA_character_]
  use.clin[overallSurvival <= 0, overallSurvival:=NA_real_]
  
  use.clin[,isDead:=vitalStatus == "Dead"]
  
  use.clin
  
}

clinical.for.summary <- function(clin){
  
  #also another version that strictly has one sample per patient for clinical-only analysis
  #choose the first sample in these cases (the ones I've seen were pb/bm with sample dates so shouldn't make a difference)
  
  clin.sum <- clin[,.SD[order(timeOfSampleCollectionRelativeToInclusion)][1],by=dbgap_subject_id]
  
  stopifnot(clin.sum[,.N,by=dbgap_subject_id][,all(N==1)])
  
  #make sure no patients overlap between the wave cohorts
  stopifnot(clin.sum[,length(unique(cohort)),by=dbgap_subject_id][,all(V1 == 1)])
  #or alternatively
  stopifnot(length(intersect(clin.sum[cohort == "Waves1+2",dbgap_subject_id], clin.sum[cohort == "Waves3+4",dbgap_subject_id]))==0)
  
  clin.sum
  
}


parse.van.galen.supp <- function(vg.supp.file){
  
  vg <- data.table(read.xlsx(vg.supp.file))
  vg <- vg[,-c(1:4, 14),with=F]
  names(vg) <- unlist(vg[1,])
  names(vg)[1:3] <- paste(names(vg)[1:3], "Combined", sep="-")
  vg <- vg[-1,]
  vg <- vg[-51,]
  vg[,rank:=.I]
  
  vg[,`GMP-like`:=`GMP-like-Combined`]
  
  melt.vg <- melt(vg, id.vars="rank", measure.vars=setdiff(colnames(vg), "rank"), variable.factor = F)
  
  names(melt.vg) <- c("vg_rank", "vg_type", "display_label")
  
  melt.vg
  
}

get.combined.exprs <- function(exprs.file, clin){
  
  exprs.dt <- fread(exprs.file, sep="\t", header=T)
  
  exprs.mat <- as.matrix(exprs.dt[,-c(1:4),with=F])
  rownames(exprs.mat) <- exprs.dt$display_label
  
  rna.clin <- clin[manuscript_rnaseq == "yes"]
  
  exprs.mat <- exprs.mat[,rna.clin$dbgap_rnaseq_sample]
  colnames(exprs.mat) <- rna.clin$ptid
  
  t(exprs.mat)
}

mutation.only.dataset <- function(mut.file, clin){
  
  muts <- fread(mut.file, sep="\t", header=T)
  
  dna.clin <- clin[manuscript_dnaseq == "yes"]
  
  muts.m <- merge(dna.clin[,.(ptid, dbgap_sample_id=dbgap_dnaseq_sample)], muts, by="dbgap_sample_id")
  muts.m[,value:=1]
    
  muts.m[,pt_fac:=factor(ptid, levels=dna.clin$ptid)]
  
  comb.mut <- reshape2::acast(pt_fac~symbol, value.var="value", data=muts.m, fun.aggregate = function(x) as.integer(sum(x)>0), fill=0, drop=F)
  
  #need to fill in the manual entries, npm1, flt3 and the three below
  
  cons.mut <- as.matrix(dna.clin[,.(NPM1=as.integer(NPM1=="positive"), `FLT3_ITD`=as.integer(`FLT3-ITD`=="positive"),
                                    RUNX1=as.integer(is.na(RUNX1)==F), ASXL1=as.integer(is.na(ASXL1)==F), 
                                    TP53=as.integer(is.na(TP53)==F ))])
  
  rownames(cons.mut) <- dna.clin$ptid
  
  comb.mut <- cbind(comb.mut,FLT3_ITD=0)
  
  comb.mut[,c("NPM1", "FLT3_ITD", "RUNX1", "ASXL1", "TP53")] <- cons.mut[rownames(comb.mut), c("NPM1", "FLT3_ITD", "RUNX1", "ASXL1", "TP53")]
  
  stopifnot(sum(is.na(comb.mut))==0)
  
  #split into waves1+2 (train) vs waves3+4 (test)
  
  train.mut <- comb.mut[dna.clin[cohort == "Waves1+2"]$ptid,]
  train.mut <- train.mut[,colSums(train.mut) > 0]
  
  #here for consistency, only keep the genes called on the custom capture panel (as a couple patients were sequenced using nextera)
  test.mut <- comb.mut[dna.clin[cohort == "Waves3+4"]$ptid,]
  test.mut <- test.mut[,colSums(test.mut) > 0]
  test.mut <- test.mut[,c(muts.m[capture_type == "Nimblegen Custom Capture",unique(symbol)], "FLT3_ITD")]
  
  #limit the combined set to only those genes in common
  
  comb.mut <- comb.mut[,intersect(colnames(train.mut), colnames(test.mut))]
  
  list(comb=comb.mut, train=train.mut, test=test.mut)
}

inhibitor.dataset <- function(inh.file, clin){
  
  inh <- fread(inh.file, sep="\t", header=T)
  
  inh[,base_samp:=ifelse(dbgap_dnaseq_sample != "", sub("D", "", dbgap_dnaseq_sample), sub("R", "", dbgap_rnaseq_sample))]
  clin[,base_samp:=ifelse(is.na(dbgap_dnaseq_sample)==F, sub("D", "", dbgap_dnaseq_sample), sub("R", "", dbgap_rnaseq_sample))]
  
  inh.clin <- clin[manuscript_inhibitor=="yes"]
  
  stopifnot(inh.clin[,.N,by=dbgap_subject_id][,all(N==1)])
  
  use.inh <- merge(inh.clin[,.(dbgap_subject_id, base_samp, ptid)], inh, by=c("dbgap_subject_id", "base_samp"))
  
  use.inh
}

get.family.list <- function(family.file){
  
  sheets <- getSheetNames(family.file)
  lapply(setNames(sheets, sheets), function(x){
    data.table(read.xlsx(family.file, sheet=x))
  })
}

predict.modules <- function(comb.exprs, clin, wgcna.maps){
  
  train.exprs <- comb.exprs[clin[manuscript_rnaseq == "yes" & cohort == "Waves1+2"]$ptid,]
  
  train.prs <- lapply(split(wgcna.maps$cur.map, by="cur_labels"), function(x){
    
    z.exprs <- scale(train.exprs[,x$display_label])
    
    av.expr <- rowMeans(z.exprs)
    
    tmp.prs <- prcomp(z.exprs, center=F, scale.=F)
    
    tmp.me <- data.table(ptid=rownames(tmp.prs$x), tmp.prs$x[,paste0("PC", 1:5)])
    
    if(sign(cor(tmp.prs$x[names(av.expr),"PC1"], av.expr)) < 0){
      
      tmp.me$PC1 <- -tmp.me$PC1
      
    }
    
    list(mes=tmp.me, pca=tmp.prs, centers=attr(z.exprs,"scaled:center"), scales=attr(z.exprs,"scaled:scale"))
    
  })
  
  train.mes <- rbindlist(lapply(train.prs, function(x){
    
    x$mes
    
  }), idcol="module")
  
  
  #project the test data as well
  
  test.exprs <- comb.exprs[clin[manuscript_rnaseq == "yes" & cohort == "Waves3+4"]$ptid,]
  
  test.mes <- rbindlist(lapply(train.prs, function(x){
    
    z.exprs <- scale(test.exprs[,names(x$centers)], center=x$centers, scale=x$scales)
    
    av.expr <- rowMeans(z.exprs)
    
    tmp.preds <- z.exprs[,rownames(x$pca$rotation)] %*% x$pca$rotation
    
    tmp.me <- data.table(ptid=rownames(tmp.preds), tmp.preds[,paste0("PC", 1:5)])
    
    if(sign(cor(tmp.preds[names(av.expr),"PC1"], av.expr)) < 0){
      
      tmp.me$PC1 <- -tmp.me$PC1
      
    }
    
    tmp.me
    
  }), idcol="module")
  
  comb.mes <- rbind(cbind(train.mes, cohort="Waves 1+2"),
                    cbind(test.mes, cohort="Waves 3+4"))
  
  comb.mes
  
}

van.galen.scores <- function(comb.exprs, clin, vg.genes){
  
  train.exprs <- comb.exprs[clin[manuscript_rnaseq == "yes" & cohort == "Waves1+2"]$ptid,]
  
  rna.clin <- split(clin[manuscript_rnaseq=="yes"], by="cohort")
  
  stopifnot(length(rna.clin) == 2)
  
  #van galen
  
  vg.genes <- vg.genes[display_label %in% colnames(train.exprs) & vg_rank %in% 1:30]
  use.vg <- vg.genes[grepl("Combined", vg_type)==F]
  
  #center and scale separately between BM and PB to mitigate distributional differences
  
  #start with waves1/2
  bm.scale.exprs <- t(scale(train.exprs[rna.clin[["Waves1+2"]][specimenType == "Bone Marrow Aspirate",ptid],]))
  pb.scale.exprs <- t(scale(train.exprs[rna.clin[["Waves1+2"]][specimenType != "Bone Marrow Aspirate",ptid],]))
  scale.exprs <- cbind(bm.scale.exprs, pb.scale.exprs[rownames(bm.scale.exprs),])
  
  split.vg <- split(use.vg, by="vg_type")
  
  eigct.list <- lapply(split.vg, function(x){
    
    tmp <- prcomp(t(scale.exprs[x$display_label,]), center=F, scale=F)
    
    #similar to WGCNA MEs, align with average expression
    
    av.expr <- colMeans(scale.exprs[x$display_label,])
    
    if(sign(cor(tmp$x[,"PC1"], av.expr)) < 0){
      
      tmp.dt <- data.table(ptid=rownames(tmp$x), PC1=-tmp$x[,"PC1"], PC2=tmp$x[,"PC2"])
      
    }else{
      tmp.dt <- data.table(ptid=rownames(tmp$x), PC1=tmp$x[,"PC1"], PC2=tmp$x[,"PC2"])
    }
    
    #roughly group
    tmp.dt[,groups:="intermediate"]
    tmp.dt[PC1 < quantile(PC1, .25), groups:="low"]
    tmp.dt[PC1 > quantile(PC1, .75), groups:="high"]
    
    #cors with eigengene
    
    kme.mat <- cor(t(scale.exprs[x$display_label,])[tmp.dt$ptid,], tmp.dt$PC1)
    kme.dt <- data.table(display_label=rownames(kme.mat), kme=kme.mat[,1])
    kme.dt <- merge(x, kme.dt, by="display_label")
    
    #also record eigenvector
    
    kme.dt <- cbind(kme.dt, ev=tmp$rotation[kme.dt$display_label,"PC1"])
    
    list(summary=tmp.dt, pcs=tmp, kme=kme.dt)
  })
  
  test.exprs <- comb.exprs[clin[manuscript_rnaseq == "yes" & cohort == "Waves3+4"]$ptid,]
  
  wv34.bm.scale.exprs <- t(scale(test.exprs[rna.clin[["Waves3+4"]][specimenType == "Bone Marrow Aspirate",ptid],], center=attr(bm.scale.exprs,"scaled:center"), scale=attr(bm.scale.exprs,"scaled:scale")))
  wv34.pb.scale.exprs <- t(scale(test.exprs[rna.clin[["Waves3+4"]][specimenType != "Bone Marrow Aspirate",ptid],], center=attr(pb.scale.exprs,"scaled:center"), scale=attr(pb.scale.exprs,"scaled:scale")))
  
  wv34.scale.exprs <- cbind(wv34.bm.scale.exprs, wv34.pb.scale.exprs[rownames(wv34.bm.scale.exprs),])
  
  wv34.eigct.list <- lapply(setNames(names(split.vg), names(split.vg)), function(x){
    
    tmp <- t(wv34.scale.exprs[rownames(eigct.list[[x]]$pcs$rotation),]) %*% eigct.list[[x]]$pcs$rotation
    
    #similar to WGCNA MEs, align with average expression
    
    av.expr <- colMeans(wv34.scale.exprs[rownames(eigct.list[[x]]$pcs$rotation),])
    
    if(sign(cor(tmp[,"PC1"], av.expr)) < 0){
      
      tmp.dt <- data.table(ptid=rownames(tmp), PC1=-tmp[,"PC1"], PC2=tmp[,"PC2"])
      
    }else{
      tmp.dt <- data.table(ptid=rownames(tmp), PC1=tmp[,"PC1"], PC2=tmp[,"PC2"])
    }
    
    #roughly group
    tmp.dt[,groups:="intermediate"]
    tmp.dt[PC1 < quantile(PC1, .25), groups:="low"]
    tmp.dt[PC1 > quantile(PC1, .75), groups:="high"]
    
    list(summary=tmp.dt)
  })
  
  #Make sure the eigencell-types are similar enough
  
  eigct.dt <- rbind(
    cbind(rbindlist(lapply(eigct.list, "[[", "summary"), idcol="vg_type"), cohort="Waves1+2"),
    cbind(rbindlist(lapply(wv34.eigct.list, "[[", "summary"), idcol="vg_type"), cohort="Waves3+4")
  )
  
  eigct.dt <- merge(eigct.dt, clin[manuscript_rnaseq=="yes",.(ptid, cohort, specimenType)], by=c("ptid", "cohort"))
  
  list(summary=eigct.dt, wv12_results=eigct.list, wv12_sexprs=scale.exprs)
  
}

mutation.expression.features <- function(vg.scores, comb.mes, mut.list, clin){
  
  #first modules
  me.dt <- dcast(ptid~module, value.var="PC1",data=comb.mes[module != "M0"])
  
  me.dt <- merge(me.dt, comb.mes[module == "M0",.(ptid, M0.PC1=PC1, M0.PC2=PC2, M0.PC3=PC3, M0.PC4=PC4, M0.PC5=PC5 )], by="ptid")
  
  #then cell-types
  
  #figure out correlations, remove those that are too correlated with cell-types (.75)
  
  test.cors <- merge(melt(me.dt, id.vars=c("ptid"), variable.factor=F), vg.scores$summary, by=c("ptid"), allow.cartesian=TRUE)
  
  #only remove those that have both a positive and negative corr at this level
  mods.to.rm <- test.cors[,.(cor(value, PC1)),by=.(variable, vg_type)][abs(V1) > .6,.(length(unique(sign(V1)))),by=.(variable)][V1 > 1]
  
  me.dt <- me.dt[,names(me.dt) %in% mods.to.rm$variable == F,with=F]
  
  #now combine with cell-type
  
  ect.mat <- dcast(ptid~vg_type, value.var="PC1",data=vg.scores$summary)
  
  cm.dt <- merge(me.dt, ect.mat, by="ptid")
  
  #consensus fusions
  
  clin[is.na(consensusAMLFusions),consensusAMLFusions:=""]
  
  #the N/A here indicates those without rnaseq or karyotypes so are universally missing
  cons.pts <- clin[consensusAMLFusions != "N/A"]
  
  uniq.pts <- cons.pts[,.(consensusAMLFusions=paste(unique(consensusAMLFusions), collapse=";"),.N), by=.(ptid)]
  uniq.pts[,ConsensusFusion:=make.names(paste("mut",sub(";", "", consensusAMLFusions), sep="."))]
  uniq.pts[,value:=1]
  uniq.pts[,ptfac:=factor(ptid)]
  
  num.fus <- uniq.pts[,.N,by=ConsensusFusion]
  
  uniq.pts <- uniq.pts[ConsensusFusion %in% num.fus[N < 10,ConsensusFusion]==F]
  
  fusion.mat <- dcast(ptfac~ConsensusFusion, value.var="value", data=uniq.pts[ConsensusFusion != "mut."], fill=0,drop=F)
  fusion.mat[,ptid:=as.character(ptfac)]
  fusion.mat[,ptfac:=NULL]
  
  cmf.dt <- merge(cm.dt, fusion.mat, by="ptid", all=T)
  
  #limit to those with 10 mutations in train.mut
  use.genes <- colnames(mut.list$train)[colSums(mut.list$train) >= 10]
  
  #make sure in test
  stopifnot(length(setdiff(use.genes, colnames(mut.list$test)))==0)
  
  comb.mut.dt <- data.table(ptid=rownames(mut.list$comb), mut=mut.list$comb[,use.genes])
  
  cmm.dt <- merge(cmf.dt, comb.mut.dt, by="ptid", all=T)
  
  names(cmm.dt) <- make.names(names(cmm.dt))
  
  cmm.dt
}

drug.family.enrichment <- function(family.list, inhib){
  
  families <- family.list$drug_family
  fam.list <- split(families$inhibitor, families$family)
  
  inhib[,perc_max:=auc/((log10(max_conc)-log10(min_conc))*100)]
  
  inhib.mat <- reshape2::acast(inhibitor~ptid, value.var="perc_max", data=inhib)
  
  ssg.dt <- rbindlist(lapply(colnames(inhib.mat), function(x){
    
    tmp.mat <- tryCatch(suppressWarnings(gsva(inhib.mat[is.na(inhib.mat[,x])==F,x,drop=F], fam.list, method="ssgsea", ssgsea.norm=F, min.sz=5)), error = function(e) NULL)
    
    if (is.null(tmp.mat)){
      
      tmp.dt <- data.table(Var1=character(0), Var2=character(0), value=numeric(0))
      
    }else{
      tmp.dt <- data.table(reshape2::melt(data=tmp.mat, as.is=T))
    }
    
    tmp.dt[,.(family=Var1, ptid=Var2, ssgsea=value)]
  }))
  
  #remove those families with few samples
  
  ssg.dt[,num_samps:=.N,by=family]
  
  ssg.dt[num_samps > 50]
  
}

pear1.rna.clinical <- function(exprs, clin){
  
  rna.clin <- clin[manuscript_rnaseq=="yes"]
  
  p1.dt <- cbind(rna.clin[,.(ptid, overallSurvival, isDead, specimenType, ELN2017, age_cut, isInitial, isDenovo, isRelapse, isTransformed, isTherapy)], PEAR1=exprs[rna.clin$ptid,"PEAR1"])
  p1.dt[,tissue:=ifelse(specimenType == "Bone Marrow Aspirate", "BM", "PB/Leuk")]
  
  p1.dt
  
}

compute.lsc17 <- function(exprs){
  
  #from https://www.nature.com/articles/nature20598#Tab8
  #(DNMT3B×0.0874) + (ZBTB46×−0.0347) + (NYNRIN×0.00865) + (ARHGAP22×−0.0138) + (LAPTM4B×0.00582) + (MMRN1×0.0258) + (DPYSL3×0.0284) + 
  #(KIAA0125×0.0196) + (CDK6×−0.0704) + (CPXM1×−0.0258) + (SOCS2×0.0271) + (SMIM24×−0.0226) + (EMP1×0.0146) + (NGFRAP1×0.0465) + 
  #(CD34×0.0338) + (AKR1C3×−0.0402) + (GPR56×0.0501). 
  #As above- and below-median scores in the training cohort were associated with adverse and favourable cytogenetic risk, respectively, a median threshold was used to discretize scores into high and low groups.
  
  lsc17.coefs <- c(DNMT3B=0.0874, ZBTB46=-0.0347, NYNRIN=0.00865, ARHGAP22=-0.0138, LAPTM4B=0.00582, MMRN1=0.0258, DPYSL3=0.0284, 
                   KIAA0125=0.0196, CDK6=-0.0704, CPXM1=-0.0258, SOCS2=0.0271, SMIM24=-0.0226, EMP1=0.0146, NGFRAP1=0.0465, CD34=0.0338, 
                   AKR1C3=-0.0402,GPR56=0.0501)
  
  common.genes <- intersect(names(lsc17.coefs), colnames(exprs))#16
  
  stopifnot(length(common.genes) == 16)
  stopifnot(all(setdiff(names(lsc17.coefs), colnames(exprs)) == "SMIM24"))
  
  #smim24 has an ensembl id of ENSG00000095932 so maps to C19orf77
  
  names(lsc17.coefs)[names(lsc17.coefs) == "SMIM24"] <- "C19orf77"
  
  common.genes <- intersect(names(lsc17.coefs), colnames(exprs))#17
  
  stopifnot(length(common.genes) == 17)
  
  lsc17.score <- exprs[,names(lsc17.coefs)] %*% lsc17.coefs
  
  #median determines low vs high
  
  lsc17.med <- median(lsc17.score)
  
  #should compute median per cohort
  lsc.dt <- data.table(ptid=rownames(lsc17.score), LSC17=lsc17.score[,1])
  
  lsc.dt
  
}


