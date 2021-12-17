mutation.freq.by.cohort <- function(mut.list){
  
  train.mut <- mut.list$train
  test.mut <- mut.list$test
  
  #wave 1/2 data
  train.dt <- data.table(reshape2::melt(train.mut, as.is=T))[value > 0]
  train.dt[Var2=="FLT3_ITD", Var2:="FLT3"]
  train.freq <- train.dt[,.(num_pats=length(unique(Var1))),by=.(symbol=Var2)]
  train.freq[,prop:=num_pats/nrow(train.mut)]
  
  #retain those > 5%
  train.5p <- train.freq[prop > .05][order(-prop)]
  
  #waves 3/4 data
  test.dt <- data.table(reshape2::melt(test.mut, as.is=T))[value > 0]
  test.dt[Var2=="FLT3_ITD", Var2:="FLT3"]
  test.freq <- test.dt[,.(num_pats=length(unique(Var1))),by=.(symbol=Var2)]
  test.freq[,prop:=num_pats/nrow(test.mut)]
  
  test.5p <- test.freq[prop > .05][order(-prop)]
  
  #now come up with overall ordering based on the parallel prop max
  gene.univ <- data.table(symbol=union(train.5p$symbol, test.5p$symbol))
  
  comb.freq <- merge(gene.univ, train.freq[,.(symbol, train_prop=prop)], by="symbol",all.x=T, all.y=F)
  
  comb.freq <- merge(comb.freq, test.freq[,.(symbol, test_prop=prop)], by="symbol",all.x=T, all.y=F)
  
  comb.freq[,max_prop:=pmax(train_prop, test_prop, na.rm=T)]
  
  #this is the overall ordering
  comb.freq <- comb.freq[order(max_prop)]
  
  comb.freq[,symb_ord:=factor(symbol, levels=symbol, ordered=T)]
  
  comb.freq[,diff:=test_prop - train_prop]
  
  comb.freq
  
}

summarize.denovo.aucs <- function(inhib, clin.dt){
  
  inhib.clin <- clin.dt[manuscript_inhibitor=="yes"]
  
  inhib.tt <- inhib[status == "train/test"]
  
  inhib.tt.clin <- merge(inhib.tt, inhib.clin[,.(ptid, cohort, isDenovo)], by="ptid")
  
  stopifnot(inhib.tt.clin[,.N] == inhib.tt[,.N])
  
  inhib.tt.clin.sum <- inhib.tt.clin[isDenovo=="TRUE",.(median=median(auc), mean=mean(auc), sd=sd(auc), n=.N), by=.(inhibitor, cohort)]
  
  inhib.tt.clin.ef <- dcast(inhibitor~cohort, value.var=c("mean","median") ,data=inhib.tt.clin.sum)
  
  inhib.tt.clin.ef
}

vg.muts.assocs <- function(features){
  
  mut.ct <- as.matrix(features[,grepl("mut\\.", names(features)) | grepl("\\.like", names(features)),with=F])
  rownames(mut.ct) <- features$ptid
  
  mut.ct <- mut.ct[complete.cases(mut.ct),]
  
  mut.cols <- colnames(mut.ct)[grep("mut\\.", colnames(mut.ct))]
  
  ct.cols <- colnames(mut.ct)[grep("\\.like", colnames(mut.ct))]
  
  mut.ct.assocs <- .run.assoc.welch(mut.ct[,mut.cols], mut.ct[,ct.cols], min.muts=5)
  
  #mut.ct.assocs[,hist(pval)]
  mut.ct.assocs[,qval:=qvalue(pval)$qvalues]
  
  mut.ct.assocs
}

write.vma <- function(mut.ct.assocs){
  
  write.xlsx(list(ct_mut_assoc=mut.ct.assocs[,.(cell_type=inhibitor, gene, num_muts, num_neg, estimate, glass_d, t, pval, qval)]),file="output_files/Table_S2.xlsx")
  
  return("output_files/Table_S2.xlsx")
  
}

drug.ct.cors <- function(inhib, vg.scores){
  
  tmp <- merge(vg.scores$summary, inhib[status %in% c("train/test","combined")], by="ptid", all=F, allow.cartesian=T)
  
  all.cors <- tmp[,tidy(cor.test(auc, PC1)),by=.(inhibitor, vg_type)]
  
  #as the pvalues seem to indicate a little too much dependency
  all.cors[,fdr:=p.adjust(p.value, method="BY")]
  
  all.cors
}

inhibitor.mutation.interactions <- function(inhib, features){
  
  inhib.mat <- reshape2::acast(ptid~inhibitor, value.var="auc", data=inhib[status %in% c("train/test","combined")])
  
  feat.mat <- as.matrix(features[,-1,with=F])
  rownames(feat.mat) <- features[[1]]
  feat.mat <- feat.mat[complete.cases(feat.mat),]
  
  inter.res <- .run.drug.mut.ct.inters(feat.mat, inhib.mat, mut.re="mut\\.", ct.re="\\.like")
  
  inter.res[,qval:=qvalue(pval)$qvalues]
  
  inter.res
}

write.cor.ints <- function(inters, cors){
  
  top.inters <- inters[,.SD[order(pval)][1],by=.(inhibitor, vg_type=sub("\\.", "-", ct))]
  
  cors.m <- merge(cors[,.(inhibitor, cell_type=vg_type, cor=estimate, cor_fdr=fdr)], top.inters[,.(inhibitor, cell_type=vg_type, inter_gene=gene, inter_est=estimate, inter_muts=num_muts, inter_qval=qval)], by=c("inhibitor", "cell_type"))
  
  write.xlsx(list(drug_cell_type=cors.m), file="output_files/Table_S3.xlsx")
  
  return("output_files/Table_S3.xlsx")
}

family.mut.assocs <- function(fam.enrich,family.list, features){
  
  fam.syns <- family.list$synonyms
  
  ssg.dt <- merge(fam.enrich, fam.syns, by="family")
  
  inhib.mat <- reshape2::acast(ptid~synonym, value.var="ssgsea", data=ssg.dt)
  
  feat.mat <- as.matrix(features[,grepl("mut\\.", names(features)),with=F])
  rownames(feat.mat) <- features$ptid
  
  feat.mat <- feat.mat[complete.cases(feat.mat),]
  
  comb.assoc <- .run.assoc.welch(feat.mat, inhib.mat, min.muts=5)
  
  comb.assoc[,qval:=qvalue(pval)$qvalues]
  
  comb.assoc
  
}

write.fma <- function(comb.assoc, family.list){
  
  comb.assoc <- merge(family.list$synonyms, comb.assoc, by.x="synonym", by.y="inhibitor")
  
  write.xlsx(list(family_mut_assoc=comb.assoc[,.(family, synonym, gene, num_muts, num_neg, estimate, glass_d, t, pval, qval)]),file="output_files/Table_S5.xlsx")
  
  return("output_files/Table_S5.xlsx")
}

family.ct.cors <- function(fam.enrich, vg.scores){
  
  tmp <- merge(vg.scores$summary, fam.enrich, by="ptid", all=F, allow.cartesian=T)
  
  all.cors <- tmp[,tidy(cor.test(ssgsea, PC1)),by=.(family, vg_type)]
  
  #as the pvalues seem to indicate a little too much dependency
  all.cors[,fdr:=p.adjust(p.value, method="BY")]
  
  all.cors
}


family.mutation.interactions <- function(fam.enrich, features){
  
  inhib.mat <- reshape2::acast(ptid~family, value.var="ssgsea", data=fam.enrich)
  
  feat.mat <- as.matrix(features[,-1,with=F])
  rownames(feat.mat) <- features[[1]]
  feat.mat <- feat.mat[complete.cases(feat.mat),]
  
  inter.res <- .run.drug.mut.ct.inters(feat.mat, inhib.mat, mut.re="mut\\.", ct.re="\\.like")
  
  inter.res[,qval:=qvalue(pval)$qvalues]
  
  inter.res
}

write.fam.cor.inter <- function(inters, cors, family.list){
  
  top.inters <- inters[,.SD[order(pval)][1],by=.(family=inhibitor, vg_type=sub("\\.", "-", ct))]
  
  cors.m <- merge(cors[,.(family, cell_type=vg_type, cor=estimate, cor_fdr=fdr)], top.inters[,.(family, cell_type=vg_type, inter_gene=gene, inter_est=estimate, inter_muts=num_muts, inter_qval=qval)], by=c("family", "cell_type"))
  
  syn.cors <- merge(family.list$synonyms, cors.m, by="family")
  
  write.xlsx(list(family_cell_type=syn.cors), file="output_files/Table_S6.xlsx")
  
  return("output_files/Table_S6.xlsx")
}

survival.trees <- function(feats, clin){
  
  rna.clin <- clin[manuscript_rnaseq == "yes"]
  
  rna.clin.m <- merge(rna.clin, feats, by="ptid")
  rna.clin.m$value <- 1
  
  age.cast <- dcast(ptid~age_cut, value.var="value", data=rna.clin.m[is.na(age_cut)==F], fun.aggregate=function(x) as.integer(length(x) > 0), fill=0L)
  
  rna.clin.m <- merge(rna.clin.m, age.cast, by="ptid")
  
  rna.clin.df <- as.data.frame(rna.clin.m[,c("ptid", "overallSurvival", "isDead", "young", "middle", "older", "oldest",
                                             setdiff(names(feats), "ptid")), with=F])
  
  rownames(rna.clin.df) <- rna.clin.df$ptid
  rna.clin.df$ptid <- NULL
  
  rna.clin.df <- rna.clin.df[complete.cases(rna.clin.df),]
  
  #works the same with/without minbucket
  s.tree <- ctree(Surv(time=overallSurvival, event=isDead, type="right") ~  ., data = rna.clin.df, control=ctree_control(mincriterion = 0.95, minbucket=20))
  
  sample.nodes <- data.table(ptid=rownames(rna.clin.df))
  sample.nodes[,node:=where(s.tree)]
  
  list(tree=s.tree, nodes=sample.nodes, data=rna.clin.df)
  
}

overall.pear1.trees <- function(p1.dt){
  
  p1.df <- as.data.frame(p1.dt[is.na(overallSurvival)==F & is.na(isDead)==F])
  
  comb.tree <- ctree(Surv(time=overallSurvival, event=isDead, type="right") ~  PEAR1, 
                     data = p1.df, 
                     control=ctree_control(mincriterion = 0.95, minbucket=20, maxdepth=1))
  
  bm.tree <- ctree(Surv(time=overallSurvival, event=isDead, type="right") ~  PEAR1, 
                   data = p1.df[p1.df$tissue == "BM",], 
                   control=ctree_control(mincriterion = 0.95, minbucket=20, maxdepth=1))
  
  pb.tree <- ctree(Surv(time=overallSurvival, event=isDead, type="right") ~  PEAR1, 
                   data = p1.df[p1.df$tissue == "PB/Leuk",], 
                   control=ctree_control(mincriterion = 0.95, minbucket=20, maxdepth=1))
  
  young.df <- p1.df[is.na(p1.df$age_cut) == F & p1.df$age_cut == "young",]
  
  ycomb.tree <- ctree(Surv(time=overallSurvival, event=isDead, type="right") ~  PEAR1, 
                      data = young.df, 
                      control=ctree_control(mincriterion = 0.95, minbucket=20, maxdepth=1))
  
  ybm.tree <- ctree(Surv(time=overallSurvival, event=isDead, type="right") ~  PEAR1, 
                    data = young.df[young.df$tissue == "BM",], 
                    control=ctree_control(mincriterion = 0.95, minbucket=20, maxdepth=1))
  
  ypb.tree <- ctree(Surv(time=overallSurvival, event=isDead, type="right") ~  PEAR1, 
                    data = young.df[young.df$tissue == "PB/Leuk",], 
                    control=ctree_control(mincriterion = 0.95, minbucket=20, maxdepth=1))
  
  list(overall=comb.tree, BM=bm.tree, PB=pb.tree, young_overall=ycomb.tree, young_BM=ybm.tree, young_PB=ypb.tree)
}

combined.kmes <- function(comb.exprs,wgcna.maps, wgcna.mes){
  
  #make a new set of kmes etc using the combined expression set
  
  me.mat <- reshape2::acast(ptid~module,value.var="PC1", data=wgcna.mes[module != "M0",.(ptid, module=sub("M", "ME", module), PC1)])
  
  kme.mat <- WGCNA::signedKME(comb.exprs[,wgcna.maps$cur.map$display_label], me.mat[rownames(comb.exprs),], corFnc="bicor", corOptions = "maxPOutliers = 0.1")
  
  melt.kme <- data.table(reshape2::melt(as.matrix(kme.mat), as.is=T))
  melt.kme[,cur_labels:=sub("kME","M", Var2)]
  
  cur.map <- merge(melt.kme[,.(display_label=Var1, cur_labels, comb_kme=value)], wgcna.maps$cur.map, by=c("display_label", "cur_labels"))
  
  cur.map
  
}

write.mod.membership <- function(exprs.file, wgcna.maps, comb.kme){
  
  exp.dt <- fread(exprs.file, sep="\t", header=T)
  
  all.kme <- merge(wgcna.maps$cur.map[,.(cur_labels, display_label, prev_modules)], comb.kme, by=c("cur_labels", "display_label", "prev_modules"), all=T)
  
  all.kme <- merge(exp.dt[,.(stable_id, display_label, description, biotype)], all.kme, by="display_label")
  
  write.xlsx(list(module_membership=all.kme[order(cur_labels, -comb_kme),.(module_labels=sub("M", "Mod", cur_labels), module_colors=prev_modules, stable_id, symbol=display_label, description, biotype, comb_kme)]), file="output_files/Table_S7.xlsx")
  
  return("output_files/Table_S7.xlsx")
}

detrans.by.vg.by.celltype <- function(clin, vg.scores, inhib, ct.ord){
  
  rna.clin <- clin[manuscript_rnaseq=="yes" & dxAtSpecimenAcquisition %in% c("ACUTE MYELOID LEUKAEMIA (AML) AND RELATED PRECURSOR NEOPLASMS", "ACUTE LEUKAEMIAS OF AMBIGUOUS LINEAGE")==T]
  
  mlt.type <- suppressWarnings(melt(rna.clin, measure.vars=c("isDenovo", "isTransformed"), id.vars=c("ptid", "cohort"), variable.factor=F))
  mlt.sum <- mlt.type[value %in% c("TRUE"),.(type=paste(unique(variable), collapse=";")), by=.(ptid, cohort)]
  mlt.sum <- mlt.sum[grepl(";", type)==F]
  
  vg.sum <- merge(vg.scores$summary[,.(ptid, vg_type, cohort, PC1)], mlt.sum, by=c("ptid", "cohort"))
  
  vg.sum[,cats:=ifelse(type=="isDenovo","Denovo", "Transformed")]
  
  vg.sum[,vg_fac:=factor(vg_type, levels=ct.ord, ordered=T)]
  
  vg.cast <- dcast(ptid+cohort+cats+type~vg_type, value.var="PC1" ,data=vg.sum)
  
  vg.cast.inh <- merge(vg.cast, inhib[,.(ptid, status, inhibitor, auc)], by=c("ptid"), allow.cartesian=T)
  
  list(vg_sum=vg.sum, vg_inh=vg.cast.inh)
  
}


