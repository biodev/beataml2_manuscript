mutation.freq.by.cohort <- function(mut.list, fus.mat){
  
  melt.fus <- melt(id.vars=c("ptid", "cohort"), data=fus.mat, variable.factor=F)
  fus.list <- split(melt.fus, by="cohort")
  
  train.mut <- mut.list$train
  train.fus <- fus.list$`Waves1+2`[value > 0]
  
  test.mut <- mut.list$test
  test.fus <- fus.list$`Waves3+4`[value > 0]
  
  #wave 1/2 data
  train.dt <- data.table(reshape2::melt(train.mut, as.is=T))[value > 0]
  train.dt[Var2=="FLT3_ITD", Var2:="FLT3"]
  train.freq <- train.dt[,.(num_pats=length(unique(Var1))),by=.(symbol=Var2)]
  train.freq[,prop:=num_pats/nrow(train.mut)]
  
  train.fus.freq <- train.fus[,.(num_pats=sum(value)),by=.(symbol=sub("\\.", "-", sub("mut.", "", variable)))]
  train.fus.freq[,prop:=num_pats/fus.mat[cohort == "Waves1+2",.N]]
  
  train.freq <- rbind(cbind(train.freq, type="mut"), cbind(train.fus.freq[,names(train.freq),with=F], type="fus"))
  
  #retain those > 5% or fusions
  train.5p <- train.freq[(prop > .05) | (type == "fus")]
  
  #waves 3/4 data
  test.dt <- data.table(reshape2::melt(test.mut, as.is=T))[value > 0]
  test.dt[Var2=="FLT3_ITD", Var2:="FLT3"]
  test.freq <- test.dt[,.(num_pats=length(unique(Var1))),by=.(symbol=Var2)]
  test.freq[,prop:=num_pats/nrow(test.mut)]
  
  test.fus.freq <- test.fus[,.(num_pats=sum(value)),by=.(symbol=sub("\\.", "-", sub("mut.", "", variable)))]
  test.fus.freq[,prop:=num_pats/fus.mat[cohort == "Waves3+4",.N]]
  
  test.freq <- rbind(cbind(test.freq, type="mut"), cbind(test.fus.freq[,names(test.freq),with=F], type="fus"))
  
  test.5p <- test.freq[(prop > .05) | (type == "fus")]
  
  #now come up with overall ordering based on the parallel prop max
  gene.univ <- data.table(symbol=union(train.5p$symbol, test.5p$symbol))
  
  comb.freq <- merge(gene.univ, train.freq[,.(symbol, train_prop=prop)], by="symbol",all.x=T, all.y=F)
  
  comb.freq <- merge(comb.freq, test.freq[,.(symbol, test_prop=prop)], by="symbol",all.x=T, all.y=F)
  
  comb.freq[,max_prop:=pmax(train_prop, test_prop, na.rm=T)]
  
  #this is the overall ordering
  comb.freq <- comb.freq[order(max_prop)]
  
  comb.freq[is.na(train_prop), train_prop:=0]
  comb.freq[is.na(test_prop), test_prop:=0]
  
  comb.freq[,symb_ord:=factor(symbol, levels=symbol, ordered=T)]
  
  comb.freq[,diff:=test_prop - train_prop]
  
  comb.freq
  
}


train.test.mut.assocs <- function(inhib, mut.list){
  
  inhib.mat <- reshape2::acast(ptid~inhibitor, value.var="auc",data=inhib[status == "train/test"])
  
  #test and train portion
  
  #only run those with >= 5 patients
  
  train.assoc <- .run.assoc.welch(mut.list$train, inhib.mat, min.muts=5)
  
  #only run those with >= 3 patients (as the ratio of train to test is 5/3)
  
  test.assoc <- .run.assoc.welch(mut.list$test, inhib.mat, min.muts=3)
  
  
  tt.assoc <- merge(train.assoc[,.(inhibitor, gene, train_est=estimate, train_se=se, train_t=t, train_pval=pval, train_num=num_inhib, train_muts=num_muts, train_glass=glass_d,
                                   train_glass_var=glass_var)], 
                    test.assoc[,.(inhibitor, gene, test_est=estimate, test_se=se, test_t=t, test_pval=pval, test_num=num_inhib, test_muts=num_muts, test_glass=glass_d,
                                  test_glass_var=glass_var)], 
                    by=c("inhibitor", "gene"), all=F)
  
  #tt.assoc[,hist(train_pval)]#~borderline
  tt.assoc[, train_fdr:=qvalue(train_pval)$qvalues]
  
  #tt.assoc[,hist(test_pval)]
  tt.assoc[, test_fdr:=qvalue(test_pval)$qvalues]
  
  tt.assoc[,sig_cat:="Neither"]
  tt.assoc[train_fdr < .05, sig_cat:="Waves 1+2"]
  tt.assoc[test_fdr < .05, sig_cat:="Waves 3+4"]
  tt.assoc[train_fdr < .05 & test_fdr < .05, sig_cat:="Both"]
  
  tt.assoc
  
}

wgcna.mods.by.drug2 <- function(clin, inhib, wgcna.mes, wgcna.maps){
  
  me.mat <- reshape2::acast(ptid~module, value.var="PC1",data=wgcna.mes)
  
  inhib <- merge(clin[manuscript_inhibitor == "yes",.(ptid, cohort)], inhib, by="ptid")
  
  inh.kme <- rbindlist(lapply(split(inhib[status == "train/test"], by="cohort"), function(x){
    
    inh.mat <- reshape2::acast(ptid~inhibitor, value.var="auc",data=x)
    
    common.pts <- intersect(rownames(inh.mat), rownames(me.mat))
    
    inh.counts <- colSums(is.na(inh.mat[common.pts,])==F)
    
    stopifnot(all(inh.counts > 30) )
    
    cor.inh <- cor(me.mat[common.pts,], inh.mat[common.pts,], use = 'pairwise.complete.obs' )
    
    tmp.melt.cor <- data.table(reshape2::melt(as.matrix(cor.inh), as.is=T, value.name="cor"))
    names(tmp.melt.cor)[1:2] <- c("cur_labels", "inhibitor")
    
    tmp.melt.cor
    
  }), idcol="cohort")
  
  cor.kme <- merge(inh.kme, wgcna.maps$mod.map[cur_labels != "M0"], by="cur_labels")
  
  cor.kme[,titles:=paste0(sub("M", "Mod", cur_labels), " (", prev_modules, ")")]
  
  cast.ck <- dcast(titles+inhibitor~cohort, value.var="cor", data=cor.kme)
  
  cast.ck[,for_col:=(`Waves1+2` + `Waves3+4`)/2]
  
  cast.ck 
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

  rna.clin <- clin[manuscript_rnaseq == "yes" & isDenovo == "TRUE" & is.na(overallSurvival)==F & is.na(isDead)==F]
  
  rna.clin.m <- merge(rna.clin, feats, by="ptid")
  
  rna.clin.m$value <- 1
  
  age.cast <- dcast(ptid~age_cut, value.var="value", data=rna.clin.m[is.na(age_cut)==F], fun.aggregate=function(x) as.integer(length(x) > 0), fill=0L)
  
  rna.clin.m <- merge(rna.clin.m, age.cast, by="ptid")
  
  rna.clin.df <- as.data.frame(rna.clin.m[,c("ptid", "overallSurvival", "isDead", "young", "middle", "older", "oldest",
                                             setdiff(names(feats), "ptid")), with=F])
  
  rownames(rna.clin.df) <- rna.clin.df$ptid
  rna.clin.df$ptid <- NULL
  
  rna.clin.df <- rna.clin.df[complete.cases(rna.clin.df),]
  
  set.seed(3252)
  
  fors.tr <- cforest(Surv(time=overallSurvival, event=isDead, type="right")~., data=rna.clin.df, 
                     control=ctree_control(teststat="quad", testtype="Univ", mincriterion = 0, saveinfo = FALSE),
                     ntree=1000,
                     perturb = list(replace = FALSE, fraction = 0.632))
  
  preds <-  predict(fors.tr, type='response', OOB=T)
  
  preds[is.infinite(preds)] <- max(preds[is.infinite(preds)==F])
  
  test.df <- rna.clin.df[,-c(1:2)]
  
  test.df <- cbind(preds, test.df)
  
  test.tree <- ctree(log10(preds) ~  ., data = test.df, control=ctree_control(mincriterion = 0.95, minbucket = 20L))
  
  overall.s.tree <- ctree(Surv(time=overallSurvival, event=isDead, type="right") ~ M3, data=rna.clin.df, control=ctree_control(mincriterion = 0.95, minbucket = 20L))
  
  list(tree=test.tree, data=rna.clin.df, pred_cor=cor( predict(test.tree), log10(test.df$preds)), overall_Mod3_split=overall.s.tree)
  
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

