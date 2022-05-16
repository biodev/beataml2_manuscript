mut.freq.by.wave.plot <- function(comb.freq){
  
  p1 <- ggplot(data=comb.freq, mapping=aes(x=symb_ord, y=-train_prop)) + geom_col(fill="#023FA5") + coord_flip() +
    scale_x_discrete(drop=F, position = "top") +
    scale_y_continuous(labels=function(x) scales::percent(-x)) +
    theme_bw() + xlab("") + ylab("Wave 1/2 Percentage") +
    theme(
      text=element_text(family="Arial", size=8),
      axis.text.y=element_blank(),
      axis.text.x = element_text(family="Arial", size=8),
      
    )
  
  p2 <- ggplot(data=comb.freq, mapping=aes(x=symb_ord, y=test_prop)) + geom_col(fill="#7D87B9") + coord_flip() +
    scale_y_continuous(labels=scales::percent) +
    theme_bw() + xlab("") + ylab("Wave 3/4 Percentage") +
    theme(
      text=element_text(family="Arial", size=8),
      axis.text.x = element_text(family="Arial", size=8),
      axis.text.y=element_text(hjust=.5, family="Arial", size=8)
    )
  
  #difference
  
  p3 <- ggplot(data=comb.freq, mapping=aes(x=symb_ord, y=diff, fill=ifelse(diff < 0, "low", "high"))) + geom_col() + coord_flip() +
    scale_y_continuous(labels=scales::percent, limits=c(-.15, .15)) +
    scale_fill_manual(values=c(low="blue", high="red"), guide="none") + 
    theme_bw() + xlab("") + ylab("Percentage Difference") +
    theme(
      text=element_text(family="Arial", size=8),
      axis.text.y=element_blank(),
      axis.text.x = element_text(family="Arial", size=8)
    )
  
  ggsave(p1 + p2 + p3, file="figures/beataml_mut_wave_comp_v1.pdf", width=(114/25.4)+.5, height=114/25.4)
  
  return("figures/beataml_mut_wave_comp_v1.pdf")
}

surv.by.wave.plot <- function(clin.sum){
  
  over.wave.dt <- clin.sum[is.na(isDead)==F & is.na(overallSurvival)==F,.(ptid, isDead, overallSurvival, cohort)]
  
  over.wave.fit <- survfit(Surv(time=overallSurvival, event=isDead, type="right") ~  cohort, data = as.data.frame(over.wave.dt)) 
  
  splot <- ggsurvplot(over.wave.fit, data = as.data.frame(over.wave.dt), pval=F, palette=c(`cohort=Waves1+2`="#023FA5", `cohort=Waves3+4`="#7D87B9"), ggtheme=theme_bw())
  
  sv.plot <- splot$plot + theme(
    text=element_text(family="Arial", size=8),
    axis.text = element_text(family="Arial", size=8)
  )
  
  ggsave(sv.plot, file="figures/survival_by_wave_v1.pdf", width=85/25.4, height=85/25.4)
  
  return("figures/survival_by_wave_v1.pdf")  
  
}

denovo.inhib.plot <- function(inhib.sum){
  
  
  mean.clin.plot <- ggplot(data=inhib.sum, mapping=aes(x=`mean_Waves1+2`, y=`mean_Waves3+4`)) + 
    geom_point() + 
    geom_smooth(method="lm", se=F, formula=y~x) +
    theme_bw() + xlab("Waves1+2 Inhibitor AUC") + ylab("Waves3+4 Inhibitor AUC") +
    theme(text=element_text(size=8, family="Arial"),
          axis.text = element_text(family="Arial", size=8))
  
  ggsave(mean.clin.plot, file="figures/beataml_wv1to4_inhib_conc_denovo_only_v1.pdf", width=85/25.4, height=85/25.4)
  
  return("figures/beataml_wv1to4_inhib_conc_denovo_only_v1.pdf")
  
}

train.test.mut.assocs.plot <- function(tt.assoc, wv.cols){
  
  full.glass <- ggplot(data=tt.assoc, mapping=aes(x=train_glass, y=test_glass, shape=sig_cat, alpha=sig_cat, color=sig_cat)) + 
    geom_point(size=3) +
    geom_vline(xintercept=0, linetype="dashed") +
    geom_hline(yintercept=0, linetype="dashed") +
    scale_alpha_manual(values=c(Neither=.25, `Waves 1+2`=.75, `Waves 3+4`=.75, Both=1), name="Significance") +
    scale_shape_manual(values=c(Neither=1, `Waves 1+2`=10, `Waves 3+4`=7, Both=19), name="Significance") +
    scale_color_manual(values=c(wv.cols, Both="black", Neither="grey"), name="Significance") +
    theme_bw() + xlab("Waves 1/2 Effect Size") + ylab("Waves 3/4 Effect Size") +
    theme(
      text=element_text(size=8, family="Arial"),
      legend.text=element_text(size=8, family="Arial"),
      axis.text=element_text(size=8, family="Arial")
    )
  
  ggsave(full.glass, file="figures/beataml_waves1to4_tt_assoc_v2.pdf", width=114/25.4, height=85/25.4)
  
  return("figures/beataml_waves1to4_tt_assoc_v2.pdf")
}

vg.exprs.plot2 <- function(wgcna.mes, wgcna.maps, vg.scores, ct.ord){
  
  ect.mat <- reshape2::acast(ptid~vg_type, value.var="PC1",data=vg.scores$summary)
  
  wgcna.mes <- merge(wgcna.mes, wgcna.maps$mod.map, by.x="module", by.y="cur_labels")
  wgcna.mes[,titles:=paste0(module, " (", prev_modules, ")")]
  
  me.mat <- reshape2::acast(ptid~titles,  value.var="PC1", data=wgcna.mes[module != "M0"])
  
  cor.mat <- t(cor(me.mat, ect.mat[rownames(me.mat),]))
  
  #reorder rows to match differentiation
  cor.mat <- cor.mat[ct.ord,]
  
  #cluster columns
  col.clust <- hclust(dist(t(cor.mat)), method="average")
  
  col.split <- cutree(col.clust, k=3)
  
  stopifnot(all(names(col.split) == colnames(cor.mat)))
  
  #relabel splits
  col.split <- factor(as.character(col.split), levels=c(1, 2, 3), ordered=T)
  
  rc.cols <- brewer.pal(n=11, name="PiYG")
  
  cfun <- circlize::colorRamp2(breaks=c(-1, 0, 1), colors=rc.cols[c(11, 6, 1)])
  
  ht <- Heatmap(cor.mat, col=cfun, cluster_rows=F, cluster_columns=T, column_split=col.split, cluster_column_slices=F,
                heatmap_legend_param = list(title_position = "leftcenter", direction = "horizontal",
                                            title_gp = gpar(fontsize = 8, fontface = "bold", fontfamily="Arial"), 
                                            labels_gp = gpar(fontsize = 8, fontfamily="Arial")), 
                rect_gp = gpar(col = "grey"), border=T, name="Pearsons Corr.", column_title=NULL,
                column_names_gp=gpar(fontfamily="Arial", fontsize=6), 
                row_names_gp=gpar(fontfamily="Arial", fontsize=6))
  
  pdf(file="figures/older_wgcna_heatmap_vg_v2.pdf", width=85/25.4, height=3)
  
  draw( ht, heatmap_legend_side = "bottom")
  
  dev.off()
  
  return("figures/older_wgcna_heatmap_vg_v2.pdf")
}

vg.muts.assoc.plot <- function(mut.ct.assocs, ct.ord, ct.colors){
  
  mut.ct.assocs[,plot_vals:=sign(estimate) * -log10(pval)]
  
  mut.ct.assocs[,fill_class:=NA_character_]
  mut.ct.assocs[estimate < 0 & qval < .05,fill_class:="sig_low"]
  mut.ct.assocs[estimate < 0 & qval >= .05,fill_class:="low"]
  mut.ct.assocs[estimate > 0 & qval >= .05,fill_class:="high"]
  mut.ct.assocs[estimate > 0 & qval < .05,fill_class:="sig_high"]
  
  #re-order
  
  mut.ct.assocs[,ct_fac:=factor(inhibitor, levels=sub("\\-", ".", ct.ord), 
                                labels=ct.ord, ordered=T)]
  
  top.genes <- mut.ct.assocs[qval < .05 & estimate > 0,cbind(.SD[order(pval)][1:min(5, .N)], ind=1:min(5, .N)), by=ct_fac]
  top.genes <- top.genes[(ct_fac %in% c("cDC-like") & ind %in% 3:5)==F]
  top.genes <- top.genes[(ct_fac %in% c("Promono-like") & ind %in% 2:5)==F]
  top.genes <- top.genes[(ct_fac %in% c("Monocyte-like") & ind %in% 4:5)==F]
  
  mut.ct.assocs <- merge(mut.ct.assocs, top.genes[,.(ct_fac, gene, to_lab=T)], by=c("ct_fac", "gene"), all.x=T, all.y=F)
  
  mut.ct.assocs[, labs:=ifelse(is.na(to_lab)==F, sub("[\\._]", "-", sub("mut\\.", "", gene)), "")]
  
  pos <- position_jitter(height=0, width=.1, seed=123)
  
  mct.plot <- ggplot(data=mut.ct.assocs, mapping=aes(x=ct_fac, y=plot_vals, fill=as.character(ct_fac), label=labs)) + 
    scale_fill_manual(values=ct.colors, guide="none") +
    geom_point(position=pos, size=3, shape=21, alpha=.75) +
    geom_text_repel(position=pos, min.segment.length = 0, seed = 6828, box.padding=.25, family="Arial", size=8/ggplot2:::.pt) +
    theme_bw() + xlab("") + ylab("signed -log10(Pvalue)") +
    theme(
      text=element_text(size=8, family="Arial"),
      axis.text=element_text(size=8, family="Arial"),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
    )
  
  ggsave(mct.plot, file="figures/mut_ct_assoc_v1.pdf", width=114/25.4, height=4)
  
  return("figures/mut_ct_assoc_v1.pdf")
  
}

drug.ct.heatmap <- function(all.cors, inh.inters, ct.ord){
  
  sig.cors <- all.cors[fdr < .05]
  inhib.map <- sig.cors[,.(max_cor=max(abs(estimate))),by=.(inhibitor)]
  
  ct.mat <- reshape2::acast(vg_type ~inhibitor, value.var="estimate", data=all.cors)
  ct.sig.mat <- ct.mat[ct.ord, inhib.map$inhibitor]
  
  inhib.map[,sig_inter:=inhibitor %in% inh.inters[qval < .05]$inhibitor]
  inhib.map[,sug_inter:=inhibitor %in% inh.inters[qval < .1]$inhibitor]
  
  ct.sig.mat <- ct.sig.mat[,inhib.map$inhibitor]
  
  hb <- HeatmapAnnotation(
    inter_sug=anno_points(rep(.8, ncol(ct.sig.mat)), pch=ifelse(inhib.map$sug_inter , 8, NA_integer_),gp=gpar(border = "grey", cex=3), ylim=c(0,1), border=F, axis=F, height=unit(6, "points")),
    inter_sig=anno_points(rep(.8, ncol(ct.sig.mat)), pch=ifelse(inhib.map$sig_inter , 8, NA_integer_),gp=gpar(border = "grey", cex=3), ylim=c(0,1), border=F, axis=F, height=unit(6, "points")),
    show_annotation_name=F, gap=0)
  
  ht <- Heatmap(ct.sig.mat, cluster_rows=F, cluster_columns=T, clustering_method_columns="average",col=circlize::colorRamp2(breaks=c(-.4, 0, .4), colors=c("red", "white","blue")),
                heatmap_legend_param = list(title_position = "leftcenter", direction = "horizontal", title_gp = gpar(fontsize = 8, fontface = "bold", fontfamily="Arial"), labels_gp = gpar(fontsize = 8, fontfamily="Arial")), 
                rect_gp = gpar(col = "grey", lwd=.5), border=T, show_column_names=T, column_title=NULL, column_split=5, name="Correlation", 
                column_names_gp=gpar(fontfamily="Arial", fontsize=6), row_names_gp=gpar(fontfamily="Arial", fontsize=6),
                bottom_annotation=hb)
  
  pdf(file="figures/fig3a.pdf", width=174*0.0393701, height=3)
  #bottom, left, top and right
  draw( ht, heatmap_legend_side = "bottom", padding=unit(c(5.5, 1, 5.5, 1), "points"))
  dev.off()
  
  return("figures/fig3a.pdf")
}

drug.mut.ct.inter.plot <- function(inh.inters, ct.ord, ct.colors){
  
  inh.inters[,plot_vals:=sign(estimate) * -log10(pval)]
  inh.inters[,ct_fac:=factor(ct, levels=sub("\\-",".", ct.ord), labels=ct.ord, ordered=T)]
  
  #note these threshold are for visualization purposes
  inh.inters[,c_levs:="none"]
  inh.inters[qval < .07,c_levs:=gene]
  
  sig.vals <- inh.inters[,min(pval),by=c_levs][order(V1)]
  
  sig.shapes <- setNames(c( 24, 25, 7, 3, 22, 8, 14,21)[seq_along(sig.vals$c_levs)], sig.vals$c_levs)
  
  inh.inters[,inh:=ifelse(qval < .07, inhibitor, "")]
  
  pos <- position_jitter(width = 0.3, height=0, seed = 2)
  
  sug.thresh <- inh.inters[qval < .1, min(-log10(pval))]
  sig.thresh <- inh.inters[qval < .05, min(-log10(pval))]
  
  inter.plot <- ggplot(data=inh.inters, mapping=aes(x=ct_fac, y=plot_vals, shape=c_levs, fill=ct_fac, alpha=qval < .05, size=qval < .05)) + 
    geom_hline(yintercept=c(sug.thresh, sig.thresh, -sug.thresh, -sig.thresh), linetype=c("dashed", "solid", "dashed", "solid")) +
    geom_point(position = pos) +
    scale_fill_manual(values=ct.colors, guide="none") +
    scale_shape_manual(values=sig.shapes, name="Sig. Interaction") +
    scale_alpha_manual(values=c(`TRUE`=1, `FALSE`=.5), guide="none") + 
    scale_size_manual(values=c(`TRUE`=2.5, `FALSE`=1), guide="none") + 
    geom_text_repel(position=pos, mapping=aes(label=inh), color="black", size=8/ggplot2:::.pt, family="Arial", min.segment.length = 0, box.padding = 0.5, point.padding=.25) +
    theme_bw() + xlab("") + ylab("Signed -log10(Pvalue)") + ylim(c(-6, 12)) +
    theme(axis.text = element_text(family="Arial", size=8),
          text=element_text(family="Arial", size=8),
          legend.position = "bottom",
          legend.text=element_text(family="Arial", size=8))
  
  ggsave(inter.plot, file="figures/inhib_ct_inters_v1.pdf", width=6.85, height=4.5)
  
  return("figures/inhib_ct_inters_v1.pdf")
  
}

drug.family.plot <- function(fam.enrich, family.list){
  
  fam.syns <- family.list$synonyms
  
  ssg.dt <- merge(fam.enrich, fam.syns, by="family")
  
  #exclude those with fewer than x patients
  
  ssg.sum <- ssg.dt[,.(med=median(ssgsea),.N),by=synonym]
  
  ssg.sum <- ssg.sum[N > 50]
  
  ssg.sum <- ssg.sum[order(med)]
  
  ssg.sum[,fam_ord:=factor(synonym, levels=synonym, ordered=T)]
  
  ssg.dt <- merge(ssg.dt, ssg.sum, by="synonym")
  
  ridge.plot <- ggplot(data=ssg.dt, mapping=aes(y=fam_ord, x=ssgsea, fill=med)) + 
    ggridges::geom_density_ridges() + 
    scale_fill_gradient2(low="red", mid="white", high="blue", guide="none") +
    geom_vline(xintercept=0, linetype="dashed") +
    #geom_text(data=ssg.sum, mapping=aes(x=40, label=paste("(n=", N, ")"))) +
    ylab("") + xlab("ssGSEA") +
    theme_bw() +
    theme(axis.text = element_text(family="Arial", size=6),
          text=element_text(family="Arial", size=6))
  
  ggsave(ridge.plot, file="figures/ssgsea_ridge_v1.pdf", width=85/25.4, height=4)
  
  return("figures/ssgsea_ridge_v1.pdf")
}

family.mut.assoc.plot <- function(comb.assoc){
  
  comb.assoc[,sig:="Neither"]
  comb.assoc[qval < .05 & glass_d > 0,sig:="Up"]
  comb.assoc[qval < .05 & glass_d < 0,sig:="Down"]
  
  high.assoc <- data.table(gene=c("mut.CBFB.MYH11", "mut.NRAS", "mut.KRAS", "mut.ASXL1", "mut.RUNX1","mut.STAG2","mut.ZRSR2", "mut.TP53","mut.TP53", "mut.IDH2", "mut.NRAS", "mut.FLT3_ITD"), 
                           inhibitor=c("AGC", "CDK", "CDK", "JAK", "PIK", "PIK","PIK","PIK","PIKK", "SRC", "STE7", "Type III RTK PDGFR"))
  
  comb.assoc <- merge(comb.assoc, high.assoc[,.(gene, inhibitor, highlight=T)], all.x=T, all.y=F)
  comb.assoc[is.na(highlight),highlight:=F]
  
  stopifnot(comb.assoc[highlight == T, max(qval) < .05] )
  
  c.volc <- ggplot(data=comb.assoc, mapping=aes(x=glass_d, y=-log10(pval), fill=sig)) + 
    geom_point(size=3, shape=21, mapping=aes(alpha=qval < .05)) +
    geom_vline(xintercept=0, linetype="dashed") +
    geom_text_repel(data=comb.assoc[highlight==T], mapping=aes(label=paste(inhibitor, "+", sub("[\\._]", "-", sub("mut\\.", "", gene)))), 
                    direction="y", min.segment.length = 0,box.padding=.5, xlim=c(-Inf, -2), family="Arial", size=8/ggplot2:::.pt) +
    scale_fill_manual(values=c(`Up`="blue", `Neither`="black", `Down`="red"), guide="none") +
    scale_alpha_manual(values=c(`TRUE`=.75, `FALSE`=.25), guide="none") +
    theme_bw() + ylab("-log10(Pvalue)") + xlab("Effect Size") + xlim(c(-5, 2)) +
    theme(axis.text = element_text(family="Arial", size=8),
          text=element_text(family="Arial", size=8))
  
  ggsave(c.volc, file="figures/beataml_waves1to4_family_assoc_v1.pdf", width=85/25.4, height=85/25.4)
  
  return("figures/beataml_waves1to4_family_assoc_v1.pdf")
  
}

family.ct.heatmap <- function(all.cors, family.list, family.inters, ct.ord){
  
  #add in synonyms
  
  fam.syns <- family.list$synonyms
  
  all.cors <- merge(all.cors, fam.syns, by="family")
  
  family.inters <- merge(family.inters, fam.syns, by.x="inhibitor", by.y="family")
  
  ct.mat <- reshape2::acast(vg_type ~synonym, value.var="estimate", data=all.cors)
  sig.mat <- reshape2::acast(vg_type ~synonym, value.var="estimate", data=all.cors[fdr < .05])
  
  ct.sig.mat <- ct.mat[ct.ord, colnames(sig.mat)]
  rownames(ct.sig.mat) <- ct.ord
  
  hb <- HeatmapAnnotation(sug_inter=anno_points(rep(.8, ncol(ct.sig.mat)), pch=ifelse(colnames(ct.sig.mat) %in% family.inters[qval < .1,synonym], 8, NA_integer_), gp=gpar(border = "grey", cex=3), ylim=c(0,1), border=F, axis=F, height=unit(6, "points")), 
                          sig_inter=anno_points(rep(.8, ncol(ct.sig.mat)), pch=ifelse(colnames(ct.sig.mat) %in% family.inters[qval < .05,synonym], 8, NA_integer_),gp=gpar(border = "grey", cex=3), ylim=c(0,1), border=F, axis=F, height=unit(6, "points")),
                          show_annotation_name=F, gap=0)
  
  ht <- Heatmap(ct.sig.mat, cluster_rows=F, cluster_columns=T, clustering_method_columns="average", col=circlize::colorRamp2(breaks=c(-.4, 0, .4), colors=c("red", "white","blue")),
                heatmap_legend_param = list(title_position = "leftcenter", direction = "horizontal",
                                            title_gp = gpar(fontsize = 8, fontface = "bold", fontfamily="Arial"), labels_gp = gpar(fontsize = 8, fontfamily="Arial")), 
                rect_gp = gpar(col = "grey"), border=T, show_column_names=T, column_title=NULL, column_split=5, name="Correlation", 
                bottom_annotation=hb, column_names_gp=gpar(fontfamily="Arial", fontsize=8), row_names_gp=gpar(fontfamily="Arial", fontsize=8))
  
  pdf(file="figures/vg_family_mat_v1.pdf", width=85/25.4, height=3.5)
  
  draw( ht, heatmap_legend_side = "bottom")
  
  dev.off()
  
  return("figures/vg_family_mat_v1.pdf")
  
}

family.mut.ct.inter.plot <- function(family.inters,family.list, ct.ord, ct.colors){
  
  family.inters[,plot_vals:=sign(estimate) * -log10(pval)]
  family.inters[,ct_fac:=factor(ct, levels=sub("\\-",".", ct.ord), labels=ct.ord, ordered=T)]
  
  #note these threshold are for visualization purposes
  family.inters[,c_levs:="none"]
  family.inters[qval < .06,c_levs:=gene]
  
  sig.vals <- c(unique(family.inters[qval < .06, c_levs]), "none")
  
  sig.shapes <- setNames(c(24, 14, 18, 3, 22, 25,  21)[seq_along(sig.vals)], sig.vals)
  
  fam.syns <- family.list$synonyms
  
  family.inters <- merge(family.inters, fam.syns, by.x="inhibitor", by.y="family")
  
  family.inters[,inh:=ifelse(qval < .06, synonym, "")]
  
  pos <- position_jitter(width = 0.3, height=0, seed = 2)
  
  sug.thresh <- family.inters[qval < .1, min(-log10(pval))]
  sig.thresh <- family.inters[qval < .05, min(-log10(pval))]
  
  inter.plot <- ggplot(data=family.inters, mapping=aes(x=ct_fac, y=plot_vals, shape=c_levs, fill=ct_fac, alpha=qval < .05, size=qval < .05)) + 
    geom_hline(yintercept=c(sug.thresh, sig.thresh, -sug.thresh, -sig.thresh), linetype=c("dashed", "solid", "dashed", "solid")) +
    geom_point(position = pos) +
    scale_fill_manual(values=ct.colors, guide="none") +
    scale_shape_manual(values=sig.shapes, name="Sig. Interaction") +
    scale_alpha_manual(values=c(`TRUE`=1, `FALSE`=.5), guide="none") + 
    scale_size_manual(values=c(`TRUE`=2.5, `FALSE`=1), guide="none") + 
    geom_text_repel(position=pos, mapping=aes(label=inh), color="black", min.segment.length = 0, 
                    box.padding = 0.5, point.padding=.5, family="Arial", size=8/ggplot2:::.pt) +
    theme_bw() + xlab("") + ylab("Signed -log10(Pvalue)") + ylim(c(-6, 11)) +
    guides(shape = guide_legend(nrow = 2, byrow = F)) +
    theme(axis.text = element_text(family="Arial", size=8),
          text=element_text(family="Arial", size=8),
          legend.text=element_text(family="Arial", size=6),
          legend.position="bottom") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  
  ggsave(inter.plot, file="figures/family_ct_inters_v1.pdf", width=84/25.4, height=4)
  
  return("figures/family_ct_inters_v1.pdf")
  
}

stree.plot <- function(s.tree){
  
  pdf(file="figures/rna_mut_ctree_results_v1.pdf", width=12, height=5)
  plot(s.tree$tree, type="simple")
  dev.off()
  
  return("figures/rna_mut_ctree_results_v1.pdf")
  
}


stree.summary.plot <- function(s.tree){
  
  tree.df <- s.tree$data
  tree.df$node <- paste0("Node", predict(s.tree$tree, newdata=tree.df, type="node"))
  
  mod.fit <- survfit(Surv(time=overallSurvival, event=isDead, type="right") ~  node, data = tree.df)
  
  y.plot <- ggsurvplot(mod.fit, data =tree.df, pval=F, ggtheme=theme_bw(), palette=wesanderson::wes_palette(name="Darjeeling1",n=6, type="continuous"))
  
  #get palette mapping
  g <- ggplot_build(y.plot$plot)
  stratas <- unique(g$plot$data$strata)
  
  color.dt <- data.table(colors = g$plot$scales$scales[[2]]$map(stratas), label = as.character(stratas)) 
  
  ggsave(y.plot$plot, file="figures/rna_mut_ctree_curves.pdf", width=4, height=3.5)
  
  #temporary version just to double-check legend
  
  ggsave(y.plot$plot, file="figures/rna_mut_ctree_curves_tmp.pdf", width=8, height=3.5)
  
  return(color.dt)
  
}

m3.pear1.hsc.plot <- function(cur.map, comb.exprs, vg.scores, wgcna.mes){
  
  #relationship between M3 and PEAR1 wrt HSC-like
  
  hsc.dt <- vg.scores$summary[vg_type == "HSC-like",.(ptid, HSC.like=PC1)]
  
  hsc.dt <- cbind(hsc.dt, comb.exprs[hsc.dt$ptid,"PEAR1",drop=F])
  
  hsc.dt <- merge(hsc.dt, wgcna.mes[module == "M3",.(ptid, M3=PC1)], by="ptid")
  
  m3.map <- cur.map[cur_labels == "M3"]
  m3.map[,label:=ifelse(display_label == "PEAR1", "PEAR1", "")]
  
  pos <- position_jitter(width = 0.1, seed = 123)
  
  m3.only <- ggplot(data=m3.map, mapping=aes(x="Mod3 (red)", y=comb_kme, label=label)) +
    geom_boxplot(outlier.shape=NA) +
    geom_point(position=pos, shape=21, alpha=.8, size=3, fill="red") + 
    geom_text_repel(position = pos, box.padding=.5, min.segment.length=0, family="Arial", size=8/ggplot2:::.pt) +
    theme_bw() + xlab("") + ylab("kME") +
    theme(axis.text = element_text(family="Arial", size=8),
          text=element_text(family="Arial", size=8))
  
  m3.pear <- ggplot(data=hsc.dt, mapping=aes(x=PEAR1, y=M3, fill=`HSC.like`)) + 
    geom_point(shape=21, size=3, alpha=.8) +
    viridis::scale_fill_viridis(option="plasma", end=.9, name="HSC-like") +
    theme_bw() + ylab("Mod3 (red)") + 
    theme(axis.text = element_text(family="Arial", size=8),
          text=element_text(family="Arial", size=8),
          legend.text=element_text(family="Arial", size=8))
  
  ggsave(m3.only + m3.pear, file="figures/pear1vm3_v1.pdf", width=174/24.5, height=4)
  
  return("figures/pear1vm3_v1.pdf")
}


clinical.vs.features.plot <- function(clin, feats, wgcna.maps){
  
  rna.clin <- clin[manuscript_rnaseq == "yes" & diseaseStageAtSpecimenCollection == "Initial Diagnosis"]
  
  #fix character encoding
  suppressWarnings(rna.clin[,perc_blasts:=as.numeric(sub(">", "", `%.Blasts.in.BM`, fixed=T))])
  
  mlt.type <- suppressWarnings(melt(rna.clin, measure.vars=c("isDenovo", "isTransformed", "isTherapy"), id.vars=c("ptid"), variable.factor=F))
  mlt.sum <- mlt.type[value %in% c("TRUE"),.(type=paste(unique(variable), collapse=";")), by=.(ptid)]
  #only a couple
  mlt.sum <- mlt.sum[grepl(";", type)==F]
  mlt.sum[,type:=gsub("is", "", type)]
  
  rna.clin <- merge(rna.clin, mlt.sum, by="ptid", all.x=T, all.y=F)
  
  
  rna.feats <- feats[,-grep("mut\\.", names(feats)), with=F]
  rna.feats <- rna.feats[complete.cases(rna.feats)]
  melt.feats <- melt(rna.feats, id.vars=c("ptid"), variable.factor=F)
  
  melt.feats <- merge(melt.feats, wgcna.maps$mod.map, by.x="variable", by.y="cur_labels", all.x=T, all.y=F)
  
  melt.feats[,variable:=ifelse(is.na(prev_modules)==F, paste0(sub("M", "Mod", variable), " (", prev_modules, ")"), sub("\\.", "-", variable))]
  
  rna.clin.f <- merge(rna.clin, melt.feats, by="ptid", allow.cartesian=T)
  
  #categorical list
  
  cat.list <- list(
    Therapy_vs_Denovo=list(column="type", levs=c("Denovo", "Therapy")),
    Transformed_vs_Denovo=list(column="type", levs=c("Denovo", "Transformed")),
    ELN2017_Adverse_vs_Favorable=list(column="ELN2017", levs=c("Favorable", "Adverse")),
    Female_vs_Male=list(column="consensus_sex", levs=c("Male", "Female")),
    BM_vs_PB=list(column="specimenType", levs=c("Peripheral Blood", "Bone Marrow Aspirate"))
  )
  
  cat.res <- rbindlist(lapply(cat.list, function(i){
    
    rbindlist(lapply(split(rna.clin.f, by="variable"), function(j){
      
      tmp.x <- as.data.frame(j[,c(i$column, "value"), with=F])
      tmp.x$class_fac <- factor(tmp.x[[1]], levels=i$levs)
      
      tmp.x <- tmp.x[complete.cases(tmp.x),]
      
      tmp.fit <- glm(class_fac~value, data=tmp.x, family=binomial(link = "logit"))
      
      as.data.table(tidy(tmp.fit))[term=="value"]
      
    }), idcol="variable")
    
  }), idcol="comparison")
  
  cat.res$type <- "Categorical"
  
  #survival
  
  surv.res <- rbindlist(lapply(split(rna.clin.f, by="variable"), function(j){
    
    tmp.x <- as.data.frame(j[is.na(overallSurvival)==F & is.na(isDead)==F])
    
    ph.test <- coxph(Surv(time=overallSurvival, event=isDead, type="right") ~ value, data=tmp.x)
    
    as.data.table(tidy(ph.test))
    
  }),idcol="variable")
  
  surv.res[,`:=`(comparison="All", type="Survival")]
  
  #survival by age
  
  surv.age.res <- rbindlist(lapply(split(rna.clin.f[is.na(ageAtDiagnosis)==F], by=c("variable", "age_cut")), function(j){
    
    tmp.x <- as.data.frame(j[is.na(overallSurvival)==F & is.na(isDead)==F])
    
    ph.test <- coxph(Surv(time=overallSurvival, event=isDead, type="right") ~ value, data=tmp.x)
    
    as.data.table(tidy(ph.test))
    
  }),idcol="var_age")
  
  surv.age.res[,c("variable", "comparison"):=tstrsplit(var_age, "\\.")]
  
  surv.age.res[,`:=`(type="Survival")]
  
  surv.age.res[,comparison:=.capwords(comparison)]
  
  #for continuous
  
  cont.list <- list(
    Age="ageAtDiagnosis",
    Percent_Monocytes="%.Monocytes.in.PB",
    Percent_Neutrophils="%.Neutrophils.in.PB",
    WBC_Count="wbcCount",
    Percent_Blasts=c("perc_blasts")
  )
  
  cont.res <- rbindlist(lapply(cont.list, function(i){
    
    rbindlist(lapply(split(rna.clin.f, by="variable"), function(j){
      
      tmp.x <- as.data.frame(j[,c(i, "value"), with=F])
      names(tmp.x) <- c("cont", "value")
      tmp.x <- tmp.x[complete.cases(tmp.x),]
      
      tmp.fit <- lm(cont~value, data=tmp.x)
      
      as.data.table(tidy(tmp.fit))[term=="value"]
      
    }), idcol="variable")
    
  }), idcol="comparison")
  
  cont.res$type <- "Continuous"
  
  all.res <- rbind(
    cat.res,
    surv.res[,names(cat.res),with=F],
    surv.age.res[,names(cat.res),with=F],
    cont.res[,names(cat.res),with=F]
  )
  
  all.mat <- reshape2::acast(comparison~variable, value.var="statistic",data=all.res)
  
  comp.cat <- unique(all.res[,.(comparison, type)])
  comp.cat[,type_fac:=factor(type, levels=c("Survival", "Continuous", "Categorical"), ordered=T)]
  
  all.mat <- all.mat[comp.cat$comparison,]
  
  v1 <- viridis(n=3)
  
  ht <- Heatmap(all.mat,  rect_gp = gpar(col = "white"), split=comp.cat$type_fac, col=circlize::colorRamp2(breaks=c(-6, 0, 6), colors=v1),
                cluster_columns=T, cluster_rows=T, clustering_method_columns="average", clustering_method_rows="average",
                clustering_distance_columns="euclidean", clustering_distance_rows="euclidean",
                cluster_row_slices = F, column_split=5, name="T or Z statistic",
                column_names_gp=gpar(fontfamily="Arial", fontsize=8), 
                row_names_gp=gpar(fontfamily="Arial", fontsize=8),
                heatmap_legend_param = list(title_position = "leftcenter", direction = "horizontal",
                                            title_gp = gpar(fontsize = 8, fontface = "bold", fontfamily="Arial"), labels_gp = gpar(fontsize = 8, fontfamily="Arial")))
  
  pdf(file="figures/clin_feature_heatmap_v1.pdf", width=174/24.5, height=174/24.5)
  
  draw( ht, heatmap_legend_side = "bottom")
  
  dev.off()
  
  return("figures/clin_feature_heatmap_v1.pdf")
  
}


pear1.eln.plot <- function(p1.rna){
  
  p1.rna[,age_cat:=as.character(age_cut)]
  p1.rna[age_cut %in% c("older","oldest"), age_cat:="older+oldest"]
  
  eln.res <- p1.rna[is.na(age_cut)==F & is.na(ELN2017)==F & ELN2017 != "Intermediate",cbind(as.data.table(tidy(t.test(PEAR1~ELN2017))),n=.N),by=age_cat][,.(age_cat, estimate, p.value, n)]
  
  eln.pear <- ggplot(data=p1.rna[is.na(age_cut)==F & is.na(ELN2017)==F], mapping=aes(x=factor(ELN2017, levels=c("Favorable", "Intermediate","Adverse")), y=PEAR1, fill=ELN2017)) + 
    geom_boxplot(outlier.shape=NA) +
    geom_jitter(height=0, width=.1) +
    scale_fill_discrete(guide="none") +
    facet_grid(.~factor(age_cat, levels=c("young", "middle", "older+oldest"))) +
    annotate("segment", x=1, xend=3, y=8, yend=8) +
    annotate("segment", x=1, xend=1, y=7.75, yend=8.25) +
    annotate("segment", x=3, xend=3, y=7.75, yend=8.25) +
    geom_text(data=eln.res, mapping=aes(x=2, y=8.75, label=paste("Pvalue: ", signif(p.value, digits=3), ifelse(p.value >= .05, " (N.S.)", ""))), inherit.aes=F) +
    theme_bw() + xlab("") +
    theme(axis.text = element_text(family="Arial", size=8),
          text=element_text(family="Arial", size=8),
          strip.text=element_text(family="Arial", size=8))
  
  
  ggsave(eln.pear, file="figures/pear1_v_eln2017.pdf", width=174/24.5, height=3)
  
  return("figures/pear1_v_eln2017.pdf")
}

pear1vslsc17 <- function(p1.rna, lsc.dt){
  
  p1.rna <- merge(p1.rna, lsc.dt, by="ptid")
  
  rna.clin.df <- as.data.frame(p1.rna[age_cut == "young"])
  
  lsc.plot <- .surv.plot.from.tree(rna.clin.df,var.name="LSC17", subset.descr="")
  
  lsc.plot.p <- lsc.plot$plot + theme(axis.text = element_text(family="Arial", size=8),
                                      text=element_text(family="Arial", size=8),
                                      strip.text=element_text(family="Arial", size=8))
  
  p1.plot <- .surv.plot.from.tree(rna.clin.df,var.name="PEAR1", subset.descr="")
  
  p1.plot.p <-  p1.plot$plot + theme(axis.text = element_text(family="Arial", size=8),
                                     text=element_text(family="Arial", size=8),
                                     strip.text=element_text(family="Arial", size=8))
  
  ggsave(p1.plot.p + lsc.plot.p, file="figures/baml_pear1_young_v1.pdf", width=114/24.5, height=3)
  
  return("figures/baml_pear1_young_v1.pdf")
}

pear1vslsc17.hr <- function(p1.rna, lsc.dt){
  
  all.dt <- merge(p1.rna, lsc.dt, by="ptid")
  
  ctree.res.dt <- rbindlist(list(
    #all patients
    PEAR1.all.full=.hazard.ratio.cat(all.dt, split.col="PEAR1", type="tree"),
    LSC17.all.full=.hazard.ratio.cat(all.dt, split.col="LSC17", type="tree"),
    #young patients
    PEAR1.all.young=.hazard.ratio.cat(all.dt[age_cut == "young",], split.col="PEAR1", type="tree"),
    LSC17.all.young=.hazard.ratio.cat(all.dt[age_cut == "young",], split.col="LSC17", type="tree"),
    #BM
    PEAR1.BM.full=.hazard.ratio.cat(all.dt[specimenType == "Bone Marrow Aspirate",], split.col="PEAR1", type="tree"),
    LSC17.BM.full=.hazard.ratio.cat(all.dt[specimenType == "Bone Marrow Aspirate",], split.col="LSC17", type="tree"),
    #young BM
    PEAR1.BM.young=.hazard.ratio.cat(all.dt[age_cut == "young" & specimenType == "Bone Marrow Aspirate",], split.col="PEAR1", type="tree"),
    LSC17.BM.young=.hazard.ratio.cat(all.dt[age_cut == "young" & specimenType == "Bone Marrow Aspirate",], split.col="LSC17", type="tree")
    
  ), idcol="version")
  
  
  ctree.res.dt[,c("pred", "tissue", "cohort"):=tstrsplit(version, "\\.")]
  
  #median
  
  med.res.dt <- rbindlist(list(
    #all patients
    PEAR1.all.full=.hazard.ratio.cat(all.dt, split.col="PEAR1", type="median"),
    LSC17.all.full=.hazard.ratio.cat(all.dt, split.col="LSC17", type="median"),
    #young patients
    PEAR1.all.young=.hazard.ratio.cat(all.dt[age_cut == "young",], split.col="PEAR1", type="median"),
    LSC17.all.young=.hazard.ratio.cat(all.dt[age_cut == "young",], split.col="LSC17", type="median"),
    #BM
    PEAR1.BM.full=.hazard.ratio.cat(all.dt[specimenType == "Bone Marrow Aspirate",], split.col="PEAR1", type="median"),
    LSC17.BM.full=.hazard.ratio.cat(all.dt[specimenType == "Bone Marrow Aspirate",], split.col="LSC17", type="median"),
    #young BM
    PEAR1.BM.young=.hazard.ratio.cat(all.dt[age_cut == "young" & specimenType == "Bone Marrow Aspirate",], split.col="PEAR1", type="median"),
    LSC17.BM.young=.hazard.ratio.cat(all.dt[age_cut == "young" & specimenType == "Bone Marrow Aspirate",], split.col="LSC17", type="median")
    
  ), idcol="version")
  
  
  med.res.dt[,c("pred", "tissue", "cohort"):=tstrsplit(version, "\\.")]
  
  comb.res.dt <- rbind(cbind(ctree.res.dt, type="ctree"), cbind(med.res.dt[,names(ctree.res.dt),with=F], type="median"))
  comb.res.dt[,worked:=T]
  
  comb.res.dt[is.na(exp.coef.), `:=`(exp.coef.=2.5, lower..95=1, upper..95=1, worked=F)]
  
  p3 <- ggplot(data=comb.res.dt, mapping=aes(y=version, x=exp.coef., xmin=lower..95, xmax=upper..95, color=type, alpha=worked)) + geom_pointrange(position = position_dodge(width=.75)) +
    geom_vline(xintercept=1, linetype="dashed") +
    geom_text(data=comb.res.dt[worked==F], mapping=aes(label="N/A"), alpha=1, nudge_y=-.2, show.legend = F) +
    scale_alpha_manual(values=c(`TRUE`=1, `FALSE`=0), guide="none") + 
    scale_color_manual(values=c(ctree="black", median="darkgrey"), name="Split Type") +
    facet_grid(tissue+cohort~., scales="free_y", labeller=function(x) lapply(x, .capwords)) + 
    scale_y_discrete(labels=function(x) sapply(strsplit(x, "\\."), "[[", 1)) + 
    theme_bw() + ylab("") + xlab("Hazard Ratio (+/- 95% CI)")
  
  ggsave(p3, file="figures/pear1_vs_lsc17_hr.pdf", width=4.4, height=4)
  
  return("figures/pear1_vs_lsc17_hr.pdf") 
}

exprs.boxplots.by.cohort <- function(clin, comb.exprs, wv.cols){
  
  rna.clin <- clin[manuscript_rnaseq == "yes"]
  
  mlt.exprs <- data.table(reshape2::melt(comb.exprs, as.is=T))
  names(mlt.exprs) <- c("ptid", "gene", "norm")
  
  mlt.exprs.w.coh <- merge(mlt.exprs, rna.clin[,.(ptid, cohort)], by="ptid")
  
  mlt.exprs.w.coh[,pt_fac:=factor(ptid, levels=rna.clin[order(cohort)]$ptid, ordered=T)]
  
  mlt.exprs.w.coh[,cohort:=sub("Waves", "Waves ", cohort)]
  
  bplot <- ggplot(data=mlt.exprs.w.coh, mapping=aes(x=pt_fac, y=norm, fill=cohort)) + geom_boxplot(outlier.shape=NA) +
    scale_fill_manual(values=wv.cols, name="Cohort") + 
    theme_bw() + theme(axis.text.x=element_blank()) + xlab("Patient Samples") + ylab("Normalized Expression") +
    theme(axis.text = element_text(family="Arial", size=8),
          text=element_text(family="Arial", size=8),
          legend.text=element_text(family="Arial", size=8),
          legend.position="bottom")
  
  ggsave(bplot, file="figures/exprs_box_by_cohort_v1.pdf", width=30, height=5)
  
  "figures/exprs_box_by_cohort_v1.pdf"
}

module.waves.overlays <- function(comb.mes, wgcna.maps, wv.cols){
  
  #for this just plot the 14 beataml paper modules (best.match correspondence)
  
  comb.mes <- merge(comb.mes, wgcna.maps$mod.map[,.(module=cur_labels, prev_modules)] , by="module")
  
  comb.mes[,mod_name:=paste0(sub("M", "Mod", module), " (", prev_modules, ")")]
  
  comb.mes <- comb.mes[module != "M0"]
  
  comb.plot <- ggplot(data=comb.mes, mapping=aes(x=PC1, y=PC2, fill=cohort, color=module)) + geom_point(shape=21, alpha=.75) + 
    facet_wrap(~mod_name, ncol=7, scales="free") + 
    scale_fill_manual(values=wv.cols, name="") + 
    scale_color_manual(values=wgcna.maps$mod.map[,setNames(prev_modules, cur_labels)], guide="none") +
    theme_bw() + xlab("PC1") + ylab("PC2") +
    theme(
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      strip.text.x = element_blank() , 
      strip.background = element_blank(),
      plot.margin = unit( c(0,0,0,0) , units = "lines" ),
      text=element_text(size=8, family="Arial"),
      legend.text=element_text(size=8, family="Arial"))
  
  ggsave(comb.plot, file="figures/previous_train_test_overlay_mes.pdf", width=174/24.5, height=3.5)
  
  return("figures/previous_train_test_overlay_mes.pdf")
}

plot.wgcna.mods.by.drug2 <- function(cast.ck){
  
  drug.me.plot <- ggplot(data=cast.ck, mapping=aes(x=`Waves1+2`, y=`Waves3+4`, fill=for_col)) + geom_point(size=3, shape=21) +
    geom_hline(yintercept=0, linetype="dashed") + 
    geom_vline(xintercept=0, linetype="dashed") +
    scale_fill_gradient2(low="red", mid="white", high="blue", guide="none") +
    geom_label_repel(data=cast.ck[order(for_col)][1:3],mapping=aes(label=paste(inhibitor, titles, sep=" + ")), min.segment.length=0, box.padding=.5, fill="white", family="Arial", size=8/ggplot2:::.pt) +
    geom_label_repel(data=cast.ck[order(-for_col)][1:3],mapping=aes(label=paste(inhibitor, titles, sep=" + ")), min.segment.length=0, box.padding=.5, fill="white", family="Arial", size=8/ggplot2:::.pt) +
    theme_bw() + xlim(c(-1, 1)) + ylim(c(-1, 1)) +
    theme(axis.text = element_text(family="Arial", size=8),
          text=element_text(family="Arial", size=8),
          strip.text=element_text(family="Arial", size=8))
  
  ggsave(drug.me.plot, file="figures/drug_cor_kme_v2.pdf", width=114/25.4, height=114/25.4)
  
  return("figures/drug_cor_kme_v2.pdf")
  
}

vg.heatmap <- function(vg.scores){
  
  base.col <- circlize::colorRamp2(seq(-6, 6, length = 50), plasma(50, end=.9))
  
  colfun <- circlize::colorRamp2(seq(-2, 2, length = 50), plasma(50))
  
  #Note only for Waves 1/2
  
  tmp.pc <- vg.scores$wv12_results[["Monocyte-like"]]$summary[order(PC1)]
  
  sub.ha <- HeatmapAnnotation(PC1 = anno_barplot(tmp.pc$PC1, gp=gpar(col=base.col(tmp.pc$PC1), fill=base.col(tmp.pc$PC1)), axis_param=list(gp=gpar(fontsize=8, fontfamily="Arial"))),
                              show_legend = c(F, T, F), height=unit(1, "in"))
  
  gene.kme <- vg.scores$wv12_results[["Monocyte-like"]]$kme[order(-kme)]
  
  ht <- Heatmap(vg.scores$wv12_sexprs[gene.kme$display_label,tmp.pc$ptid], col=colfun, top_annotation=sub.ha,
                show_column_names=F, cluster_columns=F, cluster_rows=T, show_row_names=F, row_names_side="left",
                row_title_rot=0, show_row_dend=F,
                heatmap_legend_param = list(color_bar = "continuous", legend_direction = "horizontal",
                                            legend_width = unit(5, "cm"), title_position = "leftcenter",
                                            title_gp = gpar(fontsize = 8, fontface = "bold", fontfamily="Arial"), 
                                            labels_gp = gpar(fontsize = 8, fontfamily="Arial")), 
                name="Scaled Exprs", border=T, column_title = "Monocyte-like")
  
  pdf(file="figures/van_galen_clusters_pcs_ml_v1.pdf", width=114/25.4, height=3)
  
  draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legend = T)#, padding=padding)
  
  dev.off()
  
  return("figures/van_galen_clusters_pcs_ml_v1.pdf")
}

vg.by.fab.plot <- function(clin, vg.scores){
  
  rna.clin <- clin[manuscript_rnaseq=="yes"]
  
  with.clin <- merge(rna.clin[,c("ptid", "fabBlastMorphology")], vg.scores$summary, by="ptid")
  
  #roll the a's and b's up to include with M4, M5, etc. exclude  M0/M1's
  
  with.clin[,clean_fab:=gsub(" ", "", sub("[ab]", "", fabBlastMorphology))]
  with.clin[clean_fab == "M4eo", clean_fab:="M4"]
  
  with.clin <- with.clin[clean_fab %in% paste0("M", 0:7) ]
  
  with.clin[,fab_fac:=factor(clean_fab, levels=paste0("M", 0:7), ordered=T)]
  
  with.clin[,vg_fac:=factor(vg_type, levels=c("HSC-like","Progenitor-like", "GMP-like", "Promono-like", "Monocyte-like", "cDC-like" ), ordered=T)]
  
  set.seed(123)
  
  vg.fab <- ggplot(data=with.clin[vg_type %in% c("Monocyte-like",  "HSC-like")], mapping=aes(x=fab_fac, y=PC1, fill=PC1)) + 
    geom_boxplot(outlier.shape=NA) +
    geom_jitter(height=0, width=.1, shape=21, size=3) +
    geom_hline(yintercept=0, linetype="dashed") +
    facet_grid(vg_fac~.) + 
    viridis::scale_fill_viridis(option="plasma", end=.9, guide="none") +
    theme_bw() + xlab("FAB") + ylab("van Galen Signature") +
    theme(
      text=element_text(size=8, family="Arial"),
      strip.text=element_text(size=8, family="Arial"),
      axis.text=element_text(size=8, family="Arial")
    )
  
  ggsave(vg.fab, file="figures/vg_fab_v1.pdf", width=114/25.4, height=4)
  
  return("figures/vg_fab_v1.pdf")
}

ven.pano.plot <- function(inhib, features){
  
  mlt.feats <- melt(features, measure.vars=c("Monocyte.like"), id.vars=c("ptid", "mut.FLT3_ITD"), variable.factor=F)
  
  inhib.m <- merge(inhib[inhibitor %in% c("Venetoclax", "Panobinostat")], mlt.feats[is.na(mut.FLT3_ITD ) == F & is.na(value) == F], by="ptid")
  
  inhib.m[,variable:=sub("[\\._]", "-", variable)]
  
  vp.examp <- ggplot(data=inhib.m[inhibitor %in% c("Venetoclax", "Panobinostat") & variable == "Monocyte-like"], mapping=aes(y=auc, x=value, fill=value)) +
    geom_point(shape=21, size=3, alpha=.5) +
    viridis::scale_fill_viridis(option="plasma", end=.9, guide="none") +
    geom_smooth(method="lm", formula=y~x, se=F) +
    facet_wrap(~inhibitor) +
    theme_bw() + xlab("Monocyte-like") + ylab("AUC") +
    theme(axis.text = element_text(family="Arial", size=8),
          text=element_text(family="Arial", size=8),
          strip.text=element_text(family="Arial", size=8))
  
  ggsave(vp.examp, file="figures/ven_pano_examp_v1.pdf", width=3.34, height=3.34)
  
  return("figures/ven_pano_examp_v1.pdf")
  
}

sora.inter.inter.plot <- function(inhib, features, inh.inters){
  
  mlt.feats <- melt(features, measure.vars=names(features)[grepl("like", names(features))], id.vars=c("ptid", "mut.FLT3_ITD"), variable.factor=F)
  
  inhib.m <- merge(inhib[inhibitor %in% c("Sorafenib")], mlt.feats[is.na(mut.FLT3_ITD ) == F & is.na(value) == F], by="ptid")
  
  inhib.m[,variable:=sub("[\\._]", "-", variable)]
  
  sig.sor <- inh.inters[inhibitor == "Sorafenib" & gene == "mut.FLT3_ITD"][order(pval)]
  sig.sor[,dots:=""]
  sig.sor[qval < .1, dots:="*"]
  sig.sor[qval < .05, dots:=paste0(dots, "*")]
  sig.sor[,mut.FLT3_ITD:=1]
  sig.sor[,ct_fac:=factor(sub("[\\._]", "-", ct), levels=sub("[\\._]", "-", ct), ordered=T)]
  
  inhib.m[,ct_fac:=factor(variable, levels=sig.sor[,sub("[\\._]", "-", ct)], ordered=T)]
  
  s.examp <- ggplot(data=inhib.m, mapping=aes(y=auc, x=value, fill=value)) +
    geom_point(shape=21, size=3, alpha=.5) + 
    geom_smooth(method="lm", formula=y~x, se=F) +
    geom_text(data=sig.sor, mapping=aes(x=8, y=10, label=dots), size=10,  inherit.aes=F) +
    viridis::scale_fill_viridis(option="plasma", end=.9, guide="none") +
    facet_grid(ct_fac~ifelse(mut.FLT3_ITD==1, "FLT3-ITD Positive", "FLT3-ITD Negative")) +
    theme_bw()  + xlab("Cell-type Score") + ylab("Sorafenib AUC") +
    theme(axis.text = element_text(family="Arial", size=8),
          text=element_text(family="Arial", size=8),
          strip.text=element_text(family="Arial", size=8))
  
  ggsave(s.examp, file="figures/sora_inter_examp_v1.pdf", width=6.85, height=6)
  
}

pear.tissue.surv.plot <- function(p1.dt){
  
  p1.df <- as.data.frame(p1.dt)
  
  #BM
  
  bm.surv <- .surv.plot.from.tree(p1.df[p1.df$tissue == "BM",], var.name="PEAR1", subset.descr="Bone Marrow Samples")
  
  #PB
  
  pb.surv <- .surv.plot.from.tree(p1.df[p1.df$tissue == "PB/Leuk",], var.name="PEAR1", subset.descr="PB/Leuk. Samples")
  
  #now in the young subset
  
  young.df <- p1.df[is.na(p1.df$age_cut) == F & p1.df$age_cut == "young",]
  
  #BM Young
  
  ybm.surv <- .surv.plot.from.tree(young.df[young.df$tissue == "BM",], var.name="PEAR1", subset.descr="Bone Marrow Young Patients")
  
  #PB Young
  
  ypb.surv <- .surv.plot.from.tree(young.df[young.df$tissue == "PB/Leuk",], var.name="PEAR1", subset.descr="PB/Leuk. Young Patients")
  
  ggsave(
    pb.surv$plot + 
      theme(
        text=element_text(family="Arial", size=8),
        axis.text = element_text(family="Arial", size=8)
      ) +
      bm.surv$plot + 
      theme(
        text=element_text(family="Arial", size=8),
        axis.text = element_text(family="Arial", size=8)
      ) +
      ypb.surv$plot + 
      theme(
        text=element_text(family="Arial", size=8),
        axis.text = element_text(family="Arial", size=8)
      ) +
      ybm.surv$plot + 
      theme(
        text=element_text(family="Arial", size=8),
        axis.text = element_text(family="Arial", size=8)
      ) +
      plot_layout(ncol = 2), file="figures/pear1_tissue_surv.pdf", width=8, height=174/24.5)
  
  return("figures/pear1_tissue_surv.pdf")
  
}

pear1.vs.lsc17.both.median.split <- function(p1.dt, lsc17.dt){
  
  p1.lsc <- merge(p1.dt[diseaseStageAtSpecimenCollection == "Initial Diagnosis" & is.na(overallSurvival)==F & is.na(isDead)==F], lsc17.dt, by="ptid")
  
  p1.lsc[,pear1_cat:=ifelse(PEAR1 > median(PEAR1), "high", "low")]
  p1.lsc[,lsc17_cat:=ifelse(LSC17 > median(LSC17), "high", "low")]
  
  p1.lsc.df <- as.data.frame(p1.lsc)
  
  lsc.surv <- .surv.plot.from.tree(p1.lsc.df, var.name="LSC17", subset.descr="All Initial Diagnosis")
  
  p1.surv <- .surv.plot.from.tree(p1.lsc.df, var.name="PEAR1", subset.descr="All Initial Diagnosis")
  
  #now split using median same as above
  
  p.m.fit <- survfit(Surv(time=overallSurvival, event=isDead, type="right") ~  pear1_cat, data = p1.lsc.df) 
  p.m.surv <- ggsurvplot(p.m.fit, data =  p1.lsc.df, pval=T, ggtheme=theme_bw()) + ggtitle("PEAR1 (All Initial Diagnosis; median)")
  
  lsc.m.fit <- survfit(Surv(time=overallSurvival, event=isDead, type="right") ~  lsc17_cat, data = p1.lsc.df) 
  lsc.m.surv <- ggsurvplot(lsc.m.fit, data =  p1.lsc.df, pval=T, ggtheme=theme_bw()) + ggtitle("LSC17 (All Initial Diagnosis; median)")
  
  ggsave(p1.surv$plot + 
           theme(
             text=element_text(family="Arial", size=8),
             axis.text = element_text(family="Arial", size=8)
           ) +
           lsc.surv$plot + 
           theme(
             text=element_text(family="Arial", size=8),
             axis.text = element_text(family="Arial", size=8)
           ) +
           p.m.surv$plot + 
           theme(
             text=element_text(family="Arial", size=8),
             axis.text = element_text(family="Arial", size=8)
           ) +
           lsc.m.surv$plot + 
           theme(
             text=element_text(family="Arial", size=8),
             axis.text = element_text(family="Arial", size=8)
           ) +
           plot_layout(ncol=2), file="figures/pear1_lsc17_both_split_median.pdf", width=8, height=174/24.5)
  
  return("figures/pear1_lsc17_both_split_median.pdf")
  
}

pear1.assocs.plot <- function(p1.dt, features){
  
  mut.mat<- as.matrix(features[,grepl("mut\\.", names(features)),with=F])
  rownames(mut.mat) <- features$ptid
  mut.mat <- mut.mat[complete.cases(mut.mat),]
  
  p1.mat <- as.matrix(p1.dt[,.(PEAR1)])
  rownames(p1.mat) <- p1.dt$ptid
  
  mut.p1.assocs <- .run.assoc.welch(mut.mat, p1.mat, min.muts=5)
  
  mut.p1.assocs[,fdr:=p.adjust(pval, method="BY")]
  
  mut.p1.assocs[,highlight:= gene %in% paste0("mut.", c("ASXL1", "RUNX1", "TP53", "SRSF2", "GATA2.MECOM", "NPM1", "CEBPA"))]
  
  stopifnot(mut.p1.assocs[highlight==T,max(fdr)] < .05)
  
  plot2 <- ggplot(data=mut.p1.assocs, mapping=aes(x=glass_d, y=-log10(pval),  fill=fdr < .05)) + geom_point(shape=21, mapping=aes(size=num_muts)) +
    geom_text_repel(data=mut.p1.assocs[fdr < .05], mapping=aes(label=sub("\\.", "-", sub("mut\\.", "", gene))), min.segment.length=0, box.padding=.5, family="Arial", size=8/ggplot2:::.pt) +
    scale_size_continuous(name="Number Mutations") + 
    geom_vline(xintercept=0, linetype="dashed") +
    geom_hline(yintercept=-log10(.05), linetype="dashed") +
    theme_bw() + xlab("Effect Size") +
    theme(axis.text = element_text(family="Arial", size=8),
          text=element_text(family="Arial", size=8),
          legend.text=element_text(family="Arial", size=8))
  
  ggsave(plot2, file="figures/pear1_mut_assocs_v1.pdf", width=174/24.5, height=4)
  
  return("figures/pear1_mut_assocs_v1.pdf")
}

pear1.aml.type.plot <- function(clin, p1.dt, vg.scores, wv.cols){
  
  rna.clin <- clin[manuscript_rnaseq=="yes" & diseaseStageAtSpecimenCollection == "Initial Diagnosis"]
  
  rcp <- merge(p1.dt[,.(ptid, PEAR1)], rna.clin, by="ptid")
  
  mlt.type <- suppressWarnings(melt(rcp, measure.vars=c("isDenovo", "isTransformed", "isTherapy"), id.vars=c("ptid", "cohort" ,"PEAR1"), variable.factor=F))
  mlt.sum <- mlt.type[value %in% c("TRUE"),.(type=paste(unique(variable), collapse=";")), by=.(ptid, cohort, PEAR1)]
  #only a couple
  mlt.sum <- mlt.sum[grepl(";", type)==F]
  
  devtra <- rbindlist(
    list(
      dtrans=data.table(tidy(t.test(PEAR1~type,data=mlt.sum[type %in% c("isDenovo", "isTransformed")]))),
      dther=data.table(tidy(t.test(PEAR1~type,data=mlt.sum[type %in% c("isDenovo", "isTherapy")])))
    ), idcol="type"
  )
  
  type.meds <- mlt.sum[,.(meds=median(PEAR1)),by=type][order(meds)]
  
  mlt.sum[,type_fac:=factor(type, levels=type.meds$type, labels=sub("is", "", type.meds$type), ordered=T)]
  
  plot1 <- ggplot(data=mlt.sum, mapping=aes(x=type_fac, y=PEAR1)) + geom_boxplot(outlier.shape=NA) +
    geom_jitter(height=0, width=.1, shape=21, size=2, mapping=aes(fill=sub("Waves", "Waves ", cohort)), alpha=.85) +
    scale_fill_manual(values=wv.cols, name="Cohort") +
    annotate("segment", x=1, xend=3, y=8.75, yend=8.75) +
    annotate("segment", x=1, xend=1, y=8.5, yend=9) +
    annotate("segment", x=3, xend=3, y=8.5, yend=9) +
    annotate("text", x=2, y=9.25, label=paste("Pvalue <", signif(devtra[type == "dtrans"]$p.value, 3)), family="Arial", size=8/ggplot2:::.pt) +
    annotate("segment", x=1, xend=2, y=7.75, yend=7.75) +
    annotate("segment", x=1, xend=1, y=7.5, yend=8) +
    annotate("segment", x=2, xend=2, y=7.5, yend=8) +
    annotate("text", x=1.5, y=8.25, label=paste("Pvalue :", signif(devtra[type == "dther"]$p.value, 3)), family="Arial", size=8/ggplot2:::.pt) +
    theme_bw() + xlab("") +
    theme(axis.text = element_text(family="Arial", size=8),
          text=element_text(family="Arial", size=8),
          strip.text=element_text(family="Arial", size=8),
          legend.text=element_text(family="Arial", size=8),
          legend.position = "bottom")
  
  vg.sum <- merge(vg.scores$summary[,.(ptid, vg_type, cohort, PC1)], mlt.sum[type %in% c("isDenovo", "isTransformed")], by=c("ptid", "cohort"))
  
  vg.sum[,cats:=ifelse(type=="isDenovo","Denovo", "Transformed")]
  
  hsc.sum <- vg.sum[vg_type %in% c("HSC-like", "Monocyte-like")]
  
  lplot <- ggplot(data=hsc.sum, mapping=aes(x=PC1,y=PEAR1, color=cats, fill=cats )) + 
    geom_point(size=2, shape=21, alpha=.25) +
    geom_smooth(method="lm", se=F, formula=y~x) + 
    facet_wrap(~vg_type) +
    scale_fill_discrete(name="") +
    scale_color_discrete(name="") +
    theme_bw() + xlab("Cell-Type Scores") +
    theme(axis.text = element_text(family="Arial", size=8),
          text=element_text(family="Arial", size=8),
          strip.text=element_text(family="Arial", size=8),
          legend.text=element_text(family="Arial", size=8),
          legend.position = "bottom")
  
  ggsave(plot1 + lplot + plot_layout(widths=c(1.25, 2)), file="figures/pear1_vs_aml_type.pdf", width=174/24.5, height=4)
  
  return("figures/pear1_vs_aml_type.pdf")
}


