library(targets)
source("R/preprocess.R")
source("R/figures.R")
source("R/summarization.R")
source("R/analysis_helpers.R")
tar_option_set(packages=c("extrafont","data.table", "openxlsx", "ggplot2", "ggrepel", "patchwork", "ComplexHeatmap", "survival", "survminer", "partykit", "broom", "qvalue", "viridis", "RColorBrewer", "GSVA"))

list(
  
  #files
  
  tar_target(
    clin_wkst,
    "data/beataml_wv1to4_clinical.xlsx",
    format = "file"
  ),
  
  tar_target(
    inhibitor_file,
    "data/beataml_probit_curve_fits_v4_dbgap.txt",
    format = "file"
  ),
  
  tar_target(
    exprs_file,
    "data/beataml_waves1to4_norm_exp_dbgap.txt",
    format = "file"
  ),
  
  tar_target(
    inhibitor_family_file,
    "data/beataml_drug_families.xlsx",
    format = "file"
  ),
  
  tar_target(
    muts_file,
    "data/beataml_wes_wv1to4_mutations_dbgap.txt",
    format = "file"
  ),
  
  tar_target(
    vg_supp_file,
    "data/1-s2.0-S0092867419300947-mmc3.xlsx",
    format = "file"
  ),
  
  tar_target(
    wgcna_maps,
    local({mget(load("data/wgcna/merged_older_wgcna_kme.RData"))})
  ),
  
  #settings
  
  tar_target(
    wv_colors,
    c(`Waves 1+2`="#023FA5", `Waves 3+4`="#7D87B9")
  ),
  
  tar_target(
    ct_order,
    c("HSC-like", "Progenitor-like", "GMP-like","Promono-like" , "Monocyte-like", "cDC-like")
  ),
  
  tar_target(
    ct_colors,
    setNames(as.character(wesanderson::wes_palette("IsleofDogs1", n=6, type="continuous")), ct_order)
  ),
  
  #pre-processing
  
  tar_target(
    clin_data,
    clinical.for.genomic.analysis(clin_wkst)
  ),
  
  tar_target(
    clin_sum,
    clinical.for.summary(clin_data)
  ),
  
  tar_target(
    comb_exprs,
    get.combined.exprs(exprs_file, clin_data)
  ),
  
  tar_target(
    mut_list,
    mutation.only.dataset(muts_file, clin_data)
  ),
  
  tar_target(
    fus_mat,
    fusion.only.dataset(clin_data)
  ),
  
  tar_target(
    inhib_data,
    inhibitor.dataset(inhibitor_file, clin_data)
  ),
  
  tar_target(
    vg_genes,
    parse.van.galen.supp(vg_supp_file)
  ),
  
  tar_target(
    family_list,
    get.family.list(inhibitor_family_file)
  ),
  
  tar_target(
    wgcna_mes,
    predict.modules(comb_exprs, clin_data, wgcna_maps)
  ),
  
  tar_target(
    vg_scores,
    van.galen.scores(comb_exprs, clin_data, vg_genes)
  ),
  
  tar_target(
    feature_data,
    mutation.expression.features(vg_scores, wgcna_mes, mut_list, fus_mat, clin_data)
  ),
  
  tar_target(
    family_enrich,
    drug.family.enrichment(family_list, inhib_data)
  ),
  
  tar_target(
    p1_rna,
    pear1.rna.clinical(comb_exprs, clin_data)
  ),
  
  tar_target(
    p1_denovo_surv,
    pear1.rna.denovo.survival(p1_rna)
  ),
  
  tar_target(
    lsc17,
    compute.lsc17(comb_exprs)
  ),
  
  #summarization
  
  tar_target(
    mut_freqs,
    mutation.freq.by.cohort(mut_list, fus_mat)
  ),
  
  tar_target(
    inhib_sum,
    summarize.denovo.aucs(inhib_data, clin_data)
  ),
  
  tar_target(
    vg_mut_assocs,
    vg.muts.assocs(feature_data)
  ),
  
  tar_target(
    inhib_cors,
    drug.ct.cors(inhib_data, vg_scores),
  ),
  
  tar_target(
    inhib_inters,
    inhibitor.mutation.interactions(inhib_data, feature_data)
  ),
  
  tar_target(
    family_mut_assocs,
    family.mut.assocs(family_enrich,family_list, feature_data)
  ),
  
  tar_target(
    family_cors,
    family.ct.cors(family_enrich, vg_scores)
  ),
  
  tar_target(
    family_inters,
    family.mutation.interactions(family_enrich, feature_data)
  ),
  
  tar_target(
    ctree_results,
    survival.trees(feature_data, clin_data)
  ),
  
  tar_target(
    combined_kme_map,
    combined.kmes(comb_exprs,wgcna_maps, wgcna_mes)
  ),
  
  #main figures
  
  tar_target(
    figure1a,
    mut.freq.by.wave.plot(mut_freqs),
    format = "file"
  ),
  
  tar_target(
    figure1b,
    surv.by.wave.plot(clin_sum),
    format = "file"
  ),
  
  tar_target(
    figure1c,
    denovo.inhib.plot(inhib_sum),
    format = "file"
  ),
  
  tar_target(
    figure1d,
    train.test.mut.assocs(inhib_data, mut_list, wv_colors),
    format = "file"
  ),
  
  tar_target(
   figure2a,
   vg.exprs.plot2(wgcna_mes, wgcna_maps, vg_scores, ct_order),
   format = "file"
  ),
  
  tar_target(
    figure2b,
    vg.muts.assoc.plot(vg_mut_assocs, ct_order, ct_colors),
    format = "file"
  ),
  
  tar_target(
    tableS2,
    write.vma(vg_mut_assocs),
    format = "file"
  ),
  
  tar_target(
    figure3a,
    drug.ct.heatmap(inhib_cors, inhib_inters, ct_order),
    format = "file"
  ),
  
  tar_target(
    figure3b,
    drug.mut.ct.inter.plot(inhib_inters, ct_order, ct_colors),
    format = "file"
  ),
  
  tar_target(
    tableS3,
    write.cor.ints(inhib_inters, inhib_cors),
    format = "file"
  ),
  
  tar_target(
    figure4a,
    drug.family.plot(family_enrich, family_list),
    format = "file"
  ),
  
  tar_target(
    figure4b,
    family.mut.assoc.plot(family_mut_assocs),
    format = "file"
  ),
  
  tar_target(
    tableS5,
    write.fma(family_mut_assocs, family_list),
    format = "file"
  ),
  
  tar_target(
    figure5a,
    family.ct.heatmap(family_cors, family_list, family_inters, ct_order),
    format = "file"
  ),
  
  tar_target(
    figure5b,
    family.mut.ct.inter.plot(family_inters, family_list, ct_order, ct_colors),
    format = "file"
  ),
  
  tar_target(
    tableS6,
    write.fam.cor.inter(family_inters, family_cors, family_list),
    format = "file"
  ),
  
  tar_target(
    figure6a,
    clinical.vs.features.plot(clin_data, feature_data, wgcna_maps),
    format = "file"
  ),
  
  tar_target(
    figure6b,
    m3.pear1.hsc.plot(combined_kme_map, comb_exprs, vg_scores, wgcna_mes),
    format = "file"
  ),
  
  tar_target(
    figure6c1,
    stree.plot(ctree_results),
    format = "file"
  ),
  
  tar_target(
    colors_figure6c,
    stree.summary.plot(ctree_results),
  ),
  
  tar_target(
    figure7a,
    pear1.eln.plot(p1_rna),
    format = "file"
  ),
  
  tar_target(
    figure7b,
    pear1vslsc17(p1_denovo_surv,lsc17),
    format = "file"
  ),
  
  tar_target (
   figure7c,
   pear1vslsc17.hr(p1_denovo_surv, lsc17),
   format = "file"
  ),
  
  tar_target(
    tableS7,
    write.mod.membership(exprs_file, wgcna_maps, combined_kme_map),
    format = "file"
  ),
  
  #supplemental figures/analyses
  
  tar_target(
    figureS1a,
    exprs.boxplots.by.cohort(clin_data, comb_exprs, wv_colors),
    format = "file"
  ),
  
  tar_target(
    figureS1b,
    module.waves.overlays(wgcna_mes, wgcna_maps, wv_colors),
    format = "file"
  ),
  
  tar_target(
    figureS1c,
    wgcna.mods.by.drug2(clin_data, inhib_data, wgcna_mes, wgcna_maps),
    format = "file"
  ),
  
  tar_target(
    figureS2a,
    vg.heatmap(vg_scores),
    format = "file"
  ),
  
  tar_target(
    figureS2b,
    vg.by.fab.plot(clin_data, vg_scores),
    format = "file"
  ),
  
  tar_target(
    figureS3a,
    ven.pano.plot(inhib_data, feature_data),
    format = "file"
  ),
  
  tar_target(
    figureS3b,
    sora.inter.inter.plot(inhib_data, feature_data, inhib_inters),
    format = "file"
  ),
  
  tar_target(
    figureS4abcd,
    pear.tissue.surv.plot(p1_denovo_surv),
    format = "file"
  ),
  
  tar_target(
    figureS5ab,
    pear1.aml.type.plot(clin_data, p1_rna, vg_scores, wv_colors),
    format = "file"
  ),
  
  tar_target(
    figureS5c,
    pear1.assocs.plot(p1_rna, feature_data),
    format = "file"
  ),
  
  tar_target(
    figureS6abcd,
    pear1.vs.lsc17.both.median.split(p1_rna, lsc17),
    format = "file"
  )
  
)
  