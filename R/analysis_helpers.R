#from ?toupper
.capwords <- function(s, strict = FALSE) {
  cap <- function(s) paste(toupper(substring(s, 1, 1)),
                           {s <- substring(s, 2); if(strict) tolower(s) else s},
                           sep = "", collapse = " " )
  sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}

.run.assoc.welch <- function(mut.mat, inhib.mat, min.muts=5){
  
  mut.mat <- mut.mat[intersect(rownames(mut.mat), rownames(inhib.mat)),]
  
  tmp.inhib <- inhib.mat[rownames(mut.mat),,drop=F]
  
  mut.t <- rbindlist(lapply(setNames(colnames(mut.mat), colnames(mut.mat)), function(gene){
    
    if (sum(mut.mat[,gene]) >= min.muts){
      
      inhib.t <- rbindlist(lapply(setNames(colnames(tmp.inhib), colnames(tmp.inhib)), function(inh){
        
        tmp.inp <- data.frame(inh=tmp.inhib[,inh], gene=1-mut.mat[,gene])
        
        tmp.inp <- tmp.inp[complete.cases(tmp.inp),]
        
        if (length(unique(tmp.inp$gene)) == 2 && sum(1-tmp.inp$gene) >= min.muts){
          
          fit <- t.test(inh~gene, data=tmp.inp, var.equal = FALSE, alternative = "two.sided")
          
          data.table(estimate=diff(rev(fit$estimate)), 
                     num_inhib=nrow(tmp.inp),
                     num_muts=sum(1-tmp.inp$gene),
                     sd_neg=sd(tmp.inp[tmp.inp$gene == 1,"inh"]), 
                     sd_pos=sd(tmp.inp[tmp.inp$gene == 0,"inh"]),
                     t=fit$statistic, se=fit$stderr, pval=fit$p.val)
        }else{
          
          NULL
        }
        
      }), idcol="inhibitor")
      
    }else{
      NULL
      
    }
    
  }), idcol="gene")
  
  #secondary filter
  tmp.mut.t <- mut.t[num_muts >= min.muts & is.na(estimate)==F]
  
  #from a variety of sources including https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6562325/
  #estimate should already be pos - neg so keeping neg as the control:
  tmp.mut.t[,glass_d:=estimate/sd_neg]
  
  #variance also from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6562325/
  tmp.mut.t[,num_neg:=num_inhib - num_muts]
  tmp.mut.t[,glass_var:=(num_inhib/(num_muts*num_neg))+((glass_d^2)/(2*(num_neg-1)))]
  
  tmp.mut.t
}

.run.drug.mut.ct.inters <- function(mut.ct.mat, inhib.mat, mut.re="mut\\.", ct.re="\\.like", min.muts=10){
  
  mut.cols <- colnames(mut.ct.mat)[grep(mut.re, colnames(mut.ct.mat))]
  ct.cols <- colnames(mut.ct.mat)[grep(ct.re, colnames(mut.ct.mat))]
    
  mut.ct.mat <- mut.ct.mat[intersect(rownames(mut.ct.mat), rownames(inhib.mat)),]
  
  tmp.inhib <- t(inhib.mat)[,rownames(mut.ct.mat)]
  
  inter.dt <- rbindlist(lapply(setNames(ct.cols, ct.cols), function(ct){
    
    mut.t <- rbindlist(lapply(setNames(mut.cols, mut.cols), function(gene){
      
      if (sum(mut.ct.mat[,gene]) >= min.muts){
        
        design <- model.matrix(~gene*ct, data=data.frame(gene=mut.ct.mat[,gene], ct=mut.ct.mat[,ct]))
        
        fit <- suppressWarnings(limma::eBayes(limma::lmFit(tmp.inhib[,rownames(design)], design)))
        
        #from eBayes docs--ordinary.t
        tmp.dt <- data.table(inhibitor=rownames(fit), estimate=fit$coef[,"gene:ct"], se=fit$stdev.unscaled[,"gene:ct"] * fit$sigma)
        
        tmp.dt <- cbind(tmp.dt, num_inhib=rowSums(is.na(tmp.inhib)==F), num_muts=rowSums(is.na(tmp.inhib)==F & matrix(mut.ct.mat[,gene] == 1, nrow=nrow(tmp.inhib), ncol=ncol(tmp.inhib), byrow=T)))
        
        tmp.dt
        
      }else{
        NULL
        
      }
      
    }), idcol="gene")
    
  }), idcol="ct")
  
  inter.dt[,t:=estimate/se]
  inter.dt[,pval:=2*pt(-abs(t), num_inhib-2)]
  
  #secondary filter
  inter.dt <- inter.dt[num_muts >= min.muts & is.na(estimate)==F]
  
  inter.dt
}

.hazard.ratio <- function(use.dt, surv.col){
  
  use.df <- as.data.frame(use.dt)
  
  use.df$var <- use.df[[surv.col]]
  
  ph.sum <- summary(coxph(Surv(time=overallSurvival, event=isDead, type="right") ~  var, data=use.df))
  
  ph.res <- as.data.table(ph.sum$conf.int)
  names(ph.res) <- make.names(names(ph.res))
  
  ph.res
}

.hazard.ratio.cat <- function(use.df, split.col, type=c("tree", "median"), minc=.95, minb=20){
  
  type <- match.arg(type)
  
  if (type == "tree"){
    
    use.form <- as.formula(paste("Surv(time=overallSurvival, event=isDead, type='right') ~  ", split.col))
    
    use.tree <- ctree(use.form, 
                      data = use.df, 
                      control=ctree_control(mincriterion = minc, minbucket=minb, maxdepth=1))
    
    if (depth(use.tree) == 0){
      return(data.table(exp.coef.=NA_real_, exp..coef.=NA_real_, lower..95=NA_real_, upper..95=NA_real_, n_high=NA_real_, n_low=NA_real_))
    }
    
    use.df$split <- factor(ifelse(use.df[[split.col]] <=split_node(node_party(use.tree[[1]]))$breaks, "low", "high"), levels=c("low", "high"))
    
    stopifnot(sum(diag(table(use.df$split, predict(use.tree, type="node")))) == nrow(use.df))
    
  }else if (type == "median"){
    
    use.med <- median(use.df[[split.col]], na.rm=T)
    
    use.df$split <- factor(ifelse(use.df[[split.col]] <= use.med, "low", "high"), levels=c("low", "high"))
    
  }
  
  ph.sum <- summary(coxph(Surv(time=overallSurvival, event=isDead, type="right") ~  split, data=use.df))
  
  ph.res <- as.data.table(ph.sum$conf.int)
  names(ph.res) <- make.names(names(ph.res))
  
  ph.res <- cbind(ph.res, n_high=sum(use.df$split == "high"), n_low=sum(use.df$split == "low"))
  
  return(ph.res)
}

.surv.plot.from.tree <- function(use.df, var.name="PEAR1", subset.descr=""){
  
  use.df$split <- paste("All", var.name)
  split.point <- ""
  
  use.form <- as.formula(paste("Surv(time=overallSurvival, event=isDead, type='right') ~  ", var.name))
  
  tree.fit <- ctree(use.form, 
                   data = use.df, 
                   control=ctree_control(mincriterion = 0.95, minbucket=20, maxdepth=1))
  
  if (depth(tree.fit) > 0){
    use.df$split <- ifelse(use.df[[var.name]] <=split_node(node_party(tree.fit[[1]]))$breaks, "low", "high")
    split.point <- paste("; Split:", round(split_node(node_party(tree.fit[[1]]))$breaks, digits=3))
  }
  
  use.fit <- survfit(Surv(time=overallSurvival, event=isDead, type="right") ~  split, data = use.df) 
  
  use.surv <- suppressWarnings(ggsurvplot(use.fit, data = use.df, pval=T, ggtheme=theme_bw(), conf.int = F)) + ggtitle(paste(var.name, "(",paste0(subset.descr,split.point),")"))
  
  use.surv
}
