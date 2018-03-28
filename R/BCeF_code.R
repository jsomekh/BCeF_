
# biocLite("pROC")
#library("pROC")

#' A plot.roc.AUC Function
#'
#' This function allows you to plot the ROC curves and AUC.
#' @param gold.standard.bool is a boolean vector representing the the gold standard.
#' @param predictor.df
#' @keywords fff
#' @export
#' @examples
#' BCeF()
plot.roc.AUC<-function(gold.standard.bool, predictor.df, main.text, legend.to.use, color.to.use)
{
  # browser()
  if(length(gold.standard.bool) != nrow(predictor.df))
  {
    print("error gold standard and predictor has different length")
    return(0)
  }
  cases.gold = length(which(gold.standard.bool==1))
  
  roc_obj_linear = roc(gold.standard.bool, predictor.df[,1])
  l.auc = auc(roc_obj_linear)
  plot(roc_obj_linear, 
       main = paste0(main.text, ",cases==1:",cases.gold),
       col = color.to.use[1],
       cex.main=0.8)
  auc.vec<-l.auc
  final.legend.to.use = paste0(legend.to.use[1], " (auc=",round(l.auc,3),")")
  
  for(i in 2:ncol(predictor.df))
  {
    roc_obj_tmp = roc(gold.standard.bool, predictor.df[,i])
    tmp.auc = auc(roc_obj_tmp)
    lines(roc_obj_tmp, col = color.to.use[i])
    final.legend.to.use = c(final.legend.to.use, paste0(legend.to.use[i], " (auc=",round(tmp.auc,3),")"))
  }
  
  legend("bottomright",cex=1,
         legend = final.legend.to.use, 
         col = color.to.use, 
         lwd = 0.8, 
         lty = c(1,1,1,1))#type of line
  
}
#########################################################################################
# the function generates performance (ROC and AUC) of adjusted data against raw data
# input:
#   input.edata, input.covariates.df, input.edata.description, gold.standard
# giant.all.gold -  a dataframe, first column is gene 1, second column is gene2, third columns is 1 for true associations and 0 for false associations
#                   the gene names or symbols should match between the edata and the gold standard
# input.covariates.df - this can include, known or unknown covariates, each as a vector in the pheno dataframe
# in the ourput plot the raw data is ROC curve is colored in black and the adjusted data roc curve is colored in red
#########################################################################################
#' A BCeF Function
#'
#' This function allows you to express your love of cats.
#' @param input.edata the expression data set to be adjusted.
#' @keywords fff
#' @export
#' @examples
#' BCeF()
BCeF<-function(input.edata, input.covariates.df, input.gold.standard, input.edata.description = "", input.adjustment.method.description = "")
{
  #load the normalized data per tissue
  raw.edata = as.matrix(input.edata) #file.name.to.load(ts.list[i])
  ##adjust with linear regression
  factors             = sapply(1:ncol(input.covariates.df), function(x) paste0("pheno.f[,", x,"]"))
  my.formula          = reformulate(termlabels = factors, response = 'raw.edata')
  lm.fitted.edata     = lm(my.formula)
  adjusted.reads      = t(lm.fitted.edata$residuals)#the residuals aftera djustment
  
  #generate the ROC, AUC comparing raw and adjusted data
  gold.standard.to.use            = input.gold.standard#giant.adipose.net#giant.all.net#giant.adipose.net#giant.all.gold#giant.adipose.gold#giant.adipose.net#giant.all.gold#giant.adipose.gold#giant.all.gold#giant.all.net#giant.adipose.net#rbind(giant.all.CAV2.PLS3.net.cutoff.0.5, giant.all.CAV2.PLS3.net.cutoff.0.05)#giant.adipose.CAV2.PLS3.net.cutoff.0.5, giant.adipose.CAV2.PLS3.net.cutoff.0.05)#giant.adipose.CAV2.PLS3.net.cutoff.0.5#giant.adipose.INSR.net.cutoff.0.5#giant.adipose.CAV2.PLS3.net.cutoff.0.5
  colnames(giant.all.gold)        = c("Gene1","Gene2","Confidence")
  rownames(gold.standard.to.use)  = c(1:nrow(gold.standard.to.use))
  
  bins.cors.df                  = gold.standard.to.use
  bins.cors.df$binAll.raw.r     = NA#add a zero column
  bins.cors.df$binAll.raw.pval  = NA
  bins.cors.df$binAll.r         = NA#add a zero column
  bins.cors.df$binAll.pval      = NA
  
  for(j in 1:nrow(gold.standard.to.use))
  {
    gene1.ensmbl = get.ensmbl(gold.standard.to.use$Gene1[j])[1]
    gene2.ensmbl = get.ensmbl(gold.standard.to.use$Gene2[j])[1]
    
    if( (gene1.ensmbl %in% rownames(raw.edata)) & (gene2.ensmbl %in% rownames(raw.edata)))
    {
      gene1.bin.all.raw.vec   = gene2.bin.all.raw.vec = 0
      gene1.bin.all.raw.vec   = raw.edata[gene1.ensmbl,]
      gene2.bin.all.raw.vec   = raw.edata[gene2.ensmbl,]
      result.all.raw          = cor.test(gene1.bin.all.raw.vec, gene2.bin.all.raw.vec, method = "spearman")
      
      gene1.bin.all.vec       = gene2.bin.all.vec = 0
      gene1.bin.all.vec       = adjusted.reads[gene1.ensmbl,]
      gene2.bin.all.vec       = adjusted.reads[gene2.ensmbl,]
      result.all              = cor.test(gene1.bin.all.vec, gene2.bin.all.vec, method = "spearman")
      
      bins.cors.df[j,"binAll.r"] = result.all$estimate 
      bins.cors.df[j,"binAll.raw.r"] = result.all.raw$estimate
      bins.cors.df[j,"binAll.raw.pval"] = result.all.raw$p.value
      bins.cors.df[j,"binAll.pval"] = result.all$p.value
    }else{
      print(paste0(j," is index of not found gene in edata"))
      j = j+1
    }
  }
  
  bins.cors.df            = bins.cors.df[-c(which(is.na(bins.cors.df$binAll.raw.r))),]
  rownames(bins.cors.df)  = c(1:nrow(bins.cors.df))
  gold.standard.bool      = bins.cors.df$Confidence
  pvals.df.log = 0
  pvals.df.log = -log(pvals.df,10)
  tmp = pvals.df.log
  tmp[is.infinite(tmp)]<-0
  max.val = max(tmp)
  pvals.df.log[is.infinite(pvals.df.log)]<-(max.val+10)#change Inf to max value
  pvals.df.log  = pvals.df.log[,c(input.adjustment.method.description, "Raw")]
  legend.to.use = c(input.adjustment.method.description, "Raw")
  color.to.use  = c("red","black")
  
  plot.roc.AUC(gold.standard.bool, 
                   pvals.df.log, 
                   input.edata.description, 
                   legend.to.use, 
                   color.to.use)
  
}


