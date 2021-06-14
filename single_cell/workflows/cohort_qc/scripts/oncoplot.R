#!/usr/bin/env Rscript

library(maftools)


oncoplot = function(read_maf, oncoplot_path, genes){
    png(filename=oncoplot_path, units="px", width=1600, height=1600, res=300)
    
    maftools::oncoplot(maf=read_maf,showTumorSampleBarcodes=TRUE,genes=genes)
    dev.off()
}



main = function(){
    args = commandArgs(trailingOnly=TRUE)
    genes=c("PPM1D", "TP53",  "BRCA1", "BRCA2", "MECOM", "RB1", "PTEN", "PALB2","ERBB2", "CDK12", "PIK3CA", "KRAS", "CCNE1", "MYC")

    maf_file = args[1]
    vcNames=args[2]
    cn=args[3]
    oncoplot_path = args[4]


    vcNames=read.table(vcNames,header=TRUE)$Variant_Classification

    maf = maftools::read.maf(maf=maf_file, cnTable=cn, vc_nonSyn=vcNames)
    
    oncoplot(maf, oncoplot_path, genes)


}


main()


