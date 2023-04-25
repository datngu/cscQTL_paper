setwd("/Users/datn/github/circQTL_analysis")

require(data.table)
library(ggpubr)
require(ggplot2)
library(scales)
library("ggplotify")
library(locuscomparer)




my_plotter <- function(df){
  gwas_fn = eqtl_fn = "tem.txt"
  marker_col = "rsID"
  gwas_pvalue = "gwas_pvalue"
  eqtl_pvalue = "eqtl_pvalue"
  fwrite(df, file = gwas_fn, row.names = F, sep = "\t")
  res = locuscompare(in_fn1=gwas_fn, in_fn2=eqtl_fn, title1="GWAS", title2="eQTL", marker_col1= marker_col, pval_col1=gwas_pvalue, marker_col2=marker_col, pval_col2=eqtl_pvalue)
  return(res)
}
#my_plotter(df)



plot_all <- function(data){
  locuscompare_res = list()
  load(data)
  ids = as.character( qtl_hit$V1[qtl_hit$PP.H4.abf >= 0.5])
  if(length(ids> 0)){
    for(id in ids){
      df = res[[id]]
      locuscompare_res[[id]] = my_plotter(df)
    }
  }
  return(locuscompare_res)
}


read_normial_qtl <- function(dir, chr_list = as.character(c(1:22))){
	res = fread(paste0(dir,"/chr_1_nominal_pass_best_peer_factors.tsv.gz"))
	for(i in 2:length(chr_list)){
		chr = chr_list[i]
		tem = fread(paste0(dir,"/chr_", chr ,"_nominal_pass_best_peer_factors.tsv.gz"))
		res = rbind(res, tem)
	}
	return(res)
}










t1d_recount= plot_all("consensus_3/coloc_recount/coloc_result_T1D.tsv.gz.Rdata")
cd_recount= plot_all("consensus_3/coloc_recount/coloc_result_CD.tsv.gz.Rdata")
ibd_recount= plot_all("consensus_3/coloc_recount/coloc_result_IBD.tsv.gz.Rdata")
uc_recount= plot_all("consensus_3/coloc_recount/coloc_result_UC.tsv.gz.Rdata")









####### Fig S.2,3,4,5



A <- as.grob(t1d_recount$`6__90206569__90271941__ENSG00000112182`)

A1 <- as.grob(cd_recount$`6__90206569__90271941__ENSG00000112182`)

A2 <- as.grob(ibd_recount$`6__90206569__90271941__ENSG00000112182`)

A3 <- as.grob(cd_recount$`1__155676548__155679512__ENSG00000163374`)

A4 <- as.grob(ibd_recount$`1__155676548__155679512__ENSG00000163374`)

S2 = ggarrange(A, ncol = 1, nrow = 1) + annotate(geom="text", x=0.8, y=0.55, label= "rs60066732",color="black") + annotate(geom="text", x=0.15, y=0.97, label= "T1D GWAS - circBACH2",color="black")

S3 = ggarrange(A1, ncol = 1, nrow = 1) + annotate(geom="text", x=0.15, y=0.97, label= "CD GWAS - circBACH2",color="black")

S4 = ggarrange(A2, ncol = 1, nrow = 1) + annotate(geom="text", x=0.15, y=0.97, label= "IBD GWAS - circBACH2",color="black")

S5 = ggarrange(A3, ncol = 1, nrow = 1) + annotate(geom="text", x=0.15, y=0.97, label= "CD GWAS - circYY1AP1",color="black")

S6 = ggarrange(A4, ncol = 1, nrow = 1) + annotate(geom="text", x=0.15, y=0.97, label= "IBD GWAS - circYY1AP1",color="black")




pdf( file= "output_paper/S2_new.pdf",  width= 10, height= 6)
S2
dev.off()


pdf( file= "output_paper/S3_new.pdf",  width= 10, height= 6)
S3
dev.off()


pdf( file= "output_paper/S4_new.pdf",  width= 10, height= 6)
S4
dev.off()


pdf( file= "output_paper/S5_new.pdf",  width= 10, height= 6)
S5
dev.off()

pdf( file= "output_paper/S6_new.pdf",  width= 10, height= 6)
S6
dev.off()

