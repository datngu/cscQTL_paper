setwd("/Users/datn/github/circQTL_analysis")

require(data.table)
library(ggpubr)
require(ggplot2)
library(scales)
library("ggplotify")





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


get_gen_exp_df <- function(chr, gene, snp, dir){
	# chr = "6"
	# gene = "6__90206569__90271941__ENSG00000112182"
	# snp = "rs60066732"
	# dir = "consensus_3/recount_qtl_input"
	
	vcffn = paste0(dir, "/", chr, ".vcf.gz")
	vcf = fread(vcffn)
	mysnp = vcf[vcf$ID == snp,]
	ref = paste0(mysnp$REF, mysnp$REF)
	het = paste0(mysnp$REF, mysnp$ALT)
	alt = paste0(mysnp$ALT, mysnp$ALT)
	mysnp[mysnp == "0/0"] = ref
	mysnp[mysnp == "1/0" | mysnp == "0/1"] = het
	mysnp[mysnp == "1/1"] = alt
	mysnp2 = as.character(mysnp[1,-c(1:9)])

	exp_fn = paste0(dir, "/", chr, ".bed.gz")
	exp = fread(exp_fn)

	myexp = exp[exp$TargetID == gene,]
	myexp2 = as.numeric(myexp[1,-c(1:4)])

	df = data.frame(Genotype = mysnp2, Expression = myexp2)
	df$Genotype = factor(mysnp2, levels = c(ref, het, alt))
	return(df)
	#ggplot(df, aes(x = Genotype, y = Expression)) + geom_boxplot(fill = "grey80", color = "black") + geom_jitter(width = 0.2, height = 0, color = "blue") + labs(title = "My Box Plot with Data Points", x = "Genotype", y = "Expression") + theme_classic()
}


plot_gen_exp <- function(chr, gene, snp, dir, title = NULL){
	# chr = "6"
	# gene = "6__90206569__90271941__ENSG00000112182"
	# snp = "rs60066732"
	# dir = "consensus_3/recount_qtl_input"
	
	vcffn = paste0(dir, "/", chr, ".vcf.gz")
	vcf = fread(vcffn)
	mysnp = vcf[vcf$ID == snp,]
	ref = paste0(mysnp$REF, mysnp$REF)
	het = paste0(mysnp$REF, mysnp$ALT)
	alt = paste0(mysnp$ALT, mysnp$ALT)
	mysnp[mysnp == "0/0"] = ref
	mysnp[mysnp == "1/0" | mysnp == "0/1"] = het
	mysnp[mysnp == "1/1"] = alt
	mysnp2 = as.character(mysnp[1,-c(1:9)])

	exp_fn = paste0(dir, "/", chr, ".bed.gz")
	exp = fread(exp_fn)

	myexp = exp[exp$TargetID == gene,]
	myexp2 = as.numeric(myexp[1,-c(1:4)])

	df = data.frame(Genotype = mysnp2, Expression = myexp2)
	df$Genotype = factor(mysnp2, levels = c(ref, het, alt))
	#return(df)
	p = ggplot(df, aes(x = Genotype, y = Expression, fill = Genotype,)) + geom_boxplot(color = "black", alpha = 0.9, width = 0.5) + labs(title = title, x = paste0("Genotype - ", snp), y = "Normalized expression") + theme_classic()
	p = p + geom_jitter(width = 0.3, height = 0) 
	p = p + theme(legend.position = c(0.8, 0.9))
	#p = p + geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(1))
	return(p)
}


################
t1d_recount= plot_all("consensus_3/coloc_recount/coloc_result_T1D.tsv.gz.Rdata")
cd_recount= plot_all("consensus_3/coloc_recount/coloc_result_CD.tsv.gz.Rdata")
ibd_recount= plot_all("consensus_3/coloc_recount/coloc_result_IBD.tsv.gz.Rdata")
uc_recount= plot_all("consensus_3/coloc_recount/coloc_result_UC.tsv.gz.Rdata")





####### Fig 4. T1D-circBACH2

load("consensus_3/coloc_recount/coloc_result_T1D.tsv.gz.Rdata")

###### adding effect size and p-values

crna = 	qtl_hit
pick = crna$V7 == "rs60066732" & crna$V1 == "6__90206569__90271941__ENSG00000112182"
slope = round(crna[pick,]$V14, 3)
pval = format(crna[pick,]$V6, scientific = T)
lab_b = paste0("estimated slope=", slope, "\np-value=", pval)

mrna = read_normial_qtl("consensus_3/salmon_qtl_mapping_nominal")
pickc = mrna$V8 == "rs60066732" & mrna$V1 == "ENSG00000112182"
slopec = round(mrna[pickc,]$V14, 3)
pvalc = format(mrna[pickc,]$V12, scientific = T)
lab_c = paste0("estimated slope=", slopec, "\np-value=", pvalc)




A <- as.grob(t1d_recount$`6__90206569__90271941__ENSG00000112182`)
A + annotate(geom="text", x=1, y=1, label= "rs60066732",color="black")

B = plot_gen_exp(chr = 6, gene = "6__90206569__90271941__ENSG00000112182", snp = "rs60066732", dir = "consensus_3/recount_qtl_input", title = "circBACH2")

C = plot_gen_exp(chr = 6, gene = "ENSG00000112182", snp = "rs60066732", dir = "consensus_3/salmon_qtl_input", title = "BACH2 gene")

B = B + annotate(geom="text", x=1, y=-1.2, label= lab_b,color="black")
C = C + annotate(geom="text", x=1, y=-1.5, label= lab_c,color="black")


## merge figures

r1 = ggarrange(A, ncol = 1, nrow = 1) + annotate(geom="text", x=0.8, y=0.55, label= "rs60066732",color="black") + annotate(geom="text", x=0.15, y=0.97, label= "T1D GWAS - circBACH2",color="black")

r2 = ggarrange(B, C, ncol = 2, labels = c("B.", "C."))

m = ggarrange(r1, r2,  heights= c(1, 0.8), nrow = 2, align = "v", labels = "A." )

pdf( file= "output_paper/Fig4.pdf",  width= 10, height= 10)
m
dev.off()


#######

# crna = read_normial_qtl("consensus_3/recount_qtl_mapping_nominal")

# mrna = read_normial_qtl("consensus_3/salmon_qtl_mapping_nominal")




####### Fig S.2,3,4,5

A1 <- as.grob(cd_recount$`6__90206569__90271941__ENSG00000112182`)

A2 <- as.grob(ibd_recount$`6__90206569__90271941__ENSG00000112182`)

A3 <- as.grob(cd_recount$`1__155676548__155679512__ENSG00000163374`)

A4 <- as.grob(ibd_recount$`1__155676548__155679512__ENSG00000163374`)

S2 = ggarrange(A1, ncol = 1, nrow = 1) + annotate(geom="text", x=0.15, y=0.97, label= "CD GWAS - circBACH2",color="black")

S3 = ggarrange(A2, ncol = 1, nrow = 1) + annotate(geom="text", x=0.15, y=0.97, label= "IBD GWAS - circBACH2",color="black")

S4 = ggarrange(A3, ncol = 1, nrow = 1) + annotate(geom="text", x=0.15, y=0.97, label= "CD GWAS - circYY1AP1",color="black")

S5 = ggarrange(A4, ncol = 1, nrow = 1) + annotate(geom="text", x=0.15, y=0.97, label= "IBD GWAS - circYY1AP1",color="black")




pdf( file= "output_paper/S2.pdf",  width= 10, height= 6)
S2
dev.off()


pdf( file= "output_paper/S3.pdf",  width= 10, height= 6)
S3
dev.off()


pdf( file= "output_paper/S4.pdf",  width= 10, height= 6)
S4
dev.off()


pdf( file= "output_paper/S5.pdf",  width= 10, height= 6)
S5
dev.off()



###########


t1d_recount= plot_all("consensus_3/coloc_recount/coloc_result_T1D.tsv.gz.Rdata")
cd_recount= plot_all("consensus_3/coloc_recount/coloc_result_CD.tsv.gz.Rdata")
ibd_recount= plot_all("consensus_3/coloc_recount/coloc_result_IBD.tsv.gz.Rdata")
uc_recount= plot_all("consensus_3/coloc_recount/coloc_result_UC.tsv.gz.Rdata")





####### Fig 4. T1D-circBACH2

#load("consensus_3/coloc_recount/coloc_result_T1D.tsv.gz.Rdata")













