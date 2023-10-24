
require(data.table)
library(ggpubr)
require(ggVennDiagram)
require(ggplot2)



read_circ_id <- function(tool, dir = "40_tcells/consensus_3"){
  df = fread(paste0(dir, "/quantile_norm/", tool, "_quantile_norm.tsv"))
  id = df$TargetID
  id = id[grepl("EN", id)]
  idx = sapply(id, FUN = function(x){unlist(gregexpr("__", x))[3]})
  id = substr(id, 1, idx-1)
  return(id)
}


read_circ_id_qtl <- function(tool, dir = "40_tcells/consensus_3"){
  df = fread(paste0(dir, "/qtl_mapping_apply_qvalue/", tool, ".tsv"))
  id = df$V1
  id = id[grepl("EN", id)]
  idx = sapply(id, FUN = function(x){unlist(gregexpr("__", x))[3]})
  id = substr(id, 1, idx-1)
  return(id)
}



### read circIDs of circQTL

circall = unique(read_circ_id_qtl("circall", dir = "40_tcells/consensus_3"))
ciri2 = unique(read_circ_id_qtl("ciri2", dir = "40_tcells/consensus_3"))
circexp2 = unique(read_circ_id_qtl("circexp2", dir = "40_tcells/consensus_3"))

all = c(circall, ciri2, circexp2)
t = table(all)

Supporting_method_1 = names(t)[t>=1]
Supporting_method_2 = names(t)[t>=2]
Supporting_method_3 = names(t)[t>=3]

## cscQTL

recount_1 = unique(read_circ_id_qtl("recount", dir = "40_tcells/consensus_1"))
recount_2 = unique(read_circ_id_qtl("recount", dir = "40_tcells/consensus_2"))
recount_3 = unique(read_circ_id_qtl("recount", dir = "40_tcells/consensus_3"))



get_circQTL_type = function(method_list, circ_list){
  res = c()
  for(m in names(method_list)){
      tem = method_list[[m]]
      res = c(res, sum(tem %in% circ_list))
  }
  names(res) = names(method_list)
  df = data.frame(Methods = names(res), Cumsum = res)
  row.names(df) = NULL
  return(df)
}

get_circQTL_all = function(method_list){
  res = c()
  for(m in names(method_list)){
      tem = method_list[[m]]
      res = c(res, length(tem))
  }
  names(res) = names(method_list)
  df = data.frame(Methods = names(res), Cumsum = res)
  row.names(df) = NULL
  return(df)
}

f1bx = list("cscQTL_1" = recount_1, "cscQTL_2" = recount_2, "cscQTL_3" = recount_3, Circall = circall, CIRI2 = ciri2, CIRCexplorer2 = circexp2)
#f1bx = list("cscQTL_1" = recount_1, "cscQTL_2" = recount_2, "cscQTL_3" = recount_3)

tem3 = get_circQTL_type(f1bx[1:3], Supporting_method_3)
tem3$Type = "Supported by all 3 single methods"
tem3$Count = tem3$Cumsum

tem2 = get_circQTL_type(f1bx[1:3], Supporting_method_2)
tem2$Type = "Supported by 2 single methods"
tem2$Count = tem2$Cumsum - tem3$Cumsum

tem1 = get_circQTL_type(f1bx[1:3], Supporting_method_1)
tem1$Type = "Supported by 1 single methods"
tem1$Count = tem1$Cumsum - tem2$Cumsum


tem_unique = get_circQTL_all(f1bx[1:3])
tem_unique$Type = "cscQTL specific"
tem_unique$Count = tem_unique$Cumsum - tem1$Cumsum


tem_unique2 = get_circQTL_all(f1bx[4:6])
tem_unique2$Type = "Single method"
tem_unique2$Count = tem_unique2$Cumsum


###### bovine_circQTL

read_circ_id_qtl_bovince <- function(tool, dir){
  df = fread(paste0(dir, "/qtl_mapping_apply_qvalue/", tool, ".tsv"))
  id = df$V1
  idx = sapply(id, FUN = function(x){unlist(gregexpr("__", x))[3]})
  id = substr(id, 1, idx-1)
  return(id)
}

bo_ciri2 = read_circ_id_qtl_bovince('ciri2', "40_tcells/bovince_circQTL/qtl_results/ciri2")
bo_circQTL_circexp2 = read_circ_id_qtl_bovince('circexp2', "40_tcells/bovince_circQTL/qtl_results/circexp2")
bo_circRNA_finder = read_circ_id_qtl_bovince('circRNA_finder', "40_tcells/bovince_circQTL/qtl_results/circRNA_finder")


bo_merge = unique(c(bo_ciri2, bo_circQTL_circexp2, bo_circRNA_finder))

bv = list("BOVINCE_circQTL" = bo_merge, "BV_CIRI2" = bo_ciri2, "BV_CIRCexplorer2" = bo_circQTL_circexp2, "BV_circRNA_finder" = bo_circRNA_finder)

tem_bv = get_circQTL_all(bv)
tem_bv$Type = "Single method with BOVINCE_circQTL setting"
tem_bv$Type[1] =  "BOVINCE_circQTL"
tem_bv$Count = tem_bv$Cumsum



dfb = rbind(tem3, tem2, tem1, tem_unique, tem_unique2, tem_bv)

lv = c('cscQTL_1', 'cscQTL_2', 'cscQTL_3','Circall', 'CIRI2', 'CIRCexplorer2', 'BOVINCE_circQTL', 'BV_CIRI2','BV_circRNA_finder', 'BV_CIRCexplorer2')

dfb$Methods = factor(dfb$Methods, levels = lv)
#dfb$Type = factor(dfb$Type, levels = unique(dfb$Type))

A = ggplot(data= dfb, aes(x=Methods, y=Count, fill=Type)) + geom_bar(stat="identity", width=0.5) + geom_text(aes(y = Cumsum, label = Count), vjust=1.6, color="white", size=3.5)

A = A + theme_classic() + theme(legend.position="left", axis.title.x=element_blank()) + theme(legend.position = c(0.8, 0.8)) + theme(legend.text=element_text(size=6)) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

nice_colors <- c("#4169E1", "#DC143C", "#228B22", "#00BFFF", "#DAA520", "#9370DB", "#8FBC8F", "#FF6347")

A = A + scale_fill_manual(values = nice_colors)



## B


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
  res = locuscompare(in_fn1=gwas_fn, in_fn2=eqtl_fn, title1="GWAS", title2="eQTL", marker_col1= marker_col, pval_col1=gwas_pvalue, marker_col2=marker_col, pval_col2=eqtl_pvalue, genome = 'hg38')
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


t1d_recount= plot_all("40_tcells/consensus_3/coloc_recount/coloc_result_T1D.tsv.gz.Rdata")

B <- as.grob(t1d_recount$`6__90206569__90271941__ENSG00000112182`)
B = ggarrange(B, ncol = 1, nrow = 1) + annotate(geom="text", x=0.8, y=0.55, label= "rs60066732",color="black") + annotate(geom="text", x=0.15, y=0.97, label= "T1D GWAS - circBACH2",color="black")

pdf( file= "output_paper/Main_fig2.pdf",  width= 10, height= 10)
ggarrange(A, B,  heights= c(1, 0.8), nrow = 2, labels = c("A.", "B."))
dev.off()

