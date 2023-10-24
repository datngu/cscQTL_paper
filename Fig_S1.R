# setwd("/Users/datn/github/circQTL_analysis")


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


ggven <- function(x){
  venn <- Venn(x)
  data <- process_data(venn)
  df = venn_region(data)
  df$percent = round(df$count/ sum(df$count),4)*100
  df$my_lab = paste0(df$count, "\n(", df$percent, "%)")

  p = ggplot() + geom_sf(aes(fill=id, alpha = 0.5), data = df) + geom_sf(size = 2, lty = "solid", color = "gray", data = venn_setedge(data), show.legend = F) + geom_sf_text(aes(label = name), data = venn_setlabel(data)) + geom_sf_label(aes(label=my_lab), size = 2.5, data = df, alpha = 0, label.size  = NA) + theme_void() + theme(legend.position = "none")
  #p = p + scale_color_brewer(palette="Dark2") 
  p = p + scale_x_continuous(expand = expansion(mult = .2))
  return(p)
}


#########
## single method overlaping

### read circIDs of circQTL

circall = unique(read_circ_id_qtl("circall", dir = "40_tcells/consensus_3"))
ciri2 = unique(read_circ_id_qtl("ciri2", dir = "40_tcells/consensus_3"))
circexp2 = unique(read_circ_id_qtl("circexp2", dir = "40_tcells/consensus_3"))

all = c(circall, ciri2, circexp2)
t = table(all)

Supporting_method_1 = names(t)[t>=1]
Supporting_method_2 = names(t)[t>=2]
Supporting_method_3 = names(t)[t>=3]


f1ax = list(Circall = circall, CIRI2 = ciri2, CIRCexplorer2 = circexp2, Supporting_method_1 = Supporting_method_1, Supporting_method_2 = Supporting_method_2, Supporting_method_3 = Supporting_method_3)


## fig 1 a
A = ggven(f1ax[1:3])

pdf( file= "output_paper/S1_revise.pdf",  width= 10, height= 12)
A
dev.off()



######## Write table

methods = c("Circall", "CIRI2", "CIRCexplorer2", "Supporting_method_1", "Supporting_method_2", "Supporting_method_3")
tab1 = data.frame(Methods = methods)
# tab1$Count_CircAtlas = c(y1, y2, y3, y4, y5, y6)
tab1$circQTL = lengths(f1ax)

fwrite(tab1, "output_paper/table1.csv")


















