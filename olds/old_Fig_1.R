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





# ggven <- function(x){
#   venn <- Venn(x)
#   data <- process_data(venn)
#   df = venn_region(data)
#   df$percent = round(df$count/ sum(df$count),4)*100
#   df$my_lab = paste0(df$count, "\n(", df$percent, "%)")

#   p = ggplot() + geom_sf(aes(fill=count), data = df) + geom_sf(size = 2, lty = "solid", color = "gray", data = venn_setedge(data), show.legend = F) + geom_sf_text(aes(label = name), data = venn_setlabel(data)) + geom_sf_label(aes(label=my_lab), size = 3.5, data = df, alpha = 0, label.size  = NA) + theme_void() + theme(legend.position = "none")
#   p = p + scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")
#   p = p + scale_x_continuous(expand = expansion(mult = .2))
#   return(p)
# }




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
F1a = ggven(f1ax[1:3])


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


dfb = rbind(tem3, tem2, tem1, tem_unique, tem_unique2)

dfb$Methods = factor(dfb$Methods, levels = names(f1bx))
dfb$Type = factor(dfb$Type, levels = rev(c("Supported by all 3 single methods", "Supported by 2 single methods", "Supported by 1 single methods",  "cscQTL specific", "Single method")))



# tem3 = get_circQTL_type(f1bx, Supporting_method_3)
# tem3$Type = "Supported by all 3 single methods"
# tem3$Count = tem3$Cumsum

# tem2 = get_circQTL_type(f1bx, Supporting_method_2)
# tem2$Type = "Supported by 2 single methods"
# tem2$Count = tem2$Cumsum - tem3$Cumsum

# tem1 = get_circQTL_type(f1bx, Supporting_method_1)
# tem1$Type = "Supported by 1 single methods"
# tem1$Count = tem1$Cumsum - tem2$Cumsum


# tem_unique = get_circQTL_all(f1bx)
# tem_unique$Type = "Method specific"
# tem_unique$Count = tem_unique$Cumsum - tem1$Cumsum


# # tem_unique2 = get_circQTL_all(f1bx[4:6])
# # tem_unique2$Type = "Method specific"
# # tem_unique2$Count = tem_unique2$Cumsum


# dfb = rbind(tem3, tem2, tem1, tem_unique)

# dfb$Methods = factor(dfb$Methods, levels = names(f1bx))
# dfb$Type = factor(dfb$Type, levels = rev(c("Supported by all 3 single methods", "Supported by 2 single methods", "Supported by 1 single methods",  "Method specific")))

F1b = ggplot(data= dfb, aes(x=Methods, y=Count, fill=Type)) + geom_bar(stat="identity", width=0.5) + geom_text(aes(y = Cumsum, label = Count), vjust=1.6, color="white", size=3.5)

F1b = F1b + theme_classic() + theme(legend.position="left", axis.title.x=element_blank()) + theme(legend.position = c(0.8, 0.8)) + theme(legend.text=element_text(size=6))





pdf( file= "output_paper/S1_new.pdf",  width= 10, height= 6)
ggarrange(F1a, F1b,  widths= c(1, 1), ncol = 2, labels = c("A.", "B."))
dev.off()




######## Write table

methods = c("Circall", "CIRI2", "CIRCexplorer2", "Supporting_method_1", "Supporting_method_2", "Supporting_method_3")
tab1 = data.frame(Methods = methods)
# tab1$Count_CircAtlas = c(y1, y2, y3, y4, y5, y6)
tab1$circQTL = lengths(f1ax)

fwrite(tab1, "output_paper/table1.csv")


















