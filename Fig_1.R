setwd("/Users/datn/github/circQTL_analysis")


require(data.table)
library(ggpubr)
require(ggVennDiagram)
require(ggplot2)



read_circ_id <- function(tool, dir = "consensus_3"){
  df = fread(paste0(dir, "/quantile_norm/", tool, "_quantile_norm.tsv"))
  id = df$TargetID
  id = id[grepl("EN", id)]
  idx = sapply(id, FUN = function(x){unlist(gregexpr("__", x))[3]})
  id = substr(id, 1, idx-1)
  return(id)
}


read_circ_id_qtl <- function(tool, dir = "consensus_3"){
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

  p = ggplot() + geom_sf(aes(fill=count), data = df) + geom_sf(size = 2, lty = "solid", color = "gray", data = venn_setedge(data), show.legend = F) + geom_sf_text(aes(label = name), data = venn_setlabel(data)) + geom_sf_label(aes(label=my_lab), size = 3.5, data = df, alpha = 0, label.size  = NA) + theme_void() + theme(legend.position = "none")
  p = p + scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")
  p = p + scale_x_continuous(expand = expansion(mult = .2))
  return(p)
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
# read circ_atlas

circ_atlas = as.data.frame(fread("databases/human_bed_v2.0.txt"))
circ_atlas$Chr = gsub("chr", "", circ_atlas$Chro)
circ_atlas$Id = paste(circ_atlas$Chr, circ_atlas$Start, circ_atlas$End, sep = "__")


## read circIDs of methods

circall = unique(read_circ_id("circall"), dir = "consensus_3")
ciri2 = unique(read_circ_id("ciri2"), dir = "consensus_3")
circexp2 = unique(read_circ_id("circexp2"), dir = "consensus_3")

all = c(circall, ciri2, circexp2)
t = table(all)

Supporting_method_1 = names(t)[t>=1]
Supporting_method_2 = names(t)[t>=2]
Supporting_method_3 = names(t)[t>=3]

table(circall %in% circ_atlas$Id)[1]/sum(table(circall %in% circ_atlas$Id))*100
table(ciri2 %in% circ_atlas$Id)[1]/sum(table(ciri2 %in% circ_atlas$Id))*100
table(circexp2 %in% circ_atlas$Id)[1]/sum(table(circexp2 %in% circ_atlas$Id))*100

y1 = table(circall %in% circ_atlas$Id)[2]
y2 = table(ciri2 %in% circ_atlas$Id)[2]
y3 = table(circexp2 %in% circ_atlas$Id)[2]
y4 = table(Supporting_method_1 %in% circ_atlas$Id)[2]
y5 = table(Supporting_method_2 %in% circ_atlas$Id)[2]
y6 = table(Supporting_method_3 %in% circ_atlas$Id)[2]

x = list(Circall = circall, CIRI2 = ciri2, CIRCexplorer2 = circexp2, Supporting_method_1 = Supporting_method_1, Supporting_method_2 = Supporting_method_2, Supporting_method_3 = Supporting_method_3)


methods = c("Circall", "CIRI2", "CIRCexplorer2", "Supporting_method_1", "Supporting_method_2", "Supporting_method_3")
df = data.frame(Methods = methods)
df$circRNA = lengths(x)
df$Count_CircAtlas = c(y1, y2, y3, y4, y5, y6)


# Fig 1.A
tem = df[1:3,]
tem1 = tem[,c(1,3)]
tem1$CircAtlas = "Yes"
names(tem1) = c("Methods", "Count", "CircAtlas")
tem1$percent = round(tem1$Count/ tem$circRNA * 100,2)
tem1$percent = paste0(tem1$percent, "%")
tem1$y_pos = tem1$Count


tem2 = tem[,1:2]
tem2[,2] = tem[,2] - tem[,3]
tem2$CircAtlas = "No"
names(tem2) = c("Methods", "Count", "CircAtlas")
tem2$percent = ""
tem2$y_pos = 0

dfa = rbind(tem1, tem2)
dfa$Methods = factor(dfa$Methods, levels = tem2$Methods)



F1a = ggplot(data= dfa, aes(x=Methods, y=Count, fill=CircAtlas)) + geom_bar(stat="identity", width=0.7) + geom_text(aes(y = y_pos, label = percent), vjust=1.6, color="white", size=3.5)

F1a = F1a + scale_fill_brewer(palette="Paired") + theme_classic() + theme(legend.position="top", axis.title.x=element_blank())


# Fig 1.B

F1b = ggven(x[1:3]) 


# Fig 1.C
### read circIDs of circQTL
circall = unique(read_circ_id_qtl("circall"))
ciri2 = unique(read_circ_id_qtl("ciri2"))
circexp2 = unique(read_circ_id_qtl("circexp2"))

all = c(circall, ciri2, circexp2)
t = table(all)

Supporting_method_1 = names(t)[t>=1]
Supporting_method_2 = names(t)[t>=2]
Supporting_method_3 = names(t)[t>=3]


x = list(Circall = circall, CIRI2 = ciri2, CIRCexplorer2 = circexp2, Supporting_method_1 = Supporting_method_1, Supporting_method_2 = Supporting_method_2, Supporting_method_3 = Supporting_method_3)

df$circQTL = lengths(x)

F1c = ggven(x[1:3])


## merge figures
c1 = F1a
c2 = ggarrange(F1b, F1c, nrow = 2, labels = c("B.", "C."))

m = ggarrange(c1, c2,  widths= c(0.7, 1), ncol = 2, nrow = 1, align = "v", labels = "A." )

pdf( file= "output_paper/Fig1.pdf",  width= 10, height= 8)
m
dev.off()



######## Write data
fwrite(df, "output_paper/table1.csv")


















