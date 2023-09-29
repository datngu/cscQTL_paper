require(data.table)
library(ggpubr)
require(ggVennDiagram)
require(ggplot2)

# setwd('/Users/datn/github/cscQTL_paper_dev')

read_circ_id_path <- function(path){
  df = fread(path)
  id = df$ID
  id = id[grepl("EN", id)]
  idx = sapply(id, FUN = function(x){unlist(gregexpr("__", x))[3]})
  id = substr(id, 1, idx-1)
  return(id)
}


#########
# read circ_atlas

circ_atlas = as.data.frame(fread("databases/human_bed_v2.0.txt"))
circ_atlas$Chr = gsub("chr", "", circ_atlas$Chro)
circ_atlas$Id = paste(circ_atlas$Chr, circ_atlas$Start, circ_atlas$End, sep = "__")


cscQTL_1 = read_circ_id_path('40_tcells/consensus_1/recount/recount_res_merged_bsj_2.tsv')
cscQTL_2 = read_circ_id_path('40_tcells/consensus_2/recount/recount_res_merged_bsj_2.tsv')
cscQTL_3 = read_circ_id_path('40_tcells/consensus_3/recount/recount_res_merged_bsj_2.tsv')

circall = read_circ_id_path('40_tcells/bovince_circQTL/circRNA_detection/circall/circall_res_merged_bsj_2.tsv')
ciri2 = read_circ_id_path('40_tcells/bovince_circQTL/circRNA_detection/ciri2/ciri2_res_merged_bsj_2.tsv')
circexp2 = read_circ_id_path('40_tcells/bovince_circQTL/circRNA_detection/circexp2/circexp2_res_merged_bsj_2.tsv')
circRNA_finder = read_circ_id_path('40_tcells/bovince_circQTL/circRNA_detection/circRNA_finder/circRNA_finder_res_merged_bsj_2.tsv')

y1 = table(cscQTL_1 %in% circ_atlas$Id)[2]
y2 = table(cscQTL_2 %in% circ_atlas$Id)[2]
y3 = table(cscQTL_3 %in% circ_atlas$Id)[2]

y4 = table(circall %in% circ_atlas$Id)[2]
y5 = table(ciri2 %in% circ_atlas$Id)[2]
y6 = table(circexp2 %in% circ_atlas$Id)[2]
y7 = table(circRNA_finder %in% circ_atlas$Id)[2]

x = list(cscQTL_1 = cscQTL_1, cscQTL_2 = cscQTL_2, cscQTL_3 = cscQTL_3, Circall = circall, CIRI2 = ciri2, CIRCexplorer2 = circexp2, circRNA_finder = circRNA_finder)


methods = names(x)
df = data.frame(Methods = methods)
df$circRNA = lengths(x)
df$Count_CircAtlas = c(y1, y2, y3, y4, y5, y6, y7)


# Fig A
tem = df
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



A = ggplot(data= dfa, aes(x=Methods, y=Count, fill=CircAtlas)) + geom_bar(stat="identity", width=0.7) + geom_text(aes(y = y_pos, label = percent), vjust=1.6, color="white", size=3)

A = A + scale_fill_brewer(palette="Paired") + theme_classic() + theme(legend.position="top", axis.title.x=element_blank()) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  ylab("#circRNA - BSJ >= 2 & expressed in >= 30% samples") 



## Fig B

## read statistics

circall = fread('40_tcells/bovince_circQTL/circRNA_detection/circall.stat')
ciri2 = fread('40_tcells/bovince_circQTL/circRNA_detection/ciri2.stat')
circexp2 = fread('40_tcells/bovince_circQTL/circRNA_detection/circexp2.stat')
circRNA_finder = fread('40_tcells/bovince_circQTL/circRNA_detection/circRNA_finder.stat')

methods = c(rep('Circall', 40), rep('CIRI2', 40), rep('CIRCexplorer2',40), rep('circRNA_finder',40))
counts = c(circall$V1, ciri2$V1, circexp2$V1, circRNA_finder$V1)

df = data.frame(Methods = methods, Counts = counts)

plot = ggplot(df, aes(x = Methods, y = Counts, fill = Methods)) +
  geom_boxplot()


# Calculate mean values by type
mean_values <- aggregate(Counts ~ Methods, data = df, FUN = mean)

# Add the mean labels to the plot
B = plot + geom_text(data = mean_values, aes(label = round(Counts)), vjust = -0.8, size=3) + theme_classic() + theme(legend.position = "none") + ylab("#circRNA with BSJ >= 2 detected") + theme(axis.text.x = element_text(angle = 45, hjust = 1))


pdf( file= "output_paper/S2_revise.pdf",  width= 10, height= 6)
#ggarrange(A, B,  heights= c(1, 1), nrow = 2, labels = c("A.", "B."))
ggarrange(A, B,  widths= c(1, 1), ncol = 2, labels = c("A.", "B."))

dev.off()