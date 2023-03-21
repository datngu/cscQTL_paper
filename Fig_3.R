setwd("/Users/datn/github/cscQTL_paper")

require(data.table)
library(ggpubr)
library(VennDiagram)
require(ggVennDiagram)
require(ggplot2)



read_circ_id_recount <- function(tool, dir = "consensus_3"){
  fn = paste0(dir, "/quantile_norm/", tool, "_quantile_norm.tsv")
  df = fread(fn)
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

  p = ggplot() + geom_sf(aes(fill=count), data = df) + geom_sf(size = 2, lty = "solid", color = "gray", data = venn_setedge(data), show.legend = F) + geom_sf_text(aes(label = name), data = venn_setlabel(data)) + geom_sf_label(aes(label=my_lab), size = 2.5, data = df, alpha = 0, label.size  = NA) + theme_void() + theme(legend.position = "none")
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



display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}



# read circ_atlas

circ_atlas = as.data.frame(fread("databases/human_bed_v2.0.txt"))
circ_atlas$Chr = gsub("chr", "", circ_atlas$Chro)
circ_atlas$Id = paste(circ_atlas$Chr, circ_atlas$Start, circ_atlas$End, sep = "__")


## read circIDs of methods

recount_1 = unique(read_circ_id_recount("recount", dir = "consensus_1"))
recount_2 = unique(read_circ_id_recount("recount", dir = "consensus_2"))
recount_3 = unique(read_circ_id_recount("recount", dir = "consensus_3"))



y1 = table(recount_1 %in% circ_atlas$Id)[2]
y2 = table(recount_2 %in% circ_atlas$Id)[2]
y3 = table(recount_3 %in% circ_atlas$Id)[2]
z1 = length(recount_1)
z2 = length(recount_2)
z3 = length(recount_3)

tem1 = data.frame(Methods = c("cscQTL_1", "cscQTL_2", "cscQTL_3"))
tem1$Count = c(y1, y2, y3)
tem1$CircAtlas = "Yes"

tem1$percent = round(tem1$Count/ c(z1,z2,z3)* 100,2)
tem1$percent = paste0(tem1$percent, "%")
tem1$y_pos = tem1$Count


tem2 = data.frame(Methods = c("cscQTL_1", "cscQTL_2", "cscQTL_3"))
tem2$Count = c(z1-y1, z2-y2, z3-y3)
tem2$CircAtlas = "No"
tem2$percent = ""
tem2$y_pos = 0

df = rbind(tem1, tem2)

S1a = ggplot(data= df, aes(x=Methods, y=Count, fill=CircAtlas)) + geom_bar(stat="identity", width=0.5) + geom_text(aes(y = y_pos, label = percent), vjust=1.6, color="white", size=3.5)

S1a = S1a + scale_fill_brewer(palette="Paired") + theme_classic() + theme(legend.position="top", axis.title.x=element_blank())






########## eQTL results

circall = unique(read_circ_id_qtl("circall", dir = "consensus_3"))
ciri2 = unique(read_circ_id_qtl("ciri2", dir = "consensus_3"))
circexp2 = unique(read_circ_id_qtl("circexp2", dir = "consensus_3"))
recount_1 = unique(read_circ_id_qtl("recount", dir = "consensus_1"))
recount_2 = unique(read_circ_id_qtl("recount", dir = "consensus_2"))
recount_3 = unique(read_circ_id_qtl("recount", dir = "consensus_3"))


#x = list(Circall = circall, CIRI2 = ciri2, CIRCexplorer2 = circexp2, ReCount_1 = recount_1, ReCount_2 = recount_2)

#y = list(ReCount_1 = recount_1, ReCount_2 = recount_2)

#display_venn(y, fill = c("#999999", "#E69F00"))

x = list("cscQTL_1" = recount_1, "cscQTL_2" = recount_2, "cscQTL_3" = recount_3, Circall = circall, CIRI2 = ciri2, CIRCexplorer2 = circexp2)

dfb = data.frame(Methods = names(x), Count = lengths(x))

dfb$Methods = factor(dfb$Methods, levels = dfb$Methods)


x0 = list("cscQTL_1" = recount_1, "cscQTL_2" = recount_2, "cscQTL_3" = recount_3)

x1 = list(Circall = circall, CIRI2 = ciri2, CIRCexplorer2 = circexp2, "cscQTL_1" = recount_1)

x2 = list(Circall = circall, CIRI2 = ciri2, CIRCexplorer2 = circexp2, "cscQTL_2" = recount_2)

x3 = list(Circall = circall, CIRI2 = ciri2, CIRCexplorer2 = circexp2, "cscQTL_3" = recount_3)





S1b = ggven(x0)

S1c = ggven(x1)

S1d = ggven(x2)

#S1e = ggven(x3)


#t = as.grob(display_venn(x0, fill = c("#E69F00", "#56B4E9", "#009E73"))

## merge figures
# c1 = S1a

# c2 = ggarrange(S1b, S1c, nrow = 2, labels = c("B.", "C."))

# m = ggarrange(c1, c2,  widths= c(0.7, 1), ncol = 2, nrow = 1, align = "v", labels = "A." )

m = ggarrange(S1a, S1b, S1c, S1d,  widths= c(1, 1), ncol = 2, nrow = 2, align = "v", labels = c("A.", "B.", "C.", "D." ))

pdf( file= "output_paper/S1.pdf",  width= 10, height= 8)
m
dev.off()

# m = ggarrange(S1a, S1b, widths= c(0.7, 1), ncol = 2, nrow = 1, align = "v", labels = c("A.", "B." ))

# pdf( file= "output_paper/S1.pdf",  width= 10, height= 4)
# m
# dev.off()


F2a = ggplot(data= dfb, aes(x=Methods, y=Count)) + geom_bar(stat="identity", fill="steelblue", width=0.7) + geom_text(aes(y = Count, label = Count), vjust=1.6, color="white", size=3.5) + guides(x = guide_axis(angle = 45))  + theme_classic() + theme(legend.position="top", axis.title.x=element_blank())
F2b = ggven(x3)

pdf( file= "output_paper/Fig3.pdf",  width= 10, height= 5)
ggarrange(F2a, F2b,  widths= c(0.7, 1), ncol = 2, labels = c("A.", "B."))
dev.off()



# pdf( file= "output_paper/F2b.pdf",   width= 7, height= 5)
# display_venn(x1, fill = c("#999999", "#E69F00", "#56B4E9", "#009E73"))
# dev.off()

# pdf( file= "output_paper/F2c.pdf",   width= 7, height= 5)
# display_venn(x2, fill = c("#999999", "#E69F00", "#56B4E9", "#009E73"))
# dev.off()

