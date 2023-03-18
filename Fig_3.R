setwd("/Users/datn/github/circQTL_analysis")

require(data.table)
library(ggpubr)
require(ggplot2)
library(scales)
library("ggplotify")
require(ggVennDiagram)


read_normial_qtl <- function(dir, chr_list = as.character(c(1:22))){
	res = fread(paste0(dir,"/chr_1_nominal_pass_best_peer_factors.tsv.gz"))
	for(i in 2:length(chr_list)){
		chr = chr_list[i]
		tem = fread(paste0(dir,"/chr_", chr ,"_nominal_pass_best_peer_factors.tsv.gz"))
		res = rbind(res, tem)
	}
	return(res)
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




csc = fread("consensus_3/qtl_mapping_apply_qvalue/recount.tsv")

# csc1 = fread("consensus_1/qtl_mapping_apply_qvalue/recount.tsv")





# Fig A. Q-Q plot of p-values

pnorm_cut = max(csc$V6[which(csc$V18 == max(csc$V18))])

norm_pass_all = norm_pass = read_normial_qtl("consensus_3/recount_qtl_mapping_nominal")

norm_pass = norm_pass[norm_pass$V1 %in% csc$V1,]
norm_pass = norm_pass[norm_pass$V12 < pnorm_cut,]

obs_all = -log10(sort(as.numeric(norm_pass_all$V12)))
# reduced size by 50
obs = obs_all[seq(1,length(obs_all),50)]

exp = -log10(ppoints(length(obs))/10)

df = data.frame(obs = obs, exp = exp)

A = ggplot(data = df, aes(x = exp, y = obs)) + geom_point(shape = 21, size = 1.2, color = "#4682B4", alpha = 1) + geom_abline(intercept = 0, slope = 1, size = 1)

log10Pe <- expression(paste("Expected -log"[10], plain(P)))
log10Po <- expression(paste("Observed -log"[10], plain(P)))

A = A + theme_classic() + ylab(log10Po) + xlab(log10Pe)






# Fig B. effect size circRNAs vs parent mRNA

esnp_circ = norm_pass
esnp_circ$geneID = sapply(esnp_circ$V1, FUN = function(x){unlist(strsplit(x, "__"))[4]})


esnp_mrna = read_normial_qtl("consensus_3/salmon_qtl_mapping_nominal")
esnp_mrna = esnp_mrna[esnp_mrna$V8 %in% esnp_circ$V8]

esnp_circ$circRNA_slope = esnp_circ$V14
esnp_circ$mRNA_slope = esnp_mrna$V14[match(esnp_circ$geneID, esnp_mrna$V1)]

esnp_circ = esnp_circ[!is.na(esnp_circ$mRNA_slope),]

B = ggplot(data = esnp_circ, aes(x = circRNA_slope, y = mRNA_slope)) + geom_point(shape = 21, size = 1.2, color = "#4682B4", alpha = 1) + ylab("Estimated slope of parent linear genes") + xlab("Estimated slope of circRNAs") 

m = lm(mRNA_slope ~ circRNA_slope, data = esnp_circ)
c = cor.test(esnp_circ$circRNA_slope, esnp_circ$mRNA_slope, method = "pearson")
print(c)
c = round(c$estimate,3)
text = paste0("R=", c, "\np-value < 2.2e-16")
B = B + geom_abline(intercept = m$coefficients[1], slope = m$coefficients[2], color="black", size=1) 
B = B + annotate(geom="text", x=-0.7, y=1.25, label= text,color="black")
B = B + theme_classic()
B = B + theme(legend.position = c(0.85, 0.125))


# Fig C. Distance distribution of eSNPs

## counting before ploting
norm_pass$dis = norm_pass$V7/1e3

cuts = seq(-1e3, 1e3, 2e2)
lab = paste0( "(", cuts[1:10], ":", cuts[2:11], "]")
df = data.frame(Distance = cuts)
df$Count = 0


for(i in 1:length(cuts)){
	dis = norm_pass$dis
	dis = dis[dis > cuts[i]]
	dis = dis[dis <= cuts[i+1]]
	df$Count[i] = length(dis)
}
df = df[1:10,]
df$Distance = lab



positions = as.character(df$Distance[1:10])




C = ggplot(norm_pass, aes(x=dis)) + geom_histogram(color="black", fill="#4682B4", binwidth = 20) + theme_classic() + ylab("Number of circSNPs") + xlab("Distance (kilobase)") + xlim(-1000,1000) + scale_x_continuous(breaks = seq(-1000, 1000, 200)) 

#
# Fig D. VEP annotation of eSNPs

### export eSNPs for VEP analyses
# fwrite(list(norm_pass$V8), file = "output_paper/cscQTL_all_eSNP.txt", sep = "\t")
# esnp_all = norm_pass$V8
vep = fread("databases/VEP_cscQTL_eSNP.txt")
dup = duplicated(vep$`#Uploaded_variation`)
vep = vep[!dup,]

count = sort(table(vep$Consequence), decreasing = T)

count_plot = count[1:9]
count_names = gsub("_variant", "", names(count_plot))
count_names = gsub("_exon", "", count_names)
count_names = gsub("_", " ", count_names)
count_names = c(count_names, "others")
count_plot = c(count_plot, sum(count[10:length(count)]))

count_pct = round(count_plot/sum(count_plot)*100,2)
count_names = paste0(count_names, " (", count_pct, "%)")

names(count_plot) = count_names

#positions = as.character(count_names)

df = data.frame(Group = count_names, Count = count_plot, Count_pct = count_pct)
df$Group = factor(df$Group, levels = count_names)
# Barplot
bp <- ggplot(df, aes(x="", y=Count, fill=Group))+
geom_bar(width = 1, stat = "identity", color = "gray", alpha = 0.5)
pie <- bp + coord_polar("y", start=0)

D = pie + theme_void() + theme(axis.text.x=element_blank())


#pie + geom_label_repel(aes(label = count_pct), size=5, show.legend = F, nudge_x = 1) + guides(fill = guide_legend(title = "Group"))

###### E. Counting gene containing eCircRNAs

csc3 = fread("consensus_3/qtl_mapping_apply_qvalue/recount.tsv")
idx = sapply(csc3$V1, FUN = function(x){unlist(gregexpr("__", x))[3]})
csc3$gene = substr(csc3$V1, idx+2, nchar(csc3$V1))

salmon = fread("consensus_3/qtl_mapping_apply_qvalue/salmon.tsv")


table(csc3$gene %in% salmon$V1)


x0 = list("eCircQTL host genes" = csc3$gene, "eQTL eGenes" = salmon$V1)

E = ggven(x0)


# r12 = ggarrange(A, D, C, B, widths= c(1, 1), ncol = 2, nrow = 2, labels = c("A.", "B.", "C.", "D."), align = "h")
# r3 = ggarrange(E, ncol = 2, widths= c(1.5, 0.5), nrow = 1, labels = "E")


r12 = ggarrange(A, D, C, B, widths= c(1, 1), ncol = 2, nrow = 2, labels = c("A.", "B.", "C.", "D."), align = "h")
r3 = ggarrange(E, ncol = 2, widths= c(1.5, 0.5), nrow = 1, labels = "E")

## merging and export figure
pdf( file= "output_paper/Fig3.pdf",  width= 10, height= 12)
#ggarrange(r12, r3, nrow = 2, heights = c(1, 0.4), align = "v")
ggarrange(A, D, C, E, B, widths= c(1, 1), ncol = 2, nrow = 3, labels = c("A.", "B.", "C.", "D.", "E"))
dev.off()










