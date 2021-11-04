library(ballgown)
library(metaMA)

data_directory = getwd()
data_directory

bg = ballgown(dataDir=data_directory, samplePattern='tablemaker', meas='all')
bg

bg_filtered <- subset(bg, "rowVars(texpr(bg)) > 1", genomesubset = TRUE)
pData(bg_filtered) = data.frame(id=sampleNames(bg_filtered), group=rep(c(rep("MIM159",3),rep("WT",3))))

structure(bg_filtered_filtered)$exon
structure(bg_filtered)$intron
structure(bg_filtered)$trans

transcript_fpkm = texpr(bg_filtered, 'FPKM')
transcript_cov = texpr(bg_filtered, 'cov')
whole_tx_table = texpr(bg_filtered, 'all')
exon_mcov = eexpr(bg_filtered, 'mcov')
junction_rcount = iexpr(bg_filtered)
whole_intron_table = iexpr(bg_filtered, 'all')
gene_expression = gexpr(bg_filtered)

sampleNames(bg_filtered)                        

exon_transcript_table = indexes(bg_filtered)$e2t
transcript_gene_table = indexes(bg_filtered)$t2g
head(transcript_gene_table)
write.csv(transcript_gene_table, "transcript_gene_table.csv")

phenotype_table = pData(bg_filtered)
write.csv(phenotype_table, "phenotype_table.csv")

plotTranscripts('Nitab4.5_0000681g0250', bg_filtered, 
                samples=c('L1_A_1_tablemaker', 'L2_A_1_tablemaker', 'L3_A_1_tablemaker',
                          'WT1_A_1_tablemaker', 'WT2_A_1_tablemaker', 'WT3_A_1_tablemaker'), 
                meas='FPKM', colorby='transcript')

plotMeans('Nitab4.5_0005563g0030', bg_filtered, groupvar='group', meas='FPKM', colorby='transcript')

stat_results = stattest(bg_filtered, feature='transcript', meas='FPKM', covariate='group', getFC=TRUE)
head(stat_results)
write.csv(stat_results, "stats_results.csv")

####################
bg_filtered_table = texpr(bg_filtered, 'all')
bg_filtered_gene_names = unique(bg_filtered_table[, 9:10])

gene_expression = as.data.frame(gexpr(bg_filtered))
head(gene_expression)

colnames(gene_expression) <- c("MIM159_1", "MIM159_2", "MIM159_3", "WT_1", "WT_2", "WT_3")
row.names(gene_expression)
dim(gene_expression)

data_colors=c("Slategray1", "Slategray2", "Slategray3", "olivedrab2", "olivedrab3", "olivedrab4")

i = row.names(gene_expression) == "Nitab4.5_0005563g0030"
gene_expression[i,]

genes_of_interest = c("Nitab4.5_0003771g0010", "Nitab4.5_0004861g0030",
                      "Nitab4.5_0008527g0010", "Nitab4.5_0005400g0020")
i = which(row.names(gene_expression) %in% genes_of_interest)
gene_expression[i,]

transcript_gene_table = indexes(bg)$t2g
head(transcript_gene_table)

length(row.names(transcript_gene_table))
length(unique(transcript_gene_table[,"g_id"]))

counts=table(transcript_gene_table[,"g_id"])
c_one = length(which(counts == 1))
c_more_than_one = length(which(counts > 1))
c_max = max(counts)
hist(counts, breaks=50, col="bisque4", xlab="Transcripts per gene", main="Distribution of transcript count per gene")
legend_text = c(paste("Genes with one transcript =", c_one), paste("Genes with more than one transcript =", c_more_than_one), paste("Max transcripts for single gene = ", c_max))
legend("topright", legend_text, lty=NULL)

full_table <- texpr(bg , 'all')
hist(full_table$length, breaks=50, xlab="Transcript length (bp)", main="Distribution of transcript lengths", col="steelblue")

min(gene_expression[,"MIM159_1"])
max(gene_expression[,"MIM159_1"])

min(gene_expression[,"MIM159_2"])
max(gene_expression[,"MIM159_2"])

min(gene_expression[,"MIM159_3"])
max(gene_expression[,"MIM159_3"])

min(gene_expression[,"WT_1"])
max(gene_expression[,"WT_1"])

min(gene_expression[,"WT_2"])
max(gene_expression[,"WT_2"])

min(gene_expression[,"WT_3"])
max(gene_expression[,"WT_3"])

min_nonzero=1

data_columns=c(1:6)
short_names=c("MIM159_1", "MIM159_2", "MIM159_3", "WT_1", "WT_2", "WT_3")

boxplot(log2(gene_expression[,data_columns]+min_nonzero), col=data_colors, names=short_names, las=2, ylab="log2(FPKM)", main="Distribution of FPKMs for all 6 libraries")

x = gene_expression[,"MIM159_1"]
y = gene_expression[,"MIM159_2"]
plot(x=log2(x+min_nonzero), y=log2(y+min_nonzero), pch=16, col="blue", cex=0.25, xlab="MIM159_1", ylab="MIM159_2", main="Comparison of expression values for a pair of replicates")
abline(a=0,b=1)
rs=cor(x,y)^2
legend("topleft", paste("R squared = ", round(rs, digits=3), sep=""), lwd=1, col="black")

x = gene_expression[,"WT_1"]
y = gene_expression[,"WT_2"]
plot(x=log2(x+min_nonzero), y=log2(y+min_nonzero), pch=16, col="blue", cex=0.25, xlab="WT_1", ylab="WT_2", main="Comparison of expression values for a pair of replicates")
abline(a=0,b=1)
rs=cor(x,y)^2
legend("topleft", paste("R squared = ", round(rs, digits=3), sep=""), lwd=1, col="black")

colors = colorRampPalette(c("white", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
smoothScatter(x=log2(x+min_nonzero), y=log2(y+min_nonzero), xlab="WT_1", ylab="WT_2", main="Comparison of expression values for a pair of replicates", colramp=colors, nbin=200)

gene_expression[,"sum"]=apply(gene_expression[,data_columns], 1, sum)
i = which(gene_expression[,"sum"] > 5)
r=cor(gene_expression[i,data_columns], use="pairwise.complete.obs", method="pearson")
r

d=1-r
mds=cmdscale(d, k=2, eig=TRUE)
par(mfrow=c(1,1))
plot(mds$points, type="n", xlab="", ylab="", main="MDS distance plot (all non-zero genes) for all libraries", xlim=c(-0.5,0.5), ylim=c(-0.5,0.5))
points(mds$points[,1], mds$points[,2], col="grey", cex=2, pch=16)
text(mds$points[,1], mds$points[,2], short_names, col=data_colors)

results_genes = stattest(bg_filtered, feature="gene", covariate="group", getFC=TRUE, meas="FPKM")
results_genes = merge(results_genes,short_names,by.x=c("id"),by.y=c("gene_id"))

sig=which(results_genes$pval<0.05)
results_genes[,"de"] = log2(results_genes[,"fc"])
hist(results_genes[sig,"de"], breaks=50, col="seagreen", xlab="log2(Fold change) MIM159 vs WT", main="Distribution of differential expression values")
abline(v=-2, col="black", lwd=2, lty=2)
abline(v=2, col="black", lwd=2, lty=2)
legend("topleft", "Fold-change > 4", lwd=2, lty=2)

gene_expression[,"MIM159"]=apply(gene_expression[,c(1:3)], 1, mean)
gene_expression[,"WT"]=apply(gene_expression[,c(3:6)], 1, mean)
x=log2(gene_expression[,"MIM159"]+min_nonzero)
y=log2(gene_expression[,"WT"]+min_nonzero)
plot(x=x, y=y, pch=16, cex=0.25, xlab="MIM159 FPKM (log2)", ylab="WT FPKM (log2)", main="MIM159 vs WT FPKMs")
abline(a=0, b=1)
xsig=x[sig]
ysig=y[sig]
points(x=xsig, y=ysig, col="magenta", pch=16, cex=0.5)
legend("topleft", "Significant (p < 0.05)", col="magenta", pch=16)

sigpi = which(results_genes[,"pval"]<0.05)
sigp = results_genes[sigpi,]
sigde = which(abs(sigp[,"de"]) >= 2)
sig_tn_de = sigp[sigde,]

o = order(sig_tn_de[,"qval"], -abs(sig_tn_de[,"de"]), decreasing=FALSE)
output = sig_tn_de[o,c("gene_name","id","fc","pval","qval","de")]
write.table(output, file="SigDE.txt", sep="\t", row.names=FALSE, quote=FALSE)
#View selected columns of the first 25 lines of output
output[1:25,c(1,4,5)]
