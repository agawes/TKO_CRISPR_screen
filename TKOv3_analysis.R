# cd /well/mccarthy/production/genetic-screens/data/EndoC_CRISPR_TKO/190819_NB502094_0138_AH5YHWBGXC/guide_counts/
# module load R/3.2.5

## analysis 
genes=read.table("TKOv3.rep1.norm_gene_counts.txt",sep="\t",h=T)
guides=read.table("TKOv3.rep1.norm_guide_counts.txt",sep="\t",h=T)

## look at distribution of counts per gene and per guide
pdf("gene_counts_distribution.pdf")
plot(density(genes$TKOvs3), col="darkgrey", lwd=2, lty=2, xlim=c(0,15000), xlab="Norm counts/gene", main="")
lines(density(genes$REF1), col="black", lwd=2, lty=2)
lines(density(genes$LOW1), col="blue", lwd=2)
lines(density(genes$MED1), col="grey", lwd=2)
lines(density(genes$HIGH1), col="red", lwd=2)
legend("topright",c("TKOvs3","REF1","LOW1","MED1","HIGH1"), col=c("darkgrey","black","blue","grey","red"),
lty=c(2,2,1,1,1), lwd=2, bty="n")
dev.off()

pdf("guide_counts_distribution.pdf")
plot(density(guides$TKOvs3), col="darkgrey", lwd=2, lty=2, xlim=c(0,6000), xlab="Norm counts/guide", main="")
lines(density(guides$REF1), col="black", lwd=2, lty=2)
lines(density(guides$LOW1), col="blue", lwd=2)
lines(density(guides$MED1), col="grey", lwd=2)
lines(density(guides$HIGH1), col="red", lwd=2)
legend("topright",c("TKOvs3","REF1","LOW1","MED1","HIGH1"), col=c("darkgrey","black","blue","grey","red"),
lty=c(2,2,1,1,1), lwd=2, bty="n")
dev.off()


## look for essential genes - genes present in the TKOvs3 but absent in REF1
## These would be guides (sgRNAs) that produced a knockout 
## of any essential gene and they've caused cell death so we won't be able to see those guides.

pdf("TKOv3_vs_REF1.scatter.pdf")
plot(genes$TKOvs3, genes$REF1, pch=16, col="black", xlab="Gene counts - TKOvs3 library",
ylab="Gene counts - REF1")
text(68500,75000, "LacZ")
text(28000,35000, "EGFP")
text(37000,30000, "luciferase")
plot(genes$TKOvs3, genes$REF1, pch=16, col="black", xlim=c(0,15000), ylim=c(0,15000),
xlab="Gene counts - TKOvs3 library", ylab="Gene counts - REF1")

plot(guides$TKOvs3, guides$REF1, pch=16, col="black", xlab="Guide counts - TKOvs3 library",
ylab="Guide counts - REF1")
dev.off()

pdf("REF1_vs_TKOvs3.FC_distribution.pdf")
plot(density(log2(guides$REF1/guides$TKOvs3)), xlab="Log2 Fold Change: REF1/TKOvs3", main="")
abline(v=c(-1,1), col="red", lty=2)
dev.off()



## comparison of samples - log2

## log2FC data frame
gene_log2FC=data.frame(Gene=genes$Gene,
 REF1_vs_TKO=log2(genes$REF1/genes$TKOvs3), HIGH1_vs_LOW1=log2(genes$HIGH1/genes$LOW1),
 HIGH1_vs_MED1=log2(genes$HIGH1/genes$MED1), MED1_vs_LOW1=log2(genes$MED1/genes$LOW1))

write.table(gene_log2FC,"genes_comparisons.log2FC.txt",sep="\t", quote=F, row.names=F)

pdf("gene.log2FC_distribution.pdf")
plot(density(gene_log2FC$REF1_vs_TKO, na.rm=T), xlab="Log2 Fold Change", main="", col="darkgrey")
lines(density(gene_log2FC$HIGH1_vs_LOW1,na.rm=T), col="darkred")
lines(density(gene_log2FC$HIGH1_vs_MED1, na.rm=T), col="red")
lines(density(gene_log2FC$MED1_vs_LOW1, na.rm=T), col="blue")
abline(v=c(-1,1), col="red", lty=2)
legend("topright", c("REF1 vs TKO","HIGH1 vs LOW1", "HIGH1 vs MED1", "MED1 vs LOW1"), 
col=c("darkgrey","darkred","red","blue"), lty=1, lwd=2, bty="n", title="Comparison")
dev.off()

## log2FC data frame
guide_log2FC=data.frame(Guide=guides$Guide, Gene=gsub("_.+","",guides$Guide),
 REF1_vs_TKO=log2(guides$REF1/guides$TKOvs3), HIGH1_vs_LOW1=log2(guides$HIGH1/guides$LOW1),
 HIGH1_vs_MED1=log2(guides$HIGH1/guides$MED1), MED1_vs_LOW1=log2(guides$MED1/guides$LOW1))

write.table(guide_log2FC,"guide_comparisons.log2FC.txt",sep="\t", quote=F, row.names=F)

pdf("guide.log2FC_distribution.pdf")
plot(density(guide_log2FC$REF1_vs_TKO, na.rm=T), xlab="Log2 Fold Change", main="", col="darkgrey")
lines(density(guide_log2FC$HIGH1_vs_LOW1,na.rm=T), col="darkred")
lines(density(guide_log2FC$HIGH1_vs_MED1, na.rm=T), col="red")
lines(density(guide_log2FC$MED1_vs_LOW1, na.rm=T), col="blue")
abline(v=c(-1,1), col="red", lty=2)
legend("topright", c("REF1 vs TKO","HIGH1 vs LOW1", "HIGH1 vs MED1", "MED1 vs LOW1"), 
col=c("darkgrey","darkred","red","blue"), lty=1, lwd=2, bty="n", title="Comparison")
dev.off()


### plot norm gene counts & log2FC for select genes
guides$Gene=gsub("_.+","",guides$Guide)

gene="INS"
pdf(paste0(gene,".gene_barplot.pdf"), height=7)
par(mfrow=c(2,1), las=2, mar=c(5.1,5.1,3,1), oma=c(4.5,0,0,0))
barplot(as.numeric(subset(genes, Gene==gene)[,2:6]), names.arg=names(genes)[2:6], 
 ylab="Norm gene counts", cex.names=1.2, cex.axis=1.2, cex.lab=1.2, 
 col=c("grey","black","blue","white","red"), main=gene)
  abline(h=0)

barplot(as.numeric(subset(gene_log2FC, Gene==gene)[,2:5]), names.arg=names(gene_log2FC)[2:5],
 ylab="Gene log2FC", cex.names=1.2, cex.axis=1.2, cex.lab=1.2, col=c("grey","red","coral","blue"))
 abline(h=0)
dev.off()

pdf(paste0(gene,".guide_barplot.pdf"), height=7)
par(mfrow=c(2,1), las=2, mar=c(5.1,5.1,3,1), oma=c(4.5,0,0,0))
m=data.matrix(subset(guides, Gene==gene)[,2:6])
rownames(m)=paste0(gene,"_",1:nrow(m))
barplot(m, beside=T,
 ylab="Norm guide counts", cex.names=1.2, cex.axis=1.2, cex.lab=1.2, 
 col=c(rep("grey",4),rep("black",4),rep("blue",4),rep("white",4),rep("red",4)), main=gene)
  abline(h=0)

m=data.matrix(subset(guide_log2FC, Gene==gene)[,3:6])
rownames(m)=paste0(gene,"_",1:nrow(m))

barplot(m, beside=T,
 ylab="Guide log2FC", cex.names=1.2, cex.axis=1.2, cex.lab=1.2, 
 col=c(rep("grey",4),rep("red",4),rep("coral",4),rep("blue",4)))
 abline(h=0)
dev.off()

