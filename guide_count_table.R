module load R/3.2.5
R

## read in all guide count files
library(data.table)

guide_files=list.files(path=".", pattern=".guides.txt")

guides=list()
for (f in guide_files){
	name=gsub(".guides.txt","",f)
	guides[[name]]=fread(f,sep="\t",h=F)
}

lapply(guides, nrow)
# $HIGH1
# [1] 74266768

# $LOW1
# [1] 86388941

# $MED1
# [1] 73295074

# $REF1
# [1] 68684250

# $TKOvs3_El
# [1] 24119218

### merge gene counts into a count table
gene_counts=lapply(guides, function(x) data.frame(table(x$V3)))

gene_counts_table=gene_counts$TKOvs3_El
## remove the no guide rows
gene_counts_table=gene_counts_table[-1,]
names(gene_counts_table)=c("Gene","TKOvs3")

for (i in 1:(length(gene_counts)-1)){
	names(gene_counts[[i]])=c("Gene",names(gene_counts[i]))
	gene_counts_table=merge(gene_counts_table, gene_counts[[i]], by="Gene", all.x=T)
}
gene_counts_table=gene_counts_table[,c(1,2,6,4,5,3)]
gene_counts_table[is.na(gene_counts_table)] <- 0
write.table(gene_counts_table,"TKOv3.rep1.gene_counts.txt",sep="\t",quote=F, row.names=F)

### normalize to total amount of reads & *10^8
norm_gene_counts=gene_counts_table
for (i in 2:6){
	norm_gene_counts[,i]=sapply(gene_counts_table[,i], function(x) x/sum(gene_counts_table[,i])*10^8)
}
write.table(norm_gene_counts,"TKOv3.rep1.norm_gene_counts.txt",sep="\t",quote=F, row.names=F)


### individual guide counts table
guide_counts=lapply(guides, function(x) data.frame(table(x$V2)))
guide_counts_table=guide_counts$TKOvs3_El
## remove the no guide rows
guide_counts_table=guide_counts_table[-1,]
names(guide_counts_table)=c("Guide","TKOvs3")

for (i in 1:(length(guide_counts)-1)){
	names(guide_counts[[i]])=c("Guide",names(guide_counts[i]))
	guide_counts_table=merge(guide_counts_table, guide_counts[[i]], by="Guide", all.x=T)
}
guide_counts_table=guide_counts_table[,c(1,2,6,4,5,3)]
guide_counts_table[is.na(guide_counts_table)] <- 0
write.table(guide_counts_table,"TKOv3.rep1.guide_counts.txt",sep="\t",quote=F, row.names=F)

### normalize to total amount of reads & *10^8
norm_guide_counts=guide_counts_table
for (i in 2:6){
	norm_guide_counts[,i]=sapply(guide_counts_table[,i], function(x) x/sum(guide_counts_table[,i])*10^8)
}
write.table(norm_guide_counts,"TKOv3.rep1.norm_guide_counts.txt",sep="\t",quote=F, row.names=F)
