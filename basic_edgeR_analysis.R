#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

# ----------------------------------------------------------------------------------------
# --- Do basic analysis of RNAseq data with edgeR
# ----------------------------------------------------------------------------------------

library(Rsubread)

# source("https://bioconductor.org/biocLite.R")
# biocLite("edgeR")
library(edgeR)

bams = list.files(path = "results/bam", pattern = "*.bam$", full.names=TRUE)
bams = bams[!grepl("plus", bams) & !grepl("min", bams)]

gtf.file = "results/new_annotation/all_transcripts.gtf"

fc = featureCounts(bams, annot.ext=gtf.file,
                   isGTFAnnotationFile=TRUE,
                   isPaired=TRUE)

save.image("fc.Rdata")

fc.dge = DGEList(counts=fc$counts, genes=fc$annotation)

temps = read.table("temp_list.notactuallyPRbutsampnums.txt")

# Replace PR, which is incorrect, with SAMP designation
row.names(fc.dge$samples) = gsub("PR", "SAMP", row.names(fc.dge$samples))

sample.temps = sapply(row.names(fc.dge$samples), function (x) {
    ind = gsub(".*SAMP(.*)\\.bam", "SAMP\\1", x)
    temp=temps[temps$V1 == ind,]$V2
})

# row.names(fc.dge$samples) == names(sample.temps)

fc.dge$samples$group = as.factor(sample.temps)

# Filter for genes with low counts across conditions

keep = rowSums(cpm(fc.dge)>1) >= 2
# table(keep)
fc.dge = fc.dge[keep, , keep.lib.sizes=FALSE]

# Normalize for library compositional bias
fc.dge.norm = calcNormFactors(fc.dge)
fc.dge.norm$samples

# Estimate dispersion the simple way
fc.dge.disp = estimateCommonDisp(fc.dge.norm)
fc.dge.disp2 = estimateTagwiseDisp(fc.dge.disp)

# Test for DE genes
et = exactTest(fc.dge.disp2)
topTags(et)

# Estimate dispersion the complicated way (using CR method)
group = as.factor(sample.temps)
design = model.matrix(~group)
fc.dge.disp3 = estimateDisp(fc.dge.norm, design)

# Test for DE genes
fit = glmFit(fc.dge.disp3, design)

# Find genes different between any of the seven groups
lrt = glmLRT(fit, coef=2:7)

# Get top tags
tt = topTags(lrt, n = 1000)$table[,c(1,7:15)]
tt.known = tt[grepl("TsM", tt$GeneID),]

# Find genes different between 37 and 56
lrt = glmLRT(fit, coef=7)

# Get top tags
tt = topTags(lrt, n = 1000)$table[,c(1,7,9:11)]
tt.known = tt[grepl("TsM", tt$GeneID),]

hsps = read.table("hsp_ids.txt")
hsp.idx = which(lrt$genes$GeneID %in% hsps$V1)

hsp.lrt = cbind(lrt$genes[hsp.idx,], lrt$table[hsp.idx,])

# --- Look at two heat and stress genes in particular

tt.known[which(tt.known$GeneID == "TsM_000869600"),]
#              GeneID     logFC       LR       PValue          FDR
# 17321 TsM_000869600 -1.210266 70.99528 3.580804e-17 4.274913e-15
# > 2^-1.210266
# [1] 0.4321889
tt.known[which(tt.known$GeneID == "TsM_001205200"),]
#              GeneID    logFC       LR       PValue         FDR
# 10738 TsM_001205200 1.017383 115.8714 5.071281e-27 3.62343e-24
# 2 ^ 1.017383
# [1] 2.024244

# --- Do some plotting

o = order(lrt$table$PValue)
ct = data.frame(counts = cpm(fc.dge.disp3)[o[which(tt.known$GeneID == "TsM_000869600")],],
	            temp = paste(as.character(group), "C"))
p = ggplot(ct, aes(temp, counts)) +
    geom_boxplot() +
    xlab("Temperature") +
    ylab("Normalized count per million reads") +
    ggtitle("TsM_000869600 - heat shock 70 kda protein") +
    theme_bw()
ggsave(p, file="boxplot.TsM_000869600.pdf")

ct = data.frame(counts = cpm(fc.dge.disp3)[o[which(tt.known$GeneID == "TsM_001205200")],],
                temp = paste(as.character(group), "C"))
p = ggplot(ct, aes(temp, counts)) +
    geom_boxplot() +
    xlab("Temperature") +
    ylab("Normalized count per million reads") +
    ggtitle("TsM_001205200 - Universal stress protein") +
    theme_bw()
ggsave(p, file="boxplot.TsM_001205200.pdf")
