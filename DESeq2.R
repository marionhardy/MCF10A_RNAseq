
## This assumes you have run the DESeq2_file_preparation script
# and that you have a counts dataframe and a coldata dataframe

library(DESeq2)
library(tidyverse)
library(cowplot)
library(biomaRt)
library(ggrepel)
library(RColorBrewer)
library(pheatmap)

## Expression analysis of

counts <- as.data.frame(read.csv("./data/MCF10A_counts.csv", row.names = 1))
coldata = read.csv("./data/coldata.csv")
counts = counts[,-1]

# Comparing all conditions
# Create the full model for comparison of samples

dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata,
                              design = ~condition+replicate) 

# Generate a linear model

dds$condition <- relevel(dds$condition, "IM")
dds <- DESeq(dds)

resultsNames(dds)

## Checking distribution

as_tibble(assay(dds)) %>%
  gather(sample, value = counts) %>%
  ggplot(aes(x = log2(counts + 1), fill = sample)) +
  geom_histogram(bins = 20) +
  facet_wrap(~ sample)

# Checking size factors and dispersion

sizeFactors(dds) # only takes into account the sequencing depth, looks ok
plotDispEsts(dds) # verifies normalization, graph looks a-ok

# Checking PCA

rld <- vst(dds)

p1 <- plotPCA(rld,intgroup="condition") + 
  geom_text_repel(max.overlaps = 15,
                  box.padding = 0.25,
                  segment.color = 'grey50',
                  fontface = "italic",
                  label = rld$condition)+
  labs(title = 'PCA per condition')

print(p1)  

# Checking sample similarity

sampleDists <- dist(t(assay(dds)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(dds$condition, sep="-")
colnames(sampleDistMatrix) <- paste(dds$condition, sep="-")
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

