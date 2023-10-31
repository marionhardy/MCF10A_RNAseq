# MCF10A_RNAseq

Albeck lab MCF10A screening of a bunch of metabolic conditions

Sequenced through Parse bio

## Experimental design

MCF10A in culture
Biological replicates per method : 2, made two weeks apart
- Live imaged
- RNA sequenced
- scRNA sequenced

Timing of spikes then collecting: 16 hours after spike

## Conditions tested

- Growth medium (horse serum, EGF, insulin, cholera toxin, HC, glucose, glutamine, pyruvate)
- Imaging medium (insulin, cholera toxin, HC, glucose, pyruvate)
- IM + glutamine
- IM + EGF
- IM - Insulin
- IM - HC
- IM - CT
- IM + AMPK activator MK8722
- IM + AKT inhibitor Ipasertib
- IM + ERK inhibitor PD-0325901
- IM + mTORC1 inhibitor Rapamycin
- IM + Oligomycin
- IM + MPC inhibitor UK5099
- IM + LDH inhibitor Galloflavin
- IM + IL6

For the concentrations, see the cell culture help chart on the data server in Albeck

## Content of the analyses

### Data processing

1. Run fastqc on Ubuntu: Check the quality of the 32 files. Either do it using the interface or the fastqc.sh script.
2. Build your GRCh38 index genome
3. Align FASTQ files to GRCh38 : Use the staralign.sh bin bash script. You might have to sudo chmod +x staralign.sh before being able to run it. That script runs quantMode to get count tables.
4. The outputs will be
    - Out.tab
    - Readspergene.out.tab
    - Log.progress.out
    - Log.out
    - Log.final.out

Log.out: main log file with a lot of detailed information about the run. This file is most useful for troubleshooting and debugging. Log.progress.out: reports job progress statistics, such as the number of processed reads, % of mapped reads etc. It is updated in 1 minute intervals.
 Log.final.out: summary mapping statistics after mapping job is complete, very useful for quality control. The statistics are calculated for each read (single- or paired-end) and then summed or averaged over all reads. Note that STAR counts a paired-end read as one read, (unlike the samtools flagstat/idxstats, which count each mate separately). Most of the information is collected about the UNIQUE mappers (unlike samtools flagstat/idxstats which does not separate unique or multi-mappers). Each splicing is counted in the numbers of splices, which would correspond to summing the counts in SJ.out.tab. The mismatch/indel error rates are calculated on a per base basis, i.e. as total number of mismatches/indels in all unique mappers divided by the total number of mapped bases.

With --quantMode GeneCounts option STAR will count number reads per gene while mapping.
A read is counted if it overlaps (1nt or more) one and only one gene. Both ends of the pairedend read are checked for overlaps. The counts coincide with those produced by htseq-count with
default parameters. This option requires annotations (GTF or GFF with –sjdbGTFfile option) used
at the genome generation step, or at the mapping step. STAR outputs read counts per gene into

ReadsPerGene.out.tab file with 4 columns which correspond to different strandedness options:
column 1: gene ID
column 2: counts for unstranded RNA-seq -> The one we will use!
column 3: counts for the 1st read strand aligned with RNA (htseq-count option -s yes)
column 4: counts for the 2nd read strand aligned with RNA (htseq-count option -s reverse)

### Quality control

Run the R scripts:
- DESeq2_file_preparation: Formats your file to easily make a dds object
- DESeq2 : makes your dds object and the coldata + does the basic QC which I redid in the following report

### MCF10A_report_1

1. QC and batch effect analysis + correction:

  - size factors
  - dispersion estimates
  - PCA of unfiltered/filtered data
  - PCA of batch effect corrected data
  - distance matrix for sample comparison

2. DEG analyses

  - EGF vs IM
  - Oligomycin vs IM
  - Oligonycin vs EGF

For each of these:

  - Volcano plot figure (in figures)
  - GSEA analysis, databases queried:
    - Gene Ontology Biological process
    - Gene Ontology Molecular function
    - Reactome
    - KEGG
    - Wikipathways
- Excel document with all the differentially expressed genes (in data_output)
- Excel document with each GSEA result (in data_output) -> !! These ones list the genes found in each significant terms


# MCF10A_report_2


## Introduction
## PCA
## Heatmap of z score
## Comparing the EGF+Oligo vs all others (excluding GM)
### Volcano plots
### GSEA analysis (logFc)
## What are the contribution of EGF and Oligomycin to the PC1 dimesion
## GSEA
 - GOBP
 - Reactome

















