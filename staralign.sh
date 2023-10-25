
#!/bin/bash

STAR_PATH="/usr/bin/STAR"

# Loop through the samples
for sample_number in {15..16}
do
    # Set timer to 0
    SECONDS=0
    # Generate sample names for both A and Z
    SAMPLE_A="${sample_number}A"
    SAMPLE_Z="${sample_number}Z"

    # Input fastq files for R1 and R2
    R1_FILE_A=./mcf10a/raw_data/MCF10A_BQ${SAMPLE_A}_L004_R1_001.fastq.gz
    R2_FILE_A=./mcf10a/raw_data/MCF10A_BQ${SAMPLE_A}_L004_R2_001.fastq.gz
    R1_FILE_Z=./mcf10a/raw_data/MCF10A_BQ${SAMPLE_Z}_L004_R1_001.fastq.gz
    R2_FILE_Z=./mcf10a/raw_data/MCF10A_BQ${SAMPLE_Z}_L004_R2_001.fastq.gz

    # Output directories for aligned files
    OUTPUT_DIR_A=./mcf10a/data/aligned_sample_${SAMPLE_A}
    OUTPUT_DIR_Z=./mcf10a/data/aligned_sample_${SAMPLE_Z}

    # Create the output directories if they don't exist
    mkdir -p "$OUTPUT_DIR_A"
    mkdir -p "$OUTPUT_DIR_Z"

    # Run STAR alignment for samples A
    $STAR_PATH \
        --runThreadN 15 \
        --readFilesCommand zcat \
        --genomeDir ./genome_ref/human_index \
        --readFilesIn $R1_FILE_A $R2_FILE_A \
        --outFileNamePrefix aligned_sample_${SAMPLE_A} \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode GeneCounts


    # Run STAR alignment for samples Z
    $STAR_PATH \
        --runThreadN 15 \
        --readFilesCommand zcat \
        --genomeDir ./genome_ref/human_index \
        --readFilesIn $R1_FILE_Z $R2_FILE_Z \
        --outFileNamePrefix aligned_sample_${SAMPLE_Z} \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode GeneCounts

    # Print status message
    echo "Alignments for Samples ${SAMPLE_A} and ${SAMPLE_Z} completed."
    # Get time elapsed
    ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
done


