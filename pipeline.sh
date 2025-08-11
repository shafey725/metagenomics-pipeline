#!/bin/bash
set -e

# ===== Step 0: Find FragGeneScan =====
echo "Searching for FragGeneScan binary..."
FGS_PATH=$(find /home/shafey/f1 -type f -name FragGeneScan 2>/dev/null | head -n 1)

if [ -z "$FGS_PATH" ]; then
    echo "ERROR: FragGeneScan binary not found under /home/shafey/f1/"
    echo "Please ensure FragGeneScan is downloaded and compiled."
    exit 1
fi

chmod +x "$FGS_PATH"
echo "Found FragGeneScan at: $FGS_PATH"

# ===== Step 1: Download required reference databases =====
mkdir -p ref-dbs

# SILVA SSU
if [ ! -d ref-dbs/silva_ssu-20200130 ]; then
    echo "Downloading SILVA SSU..."
    wget -P ref-dbs https://ftp.ebi.ac.uk/pub/databases/metagenomics/pipeline-5.0/ref-dbs/silva_ssu-20200130.tar.gz
    echo "Extracting SILVA SSU..."
    tar -xzf ref-dbs/silva_ssu-20200130.tar.gz -C ref-dbs
else
    echo "SILVA SSU already exists, skipping download."
fi

# SILVA LSU
if [ ! -d ref-dbs/silva_lsu-20200130 ]; then
    echo "Downloading SILVA LSU..."
    wget -P ref-dbs https://ftp.ebi.ac.uk/pub/databases/metagenomics/pipeline-5.0/ref-dbs/silva_lsu-20200130.tar.gz
    echo "Extracting SILVA LSU..."
    tar -xzf ref-dbs/silva_lsu-20200130.tar.gz -C ref-dbs
else
    echo "SILVA LSU already exists, skipping download."
fi

# Rfam models
mkdir -p rfam_models/{other_models,ribosomal_models}

if [ -z "$(ls -A rfam_models/other_models)" ]; then
    echo "Downloading other_models..."
    wget -r -nH --cut-dirs=7 --no-parent \
        ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/pipeline-5.0/ref-dbs/rfam_models/other_models/ \
        -P rfam_models/other_models
else
    echo "other_models already exists — skipping download."
fi

if [ -z "$(ls -A rfam_models/ribosomal_models)" ]; then
    echo "Downloading ribosomal_models..."
    wget -r -nH --cut-dirs=7 --no-parent \
        ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/pipeline-5.0/ref-dbs/rfam_models/ribosomal_models/ \
        -P rfam_models/ribosomal_models
else
    echo "ribosomal_models already exists — skipping download."
fi

# ===== Step 2: Merge paired-end reads =====
for R1 in *_1.fastq; do
    R2="${R1/_1.fastq/_2.fastq}"
    sample="${R1%_1.fastq}"

    SeqPrep -f "$R1" -r "$R2" \
        -1 "${sample}_unmapped_R1.fastq" \
        -2 "${sample}_unmapped_R2.fastq" \
        -s "${sample}_merged.fastq.gz"

    echo "Finished merging $sample"
done

# ===== Step 3: Trim merged reads =====
mkdir -p trimmed_reads

for merged in *_merged.fastq.gz; do
    sample="${merged%_merged.fastq.gz}"

    trimmomatic SE -phred33 "$merged" "${sample}_trimmed.fq.gz" \
        ILLUMINACLIP:adaptor.fa:2:30:10 \
        LEADING:3 TRAILING:3 \
        SLIDINGWINDOW:4:15 MINLEN:36

    mkdir -p trimmed_reads/"$sample"
    mv "${sample}_trimmed.fq.gz" trimmed_reads/"$sample"/

    echo "Finished trimming $sample"
done

# ===== Step 4–9: Process all samples =====
for sample_dir in trimmed_reads/*; do
    sample=$(basename "$sample_dir")
    echo "Processing sample: $sample"

    cd "$sample_dir"

    # Step 4: FastQC
    fastqc "${sample}_trimmed.fq.gz"

    # Step 5: Convert to FASTA
    seqkit fq2fa "${sample}_trimmed.fq.gz" > "${sample}.fa"

    # Step 6: cmsearch
    cmsearch --tblout "${sample}_ribo.tbl" ../../rfam_models/ribosomal_models/ribo.cm "${sample}.fa"

    # Step 7: Mask ribosomal regions
    awk '/^[^#]/ {if ($8 > $9) {print $1 "\t" $9-1 "\t" $8} else {print $1 "\t" $8-1 "\t" $9}}' \
        "${sample}_ribo.tbl" > "${sample}_ribo.bed"

    bedtools maskfasta -fi "${sample}.fa" -bed "${sample}_ribo.bed" -fo "${sample}_masked.fasta"

    # Step 8: Extract noncoding regions
    bedtools getfasta -fi "${sample}.fa" -bed "${sample}_ribo.bed" -fo "${sample}_noncoding_output.fasta"

    # Step 9: FragGeneScan
    "$FGS_PATH" \
        -s "${sample}_masked.fasta" \
        -o "${sample}_cds.fa" \
        -w 0 \
        -t complete

    # Step 10: MAPseq SSU
    mapseq "${sample}_noncoding_output.fasta" \
        ../../ref-dbs/silva_ssu-20200130/SSU.fasta \
        ../../ref-dbs/silva_ssu-20200130/slv_ssu_filtered2.txt \
        > mapseq_ssu_output.txt

    mapseq -otucounts mapseq_ssu_output.txt > "${sample}_ssu.otu"

    awk 'NR>1 {gsub(/;/, "\t", $3); print $4 "\t" $3}' "${sample}_ssu.otu" > "${sample}_ssu.txt"

    ktImportText "${sample}_ssu.txt" -o "${sample}_ssu_output.html"

    # Step 11: MAPseq LSU
    mapseq "${sample}_noncoding_output.fasta" \
        ../../ref-dbs/silva_lsu-20200130/LSU.fasta \
        ../../ref-dbs/silva_lsu-20200130/slv_lsu_filtered2.txt \
        > mapseq_lsu_output.txt

    mapseq -otucounts mapseq_lsu_output.txt > "${sample}_lsu.otu"

    awk 'NR>1 {gsub(/;/, "\t", $3); print $4 "\t" $3}' "${sample}_lsu.otu" > "${sample}_lsu.txt"

    ktImportText "${sample}_lsu.txt" -o "${sample}_lsu_output.html"

    cd ../..
done

echo "Pipeline completed for all samples."

