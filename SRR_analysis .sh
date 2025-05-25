#!/bin/bash

# 1. 设置变量
SRA_ID="SRR5354177"
REF_GENOME="reference_genome.fa"
ADAPTERS="/home/liyuxuan5354/.conda/envs/omics_analysis/share/trimmomatic-0.39-2/adapters/TruSeq3-PE.fa"
ANNOTATION="/home/liyuxuan5354/annotation.gtf"
THREADS=8

# 2. 日志记录
exec > >(tee -a "${SRA_ID}_analysis.log") 2>&1

# 3. 下载并转换SRA数据
echo "1. 下载并转换SRA数据..."
if ! prefetch "$SRA_ID"; then
    echo "Error: Failed to download SRA data."
    exit 1
fi
fasterq-dump --split-files --threads "$THREADS" "./$SRA_ID/$SRA_ID.sra"

# 4. 数据质控
echo "2. 运行FastQC进行数据质控..."
fastqc "${SRA_ID}_1.fastq" "${SRA_ID}_2.fastq"

# 5. 数据过滤
echo "3. 运行Trimmomatic进行数据过滤..."
trimmomatic PE -phred33 \
    "${SRA_ID}_1.fastq" "${SRA_ID}_2.fastq" \
    "${SRA_ID}_1_trimmed.fastq" "${SRA_ID}_1_unpaired.fastq" \
    "${SRA_ID}_2_trimmed.fastq" "${SRA_ID}_2_unpaired.fastq" \
    ILLUMINACLIP:"$ADAPTERS":2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# 6. 下载参考基因组（如果未下载）
if [ ! -f "$REF_GENOME" ]; then
    echo "4. 下载参考基因组..."
    wget -O "${REF_GENOME}.gz" \
        "https://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
    gunzip "${REF_GENOME}.gz"
fi

# 7. 构建HISAT2索引（如果未构建）
if [ ! -f "${REF_GENOME}.1.ht2" ]; then
    echo "5. 构建HISAT2索引..."
    hisat2-build "$REF_GENOME" "$REF_GENOME" || {
        echo "Error: Failed to build HISAT2 index."
        exit 1
    }
fi

# 8. 比对到参考基因组
echo "6. 比对到参考基因组..."
if ! hisat2 -x "$REF_GENOME" \
    -1 "${SRA_ID}_1_trimmed.fastq" -2 "${SRA_ID}_2_trimmed.fastq" \
    -S "${SRA_ID}_aligned.sam" --threads "$THREADS"; then
    echo "Error: Failed to align reads to the reference genome."
    exit 1
fi

# 9. 生成BAM文件并排序
echo "7. 生成BAM文件并排序..."
samtools view -bS "${SRA_ID}_aligned.sam" > "${SRA_ID}_aligned.bam" || {
    echo "Error: Failed to convert SAM to BAM."
    exit 1
}
samtools sort "${SRA_ID}_aligned.bam" -o "${SRA_ID}_aligned_sorted.bam" || {
    echo "Error: Failed to sort BAM file."
    exit 1
}

# 10. 基因定量
echo "8. 运行featureCounts进行基因定量..."
featureCounts -a "$ANNOTATION" -o "${SRA_ID}_counts.txt" -p "${SRA_ID}_aligned_sorted.bam" || {
    echo "Error: Failed to run featureCounts."
    exit 1
}

# 11. 统计比对率和定量分配比例
echo "9. 统计比对率和定量分配比例..."
alignment_rate=$(grep "overall alignment rate" "${SRA_ID}_aligned.log" | awk '{print $1}')
assigned_ratio=$(grep "Successfully assigned alignments" "${SRA_ID}_counts.txt.summary" | awk '{print $2}')

echo "比对率: $alignment_rate"
echo "定量分配比例: $assigned_ratio"

echo "组学数据分析流程已完成！"
