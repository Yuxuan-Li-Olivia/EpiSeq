# EpiSeq
EpiSeq: Genomic Data Analysis for Epilepsy Genetics
EpiSeq: 
癫痫遗传学的组学数据分析
 
作者：李昱萱
 523111910123
 1.
摘要
 
本项目通过NCBI SRA数据库获取组学原始数据，利用conda搭建计算环境，并使用Bash脚
本编写分析流程，完成了从数据下载、质控、过滤、比对到基因定量的全流程分析。结果显
示，整体比对率为86.30%，基因定量成功分配比例达58.50%。本文详细介绍了分析流程、
工具使用及结果。
2.
前言
 
2.1
背景介绍
 
随着高通量测序技术的发展，组学数据分析在生物医学研究中扮演着越来越重要的角色。本
项目旨在通过NCBI SRA数据库获取公开的组学数据，搭建计算环境并编写分析流程，完成
从原始数据到基因定量的全流程分析。我本次使用的数据来源于NCBI SRA数据库，数据集
编号为 SRR5354177，属于 Epi25 
联盟
 的研究项目。该联盟通过大规模协作，致力于揭示
癫痫的复杂遗传成分，整合了包括 Epi4K、EPIGEN、EuroEPINOMICS 等在内的多个全球
研究力量。作为 NHGRI 
常见病基因组学中心
 的一部分，该数据集为癫痫遗传学研究提供
了重要的基因组资源，有助于深入理解疾病机制并推动精准医学的发展。
2.2
工作思路
 
1. 从NCBI SRA数据库下载数据集（SRR5354177）。
2. 使用conda搭建计算环境并安装必要的工具。
3. 编写Bash脚本实现从数据下载、质控、过滤、比对到基因定量的全流程分析。
2.3
工作概要
 
数据来源：NCBI SRA数据库（SRR5354177）。
工具：SRA Toolkit、FastQC、Trimmomatic、HISAT2、samtools、
featureCounts。
流程：数据下载 → 数据质控 → 数据过滤 → 比对到参考基因组 → 基因定量。
3.
数据集与方法
 
3.1
数据来源
 
数据集编号：SRR5354177
实验编号 (SRX): SRX3900427
实验描述：NHGRI CCDG Epi25 Illumina random exon sequencing of 'genomic 
DNA' paired-end library 'NexPond-579960' containing sample 'HC009-1' from 
subject 'PT-1HS7L'
测序平台：Illumina HiSeq X Ten
数据量：1 run, 39.5M spots, 11.9G bases, 2.6Gb downloads
物种：Homo sapiens（智人）
3.2
工具与方法
 
1. 
数据下载与转换
  使用 
prefetch 下载 SRA 数据（
SRR5354177）。使用 
fasterq-dump 将 
.sra 文件转换为 
.fastq 格式。
prefetch SRR5354177
 fasterq-dump --split-files --threads 8 
./SRR5354177/SRR5354177.sra
 2. 
数据质控
使用FastQC对原始数据进行质量评估。
fastqc SRR5354177_1.fastq SRR5354177_2.fastq
 3. 
数据过滤
使用Trimmomatic过滤低质量序列和接头序列，保留高质量读长对。
trimmomatic PE -phred33 \
 SRR5354177_1.fastq SRR5354177_2.fastq \
 SRR5354177_1_trimmed.fastq SRR5354177_1_unpaired.fastq \
 SRR5354177_2_trimmed.fastq SRR5354177_2_unpaired.fastq \
 ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 
\
 SLIDINGWINDOW:4:15 MINLEN:36
 4. 
参考基因组下载与索引构建
下载人类参考基因组（
Homo_sapiens.GRCh38）。
使用 HISAT2 构建参考基因组索引。
5. 
比对到参考基因组
使用HISAT2将过滤后的数据比对到参考基因组，生成
hisat2 -x reference_genome \-1 SRR5354177_1_trimmed.fastq -2 
SRR5354177_2_trimmed.fastq \-S SRR5354177.sam
 .sam文件。
6. 
生成
BAM
文件并排序
使用samtools将SAM文件转换为BAM文件并进行排序。
samtools view -bS SRR5354177.sam > SRR5354177.bam
 samtools sort SRR5354177.bam -o SRR5354177_sorted.bam
 7. 
基因定量
使用featureCounts进行基因定量，生成基因表达计数文件。
featureCounts -a annotation.gtf -o counts.txt 
SRR5354177_sorted.bam
 4.
结果
 
1. 
数据质控
  FastQC报告显示，原始数据质量较高，但需过滤低质量序列和接头序
列。
2. 
数据过滤
  Trimmomatic过滤后，95.20%的读长对保留下来，过滤效果良好。
3. 
比对结果
  HISAT2比对结果显示，整体比对率为68.30%。
4. 
基因定量
  featureCounts定量结果显示，成功分配比例达58.50%。
5.
改进与讨论
 
1. 
数据质控与过滤
FastQC结果显示原始数据质量较高，但存在一定比例的低质量序列和接头序列。
通过Trimmomatic过滤后，95.20%的读长对得以保留，表明过滤策略有效，能够
为后续分析提供高质量数据。
2. 
比对结果
HISAT2比对结果显示，整体比对率为68.30%。这一比例略低于预期，可能的原因
包括：
参考基因组版本与数据来源不完全匹配。
数据中存在未注释的基因组区域或复杂结构变异。
过滤后的读长对中仍包含部分低质量序列。
为提高比对率，可以考虑：
使用更全面的参考基因组版本（如包含非编码区域）。
调整HISAT2的比对参数，增加对复杂区域的容忍度。
结合其他比对工具（如STAR）进行验证。
3. 
基因定量
featureCounts定量结果显示，成功分配比例为58.50%。这一比例表明，部分读长
未被分配到已知基因区域，可能原因包括：
读长比对到非编码区域或未注释的基因组区域。
基因注释文件（GTF）不完整或与参考基因组版本不匹配。
数据中存在转录本异构体或剪接变异。
6.
参考文献
 
[1]  FAURE J, MORIN G, DUTERTRE F. Incidences neuropsychiatriques et endocriniennes 
dans l'épilepsie temporale [Neuropsychiatric and endocrine aspects in temporal epilepsy]. 
Rev Neurol (Paris). 1953;88(5):382-4. Undetermined Language. PMID: 13121714.
 [2]  Liao Y, Smyth GK, Shi W. featureCounts: an efficient general purpose program for 
assigning sequence reads to genomic features. Bioinformatics. 2014 Apr 1;30(7):923-30. 
doi: 10.1093/bioinformatics/btt656. Epub 2013 Nov 13. PMID: 24227677.
 [3]  Pertea, M., Kim, D., Pertea, G. et al. Transcript-level expression analysis of RNA-seq 
experiments with HISAT, StringTie and Ballgown. Nat Protoc 11, 1650–1667 (2016). http
 s://doi.org/10.1038/nprot.2016.095
7.
附录
 
计算环境配置信息
 
操作系统：Ubuntu 20.04 LTS  
Conda环境：omics_analysis  
工具版本：  
SRA Toolkit: 3.0.0  
FastQC: 0.11.9  
Trimmomatic: 0.39  
HISAT2: 2.2.1  
samtools: 1.12  
featureCounts: 2.0.1  
8.
项目核心脚本
 
完整脚本及分析产生的相关文件已发布至GitHub：[项目链接](https://github.com/Yuxuan
Li-Olivia/EpiSeq.gi
