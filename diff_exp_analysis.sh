#!/bin/bash

# 指定工作目录
WORK_DIR=$1
cd $WORK_DIR

# 工作目录结构
# - 参考基因组文件fna和注释文件gff
# - reference_index（运行中创建）
# - 多个下载好的SRR文件夹
#   - SRA文件
#   - FILENAME_QC_report_results(运行中创建)
#   - FILENAME_QC_results（运行中创建）
#   - FILENAME_mapping_results（运行中创建）
# - mapping_results（运行中创建）

mkdir -p mapping_results

# 指定参考基因组文件和注释文件
GENOME=GCF_000001635.27_GRCm39_genomic.fna
GFF=GCF_000001635.27_GRCm39_genomic.gff

# 使用参考基因组建立HISAT2索引
if [ ! -e "reference_index.1.ht2" ]; then
    echo "Index not found, creating..."
    hisat2-build $GENOME reference_index
else
    echo "Index found, skipping index building step."
fi

# 对每个sra文件进行质量控制分析和修剪
for  FOLDERNAME in SRR*/
do  
    FILENAME="${FOLDERNAME%/}"
    cd $FILENAME

    # 创建文件夹,储存流程中得到的一些文件

    mkdir -p "${FILENAME}_QC_report_results"
    mkdir -p "${FILENAME}_QC_results"
    mkdir -p "${FILENAME}_mapping_results"


    # sra文件转换为fastq文件
    fastq-dump --split-3  "${FILENAME}.sra"
    
    # fastqc质控报告
    # fastqc  -t 12 -o "${FILENAME}_QC_report_results" $Seq1 $Seq2

    # Trimmomatic修剪，这里假定为双端测序
    trimmomatic PE -phred33 "${FILENAME}_1.fastq" "${FILENAME}_2.fastq" \
    "${FILENAME}_QC_results/${FILENAME}_1.trim.fastq" "${FILENAME}_QC_results/${FILENAME}_1.unpaired.trim.fastq" \
    "${FILENAME}_QC_results/${FILENAME}_2.trim.fastq" "${FILENAME}_QC_results/${FILENAME}_2.unpaired.trim.fastq" \
    HEADCROP:15 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

    # 使用HISAT2比对到参考基因组，得到sam文件
    hisat2 -p 8 -x "../reference_index" -1 "${FILENAME}_QC_results/${FILENAME}_1.trim.fastq" -2 "${FILENAME}_QC_results/${FILENAME}_2.trim.fastq" \
    -S "${FILENAME}_mapping_results/${FILENAME}_mapping.sam"


    cd "${FILENAME}_mapping_results"
    samtools view -bS "${FILENAME}_mapping.sam" > "${FILENAME}_mapping.bam"
    samtools sort "${FILENAME}_mapping.bam" -o "${FILENAME}_mapping_sorted.bam"
    mv "${FILENAME}_mapping_sorted.bam" ../../mapping_results

    # 回到主工作目录
    cd ../
    cd ../
done

# 使用featureCounts计算基因表达矩阵
cd $mapping_results
featureCounts -a $GFF -p -g gene -o Martrix.txt /*_mapping_sorted.bam

# 使用R语言脚本进行差异表达分析
Rscript diff_exp_analysis.R