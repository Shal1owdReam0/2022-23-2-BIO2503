FROM continuumio/anaconda3
LABEL maintainer="2022-2023-2_BIO2503_Group9"

RUN conda config --add channels bioconda

RUN apt-get update && apt-get install -y \
    samtools \
    r-base

WORKDIR /tools

# 安装sratoolkit
RUN wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz \
    && tar -zxvf sratoolkit.current-ubuntu64.tar.gz \
    && rm sratoolkit.current-ubuntu64.tar.gz

ENV PATH="/tools/sratoolkit.3.0.5-ubuntu64/bin:${PATH}"

# 安装fastqc
RUN conda install -y fastqc

# 安装trimmomatic
RUN conda install -y trimmomatic

# 安装hisat2 
RUN conda install -y hisat2

# 安装subread以使用featureCounts
RUN wget https://nchc.dl.sourceforge.net/project/subread/subread-2.0.3/subread-2.0.3-Linux-x86_64.tar.gz \ 
    && tar xzf subread-2.0.3-Linux-x86_64.tar.gz \
    && rm subread-2.0.3-Linux-x86_64.tar.gz

ENV PATH="/tools/subread-2.0.3-Linux-x86_64/bin:${PATH}"

COPY ./diff_exp_analysis.R /scripts/diff_exp_analysis.R
COPY ./diff_exp_analysis.sh /scripts/diff_exp_analysis.sh

RUN chmod +x /scripts/diff_exp_analysis.R
RUN chmod +x /scripts/diff_exp_analysis.sh

# 安装所需的R包和依赖
RUN R -e "install.packages(c('DESeq2', 'dplyr', 'ggpubr', 'pheatmap', 'ggplot2'), repos='http://cran.us.r-project.org')"

# 运行容器时，默认执行的命令
CMD ["/bin/bash"]