# README

## 上传文件说明

- 差异基因分析脚本 diff_exp_analysis.R
- 分析脚本 diff_exp_analysis.sh
- Dockerfile
- README
- 项目报告
- SRR-Acc-List.txt

## 构建镜像并使用脚本

进入Dockerfile目录下，使用如下命令构建镜像

```
docker build -t your-image-name .
```

你可以运行如下命令来使用这个镜像

```
docker run -it --workdir / -v /host/data:/data your-image-name bash
```

这条命令将会挂载宿主机的`/host/data`目录到容器的`/data`目录，你可以在容器内通过`/data`

访问到宿主机的`/host/data`目录的内容。并且在你登入这个容器后，你会位于`/`根目录下。

```
cd scripts
```

这个目录下存放了差异表达分析使用的shell脚本，你可以通过如下命令使用这个脚本

```
bash diff_exp_analysis.sh /data
```

注意，你需要将工作目录挂载到容器的`/data`目录下

## 工作目录结构

- 参考基因组文件fna和注释文件gff

- 多个下载好的SRR文件夹

  - SRA文件

  - FILENAME_QC_report_results

  - FILENAME_QC_results

  - FILENAME_mapping_results

- mapping_results