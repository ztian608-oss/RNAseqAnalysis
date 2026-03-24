# DEanalysis

## Description
`DEanalysis` 提供一个基于 **DESeq2** 的 RNA-seq 差异表达分析流程，支持：
- HTSeq count 导入 + GTF 注释映射
- DESeq2 拟合 + VST 标准化 + PCA
- pairwise/custom contrasts
- volcano 图、显著基因热图
- 多组（时间序列）k-means 表达模式聚类

---

## Project structure

- `deseq2_generic.R`：**主入口脚本**（建议统一使用）
- `R/utils_cli.R`：CLI 参数解析、帮助信息、基础校验
- `R/io.R`：输入读取与检查（sample table / gtf）
- `R/analysis.R`：DESeq2 建模、比较、显著基因提取
- `R/plots.R`：PCA 与 volcano
- `R/clustering.R`：热图与时间序列聚类
- `templates/run_deseq2_generic.sh`：可复用命令模板

---

## Input requirements

### 1) Count directory (`--count_dir`)
HTSeq count 文件所在目录。

### 2) Sample table (`--sample_table`)
CSV，必须包含：
- `sampleName`
- `sampleFile`
- `group_name`

可选：
- `group_order`（控制组别顺序，影响绘图与比较顺序）

### 3) Gene annotation (`--gtf`)
GTF/GFF，至少包含 `gene_id` 和 `gene_name`（建议含 `gene_type`）。

### 4) Reference group (`--ref_group`)
pairwise 比较基准组。

---

## CLI quick start

查看参数帮助：

```bash
Rscript deseq2_generic.R --help
```

运行模板（先修改路径；默认是 pairwise 模式）：

```bash
bash templates/run_deseq2_generic.sh
```

---

## Step-by-step commands

> 下面示例把流程拆成“实际运行步骤”，便于你后续封装 skill。

### Step 0. 准备输入文件

- `sample_table.csv`
- `annotation.gtf`
- `counts/` 目录
- （可选）`contrasts.csv`，列名为 `treat,base`

---

### Step 1. 先跑基础流程（PCA + pairwise DE）

```bash
Rscript deseq2_generic.R \
  --count_dir=counts \
  --sample_table=sample_table.csv \
  --gtf=annotation.gtf \
  --ref_group=control \
  --outdir=result_step1 \
  --comparison_mode=pairwise \
  --run_volcano=FALSE \
  --run_heatmap=FALSE
```

主要输出：
- `result_step1/vsd.csv`
- `result_step1/PCA_DESeq2.pdf`
- `result_step1/PCA_DESeq2_data.csv`
- `result_step1/*_DE.csv`

---

### Step 2. 打开 volcano 与显著基因热图

```bash
Rscript deseq2_generic.R \
  --count_dir=counts \
  --sample_table=sample_table.csv \
  --gtf=annotation.gtf \
  --ref_group=control \
  --outdir=result_step2 \
  --comparison_mode=pairwise \
  --padj_cutoff=0.01 \
  --log2_cutoff=1 \
  --run_volcano=TRUE \
  --run_heatmap=TRUE
```

新增输出：
- `result_step2/volcanoplot_*.pdf`
- `result_step2/significant_de_genes.txt`
- `result_step2/significant_de_heatmap.pdf`

---

### Step 3. 增加自定义 contrasts（treat vs base，可选）

```bash
Rscript deseq2_generic.R \
  --count_dir=counts \
  --sample_table=sample_table.csv \
  --gtf=annotation.gtf \
  --ref_group=control \
  --contrast_file=contrasts.csv \
  --comparison_mode=both \
  --outdir=result_step3 \
  --run_volcano=TRUE \
  --run_heatmap=TRUE
```

> 注意：`--comparison_mode=both/custom` 时，必须提供有效 `--contrast_file`。

新增输出（每个比较）：
- `*_DE.csv`
- `*_summary.csv`

---

### Step 4. 多组/时间序列聚类（group 数 > 2）

```bash
Rscript deseq2_generic.R \
  --count_dir=counts \
  --sample_table=sample_table.csv \
  --gtf=annotation.gtf \
  --ref_group=T0 \
  --outdir=result_step4 \
  --comparison_mode=both \
  --cluster_num=6 \
  --seed=12345 \
  --run_volcano=TRUE \
  --run_heatmap=TRUE
```

聚类输出：
- `timecourse_cluster_heatmap.pdf`
- `timecourse_cluster_result.csv`
- `timecourse_cluster_patterns.pdf`

---

## PCA variance information

PCA 输出已包含方差信息：
- X/Y 轴标签显示 PC1/PC2 解释方差
- subtitle 显示 `PC1 + PC2` 累积解释方差
- `PCA_DESeq2_data.csv` 中包含：
  - `PC1_var_percent`
  - `PC2_var_percent`
  - `PC1_PC2_total_percent`

---


## Skill 接口 JSON（推荐）

为了后续封装 skill，建议补充输入/输出 JSON schema。仓库里已提供两个模板：

- `templates/skill_input.schema.json`：定义运行参数接口（路径、阈值、模式开关等）
- `templates/skill_output.schema.json`：定义产物输出端口（文件路径、比较列表、统计摘要）

这样做的好处：
- skill 编排时可以自动校验参数
- 输出格式稳定，便于后续节点消费
- 更容易做 API 化和批处理

---

## Main options

- `--design='~ group_name'`
- `--contrast_file=contrasts.csv`
- `--comparison_mode=pairwise|custom|both`
- `--padj_cutoff=0.01`
- `--log2_cutoff=1`
- `--min_count_mean=10`
- `--cluster_num=6`
- `--seed=12345`
- `--run_volcano=TRUE|FALSE`
- `--run_heatmap=TRUE|FALSE`

---

## Notes

- 脚本会先校验 sample table 和 count 文件存在性。
- 注释映射或表达过滤后无基因时，会直接报错退出。
- 无显著基因时，热图和聚类会自动跳过。
