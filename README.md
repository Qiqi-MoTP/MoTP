# MoTP

The machine learning Multi-Omics Tumor Purity Prediction (MoTP) algorithm introduces a novel method for accurately estimating tumor purity from any single type of transcriptomic and DNA methylation data (mRNA, microRNA, Long non-coding RNA, or DNA methylation)

<img src="E:\luqiqi\tumor purity\手稿\FIG\FIG_1_page-0001.jpg" style="zoom:50%;" />



&nbsp;
&nbsp;

# Description

A novel algorithm that integrates mRNA expression, DNA methylation, miRNA expression, and lncRNA expression data within a uniform machine learning framework to predict tumor purity. Validated across TCGA pan-cancer datasets and various cancer cohorts, MoTP demonstrated superior performance compared to thirteen established algorithms, and higher accuracy than the algorithms based on a single omics data.



&nbsp;

# Details

+ The function `Preprocess()` is  preprocessing a list of omics data frames or matrices. Including impute missing values and [0,1] scaled (maximum=1, minimum=0). This step is not necessary, you can choose to process the data yourself, as long as the data has no missing values and the range is scaled to [0,1].
+ The function `MoTP()` is used to predicts tumor purity from either single-omics or multi-omics data using pre-trained models.
  + "data_list" is a collection of omics data frames or matrices. Specifically, each  element in the list represents a different omics dataset, It is essential to ensure that row names are features and col names are samples. The sample names must remain consistent across all the omics datasets provided. The list must contain at least one omics dataset.
  + "omics_list" is a character vector indicating the type of each omics data in data_list. Default is c("mRNA", "miRNA", "lncRNA", "DNA-methylation").

&nbsp;
&nbsp;

# Installation

- You can install the released version of **MoTP** with:
  &nbsp;

```R
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("WangX-Lab/MoTP")
```

&nbsp
&nbsp;

# Examples

&nbsp;
&nbsp;

## **Apply MoTP with single-omics data** 

### **Prepare data**

```R
library(ScPP)
load(system.file("data/binary.RData",package = "ScPP"))
```



**Bulk data**

```R
bulk[1:6,1:6]
```

|        | TCGA-CA-5256-01 | TCGA-AZ-6599-01 | TCGA-AA-3655-01 | TCGA-A6-6137-01 | TCGA-CK-4952-01 | TCGA-A6-5657-01 |
| ------ | --------------- | --------------- | --------------- | --------------- | --------------- | --------------- |
| HIF3A  | 3.7172          | 2.3437          | 2.0858          | 6.0759          | 1.9506          | 5.4777          |
| CAMK4  | 3.0698          | 4.9331          | 2.3709          | 4.1387          | 1.1557          | 4.1746          |
| RNF112 | 1.3702          | 2.4817          | 2.4947          | 3.5941          | 2.3486          | 4.9185          |
| SPN    | 5.5207          | 5.6704          | 6.8577          | 8.0598          | 5.0049          | 7.6076          |
| LRRTM1 | 3.2408          | 1.6031          | 0.9465          | 1.9142          | 0               | 3.2523          |
| GRIN1  | 3.0698          | 6.4944          | 4.3225          | 2.8073          | 7.346           | 4.5             |



**Binary data**

```R
head(binary)
```

|      | Sample          | Feature |
| ---- | --------------- | ------- |
| 1    | TCGA-CA-5256-01 | Tumor   |
| 2    | TCGA-AZ-6599-01 | Tumor   |
| 3    | TCGA-AA-3655-01 | Tumor   |
| 4    | TCGA-A6-6137-01 | Tumor   |
| 5    | TCGA-CK-4952-01 | Tumor   |
| 6    | TCGA-A6-5657-01 | Tumor   |



**scRNA-count**

```R
sc_count[1:6,1:6]
```

|          | SMC01.T_AAACCTGCATACGCCG | SMC01.T_AAACCTGGTCGCATAT | SMC01.T_AAACCTGTCCCTTGCA | SMC01.T_AAACGGGAGGGAAACA | SMC01.T_AAACGGGGTATAGGTA | SMC01.T_AAAGATGAGGCCGAAT |
| -------- | ------------------------ | ------------------------ | ------------------------ | ------------------------ | ------------------------ | ------------------------ |
| A1BG     | 0                        | 0                        | 0                        | 0                        | 0                        | 0                        |
| A1BG-AS1 | 0                        | 0                        | 0                        | 0                        | 0                        | 0                        |
| A1CF     | 0                        | 2                        | 0                        | 0                        | 3                        | 0                        |
| A2M      | 0                        | 0                        | 0                        | 0                        | 0                        | 0                        |
| A2M-AS1  | 0                        | 0                        | 0                        | 0                        | 0                        | 0                        |
| A2ML1    | 0                        | 0                        | 0                        | 0                        | 0                        | 0                        |


&nbsp;


### Execute ScPP to select the informative cells

```R
sc = sc_Preprocess(sc_count)
geneList = marker_Binary(bulk, binary, ref_group = "Normal")
res = ScPP(sc, geneList)
head(res$metadata)

# Phenotype+ genes
res$Genes_pos

# Phenotype- genes
res$Genes_neg

# Visualization of ScPP-identified cells
sc$ScPP = res$metadata$ScPP
Idents(sc) = "ScPP"
DimPlot(sc, group = "ScPP", cols = c("grey","blue","red"))

```

<img width="642" alt="image" src="https://github.com/WangX-Lab/ScPP/assets/54932820/1808015c-3790-4a07-ba15-411496d42d22">





&nbsp;
&nbsp;


## **Apply ScPP with continuous variables**

### Prepare data

```R
library(ScPP)
load(system.file("data/continuous.RData",package = "ScPP"))
```



**Bulk data**

```R
bulk[1:6,1:6]
```

|        | TCGA-AZ-6599-01 | TCGA-AA-3655-01 | TCGA-A6-6137-01 | TCGA-CK-4952-01 | TCGA-A6-5657-01 | TCGA-AD-6963-01 |
| ------ | --------------- | --------------- | --------------- | --------------- | --------------- | --------------- |
| HIF3A  | 2.3437          | 2.0858          | 6.0759          | 1.9506          | 5.4777          | 4.4634          |
| CAMK4  | 4.9331          | 2.3709          | 4.1387          | 1.1557          | 4.1746          | 3.2363          |
| RNF112 | 2.4817          | 2.4947          | 3.5941          | 2.3486          | 4.9185          | 1.4621          |
| SPN    | 5.6704          | 6.8577          | 8.0598          | 5.0049          | 7.6076          | 7.396           |
| LRRTM1 | 1.6031          | 0.9465          | 1.9142          | 0               | 3.2523          | 0               |
| GRIN1  | 6.4944          | 4.3225          | 2.8073          | 7.346           | 4.5             | 3.1816          |



**Continuous data**

```R
head(continuous)
```

|      | samp            | TMB_non_silent |
| ---- | --------------- | -------------- |
| 1    | TCGA-AZ-6599-01 | 178            |
| 2    | TCGA-AA-3655-01 | 65             |
| 3    | TCGA-A6-6137-01 | 91             |
| 4    | TCGA-CK-4952-01 | 206            |
| 5    | TCGA-A6-5657-01 | 63             |
| 6    | TCGA-AD-6963-01 | 67             |



**scRNA-count**

```R
sc_count[1:6,1:6]
```

|          | SMC01.T_AAACCTGCATACGCCG | SMC01.T_AAACCTGGTCGCATAT | SMC01.T_AAACCTGTCCCTTGCA | SMC01.T_AAACGGGAGGGAAACA | SMC01.T_AAACGGGGTATAGGTA | SMC01.T_AAAGATGAGGCCGAAT |
| -------- | ------------------------ | ------------------------ | ------------------------ | ------------------------ | ------------------------ | ------------------------ |
| A1BG     | 0                        | 0                        | 0                        | 0                        | 0                        | 0                        |
| A1BG-AS1 | 0                        | 0                        | 0                        | 0                        | 0                        | 0                        |
| A1CF     | 0                        | 2                        | 0                        | 0                        | 3                        | 0                        |
| A2M      | 0                        | 0                        | 0                        | 0                        | 0                        | 0                        |
| A2M-AS1  | 0                        | 0                        | 0                        | 0                        | 0                        | 0                        |
| A2ML1    | 0                        | 0                        | 0                        | 0                        | 0                        | 0                        |



&nbsp;

### Execute ScPP to select the informative cells

```R
sc = sc_Preprocess(sc_count)
geneList = marker_Continuous(bulk, continuous$TMB_non_silent)
res = ScPP(sc, geneList)
head(res$metadata)

# Phenotype+ genes
res$Genes_pos

# Phenotype- genes
res$Genes_neg

# Visualization of ScPP-identified cells
sc$ScPP = res$metadata$ScPP
Idents(sc) = "ScPP"
DimPlot(sc, group = "ScPP", cols = c("grey","blue","red"))
```

<img width="642" alt="image" src="https://github.com/WangX-Lab/ScPP/assets/54932820/6c7dac7b-138b-41e2-9ef1-af87baff06e3">



&nbsp;
&nbsp;

## **Apply ScPP with  survival data**

### Prepare data

```R
library(ScPP)
load(system.file("data/survival.RData",package = "ScPP"))
```



**Bulk data**

```R
bulk[1:6,1:6]
```

|         | TCGA-69-7978 | TCGA-62-8399 | TCGA-78-7539 | TCGA-73-4658 | TCGA-44-6775 | TCGA-44-2655 |
| ------- | ------------ | ------------ | ------------ | ------------ | ------------ | ------------ |
| HIF3A   | 4.2598       | 11.6239      | 9.1362       | 5.0288       | 4.0573       | 5.5335       |
| RTN4RL2 | 8.2023       | 5.5819       | 3.5365       | 7.4156       | 7.7107       | 5.3257       |
| HMGCLL1 | 2.7476       | 5.8513       | 3.8334       | 3.6447       | 2.9188       | 4.882        |
| LRRTM1  | 0            | 0.4628       | 4.7506       | 6.8005       | 7.7819       | 2.2882       |
| GRIN1   | 6.6074       | 5.4257       | 4.9563       | 7.351        | 3.5361       | 3.3311       |
| LRRTM3  | 1.7458       | 2.0092       | 0            | 1.4468       | 0            | 0            |



**Survival data**

```R
head(survival)
```

|              | status | time  |
| ------------ | ------ | ----- |
| TCGA-69-7978 | 0      | 4.4   |
| TCGA-62-8399 | 0      | 88.57 |
| TCGA-78-7539 | 0      | 25.99 |
| TCGA-73-4658 | 1      | 52.56 |
| TCGA-44-6775 | 0      | 23.16 |
| TCGA-44-2655 | 0      | 43.5  |



**scRNA-count**

```R
sc_count[1:6,1:6]
```

|          | SMC01.T_AAACCTGCATACGCCG | SMC01.T_AAACCTGGTCGCATAT | SMC01.T_AAACCTGTCCCTTGCA | SMC01.T_AAACGGGAGGGAAACA | SMC01.T_AAACGGGGTATAGGTA | SMC01.T_AAAGATGAGGCCGAAT |
| -------- | ------------------------ | ------------------------ | ------------------------ | ------------------------ | ------------------------ | ------------------------ |
| A1BG     | 0                        | 0                        | 0                        | 0                        | 0                        | 0                        |
| A1BG-AS1 | 0                        | 0                        | 0                        | 0                        | 0                        | 0                        |
| A1CF     | 0                        | 2                        | 0                        | 0                        | 3                        | 0                        |
| A2M      | 0                        | 0                        | 0                        | 0                        | 0                        | 0                        |
| A2M-AS1  | 0                        | 0                        | 0                        | 0                        | 0                        | 0                        |
| A2ML1    | 0                        | 0                        | 0                        | 0                        | 0                        | 0                        |


&nbsp;


### Execute ScPP to select the informative cells

```R
sc = sc_Preprocess(sc_count)
geneList = marker_Survival(bulk, survival)
res = ScPP(sc, geneList)
head(res$metadata)

# Phenotype+ genes
res$Genes_pos

# Phenotype- genes
res$Genes_neg

# Visualization of ScPP-identified cells
sc$ScPP = res$metadata$ScPP
Idents(sc) = "ScPP"
DimPlot(sc, group = "ScPP", cols = c("grey","blue","red"))


```

<img width="642" alt="image" src="https://github.com/WangX-Lab/ScPP/assets/54932820/47404d94-abe4-485c-8657-0f5e47bc62c3">


&nbsp;
# Contact

E-mail any questions to Xiaosheng Wang (xiaosheng.wang@cpu.edu.cn)
