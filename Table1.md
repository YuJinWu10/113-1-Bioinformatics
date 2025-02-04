---
title: "113-1 Bioinformatics table1" ## 標題
author: "Yu-Jin Wu"                  ## 作者
date: "`r Sys.Date()`"               ## 日期(Sys.Date() 會自動插入當前日期)
output:                              ## 輸出格式
  html_document:
    code_folding: hide               ## code_folding: hide：可折疊程式碼區塊
    toc: yes                         ## toc: yes：加入目錄
    toc_depth: 2                     ## toc_depth: 2：目錄顯示到第二層
    toc_float:
      collapsed: no                  ## toc_float:     collapsed: no：目錄不預設折疊
      smooth_scroll: no              ##                smooth_scroll: no：關閉滑動效果
  pdf_document:
    toc: yes
    toc_depth: '2'
  word_document:
    toc: yes
    toc_depth: '2'

---

# Install and Library the required packages
## Install the required packages
```{r}
# install.packages(c("dplyr", "table1", "Hmisc", "Exact", "DescTools"))
```
- 下載所需的package

## Library the required packages
```{r}
library(dplyr) # 資料處理與轉換
library(table1) # 製作 Table 1
library(Hmisc) # 醫學統計分析（含 upData()）
library(Exact) # Fisher’s exact test（適用於小樣本的列聯表檢定）
library(DescTools) # 統計工具（包含 GTest()，用於卡方檢定）
```


## Import dataset 
```{r}
# (Method 1)
setwd("D:/yuyu/master/113-1/TA/113-1 Bioinformatics/Final")  # setwd()：設定工作目錄
metadata <- read.table("phyloseq/weis_metadata_1210.tsv", sep = "\t", header = TRUE)  # read.table()：讀取 TSV（Tab-Separated Values）文件，分隔符為 \t & header = TRUE：表示第一行為欄位名稱

# (Method 2) 
# metadata <- read.table("D:/yuyu/master/113-1/TA/113-1 Bioinformatics/Final/phyloseq/weis_metadata_1210.tsv", sep = "\t", header = TRUE)
```
- 方法一：先指定路徑，之後只需要輸入檔案名稱即可。
- 方法二：直接輸入檔案路徑取得檔案。

## Specific Variables
```{r}
# colnames(metadata) ## What variables are in metadata
```

### 變數標籤設定
```{r}
data_table <- metadata %>%
  select(-ID) %>% ## Exclude the ID variable
  upData( ## Supplementary explanation of variable names
    labels = c(
      Calprotectin_greater_than_50 = "Calprotectin greater than 50 mcg/g",
      Gastrointestinal_symptoms = "Gastrointestinal symptoms",
      Disease_duration = "Disease duration (months)",
      Hoehn_Yahr_stage = "Hoehn-Yahr stage",
      Average_L_dopa_dose = "Average L-dopa dose (last 2 years)",
      Family_history_for_neurodegenerative_disorders = "Family history for neurodegenerative disorders"
    )
  )

```
－ 補充說明變數名稱
-  upData()（來自 Hmisc）：為變數添加標籤（label），以利後續 table1() 顯示友善的名稱。


## Select the variables that you want to compare (設定要比較的群組)
```{r}
compared_group <- "Entacapone"      # 比較的變數名稱（Entacapone）
compared_levels <- c("yes", "no")   # 該變數的兩個水平（yes/no）
compared_labels <- c("Yes", "No")   # 顯示標籤（Yes/No）
```
－ 這邊舉例的是對"Entacapone"做分組比較，分成兩組 "yes", "no" ，而我們希望它的組別以"Yes", "No"呈現


### 轉換比較變數為 factor
```{r}
data_table1 <- data_table
data_table1[compared_group] <- factor(data_table1 %>% pull(eval(parse(text = compared_group))),levels = c(compared_levels,Inf),labels = c(compared_labels,"P-value"))
data_table1 <- data_table1[!is.na(data_table1 %>% pull(eval(parse(text = compared_group)))),]
```
* 為了保留原始檔以便後續確認，我們會新建一個一模一樣的檔案，以新的檔案來執行分析。
* 在data_table1資料檔中，增加一行"compared_group" (將 compared_group 轉為 factor)
* 移除 NA 值，確保比較時不包含缺失值

### 定義 P 值計算函數 & 定義數值變數的顯示方式
```{r}
p_value <- function(x, ...) {
  #x <- x[-which(names(x) == "overall")]
  y <- unlist(x)
  g <- factor(rep(1:length(x), times=sapply(x, length)))
  if (is.numeric(y)) {
    p <- wilcox.test(y ~ g, exact = FALSE)$p.value ## If y is numeric
  } else if(table(y, g) %>% dim %>% identical(c(2, 2))) {
    p <- exact.test(table(y, g), method = "CSM")$p.value # If table(y, g) is a 2x2
  } else {
    p <- GTest(table(y, g), correct = "williams")$p.value
  }
  c("",sub("<", "&lt;", format.pval(p, digits = 3, eps = 0.001)))
}

my.render.cont <- function(x) {
  with(stats.apply.rounding(stats.default(x), digits=3), 
    c("", "Mean (SD)" = sprintf("%s (%s)", MEAN, SD))
  )
}

```
- 建立一個p_value函數，依數據的類型，使用不同方式得到p_value
  
  -  數值型數據：使用Wilcoxon 檢定（適用於非常態分佈資料）來檢測兩組數據是否來自相同的分佈
  -  2x2列聯表行數據(類別變數)：使用exact.test進行準確性檢驗 (Fisher’s exact test)
  -  其他類型的數據：使用G檢定（GTest()）進行似然比檢驗，並應用"Williams校正"來計算p-value
    
- p-value是一個統計指標，用來評估觀察到的差異是否只由隨機變異所造成的。
  -   如果 p-value 小於 0.05，則我們可以認為這個差異具有統計顯著性，意味著兩組在這個變數上的差異並非隨機產生的，有可能是真實存在的。
  -   如果 p-value 大於 0.05，則表示該變數在兩組之間的差異不具統計顯著性。


##  Table 1 (生成 Table 1)
為了了解 PD 組有無服用Entacapone是否與性別、年齡、抽菸與否等等的一般特徵相關，我們比較了 PD 組中有服用Entacapone ( N  = 11) 和未服用Entacapone ( N  = 13) 患者的幾個特徵。
```{r}
tb1 <- colnames(data_table1) %>%
  .[. != gsub("`", "", compared_group)] %>%
  paste0("`", ., "`", collapse = " + ") %>%
  paste0("~ ", ., " | ", compared_group) %>%
  as.formula %>%
  table1(data = data_table1, 
         overall = F, 
         render.continuous = my.render.cont,
         extra.col = list(`P-value` = p_value))
```

- 計算 p-value，以便評估每個變數在兩組之間是否存在顯著的差異。
- 我們得到的Table1上的p-value皆大於0.05，代表兩組 (Yes / No) 之間，在這些特徵（性別、年齡、抽菸情況等）上是沒有統計上的顯著差異。




--------------------------------------
--------------------------------------
## 補充

### Fisher’s Exact Test（費雪精確檢定）
Fisher’s Exact Test 是一種 統計檢定方法，用於檢測 兩個分類變數是否具有統計相關性，特別適用於 小樣本 或 稀疏數據（例如某些分類組合中的觀測次數過少）。
#### 1. 何時使用 Fisher’s Exact Test?
當研究的數據符合以下條件時，適合使用 Fisher’s Exact Test：
1. 分析兩個分類變數（例如：疾病 vs. 用藥情況）。
2. 樣本量小（一般來說，當列聯表中某些格子的期望值 < 5 時，卡方檢定可能不可靠）。
3. 資料呈現為 2×2 或較小的列聯表（contingency table）。
相較於 卡方檢定（Chi-square Test），Fisher’s Exact Test 不需要假設樣本量足夠大，適用於小樣本數據。

#### 2. 2×2 列聯表範例
假設研究某種藥物 Entacapone 是否影響 疾病發生率，數據如下：
|                 | 有疾病 (Disease) | 無疾病 (No Disease) |總計|
|---------------|-------------------|----------------|----------------|
| **使用藥物 (Entacapone = Yes)**  | 3 | 1 | 4 |
| **未使用藥物 (Entacapone = No)**  | 1 | 5 | 6 |
| **總計** | 4 | 6 | 10 |

- 我們希望知道 Entacapone 的使用是否與疾病發生有統計上的關聯。


#### 3. 計算原理
Fisher’s Exact Test 基於超幾何分佈（Hypergeometric Distribution），計算 觀察到這種數據的精確機率。

假設 行總計與列總計固定，計算所有可能的 2×2 表格組合，並比較目前這個表格的發生機率是否顯著低。
計算 P 值 來檢定這個組合的 極端程度。
Fisher’s Exact Test 計算的 P 值是：
P = \frac{\binom{a+b}{a} \binom{c+d}{c}}{\binom{n}{a+c}}



