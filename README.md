# 113-1-Bioinformatics

## Install the required packages


```{r}
# install.packages(c("dplyr", "table1", "Hmisc", "Exact", "DescTools"))
```

## Library the required packages
```{r}
library(dplyr)
library(table1)
library(Hmisc)
library(Exact)
library(DescTools)
```

## Import dataset 

```{r}
setwd("D:/yuyu/master/113-1/TA/113-1 Bioinformatics/Final")
metadata <- read.table("phyloseq/weis_metadata_1210.tsv", sep = "\t", header = TRUE)

# (Method 2) 
# metadata <- read.table("D:/yuyu/master/113-1/TA/113-1 Bioinformatics/Final/phyloseq/weis_metadata_1210.tsv", sep = "\t", header = TRUE)
```


## Specific Variables
```{r}
# colnames(metadata) ## What variables are in metadata
```

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
－ 這邊我們補充說明變數名稱


## Select the variables that you want to compare
```{r}
compared_group <- "Entacapone"
compared_levels <- c("yes", "no")
compared_labels <- c("Yes", "No")
```
－　這邊舉例的是對"Entacapone"做分組比較，分成兩組 "yes", "no" ，而我們希望它的組別以"Yes", "No"呈現



```{r}
data_table1 <- data_table
data_table1[compared_group] <- factor(data_table1 %>% pull(eval(parse(text = compared_group))),levels = c(compared_levels,Inf),labels = c(compared_labels,"P-value"))
data_table1 <- data_table1[!is.na(data_table1 %>% pull(eval(parse(text = compared_group)))),]
```
－　為了保留原始檔以便後續確認，我們會新建一個一模一樣的檔案，以新的檔案來執行分析。
－　在data_table1資料檔中，增加一行"compared_group"
－　刪除NA值

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
  
-  數值型數據：使用wilcox.test來檢測兩組數據是否來自相同的分佈
-  2x2列聯表行數據：使用exact.test進行準確性檢驗
-  其他類型的數據：使用GTest進行似然比檢驗，並應用"Williams校正"來計算p值

##  table 1

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
－　生成Table1，可以進行分析
