Wilcoxon 檢定與火山圖分析，針對 ASV（Amplicon Sequence Variant）和 Genus（屬級）層級的微生物豐度進行統計檢定與視覺化。主要包括：
1. 資料準備與清理
2. Wilcoxon 檢定與 T 檢定
3. Fold Change 計算
4. 火山圖（Volcano Plot）
5. 屬級（Genus Level）的 Wilcoxon 檢定


差異豐富度分析判斷兩個組別的豐富度是否有差異，哪些菌有顯著增加?哪些菌有顯著減少?
我們這邊用兩組不同檢定方式的Volcano plot
為了確定與服用Entacapone相關的關鍵腸道菌，對服用Entacapone和未服用Entacapone患者之間的微生物相對豐富度差異進行了比較

# Library the required packages
```{r} 
library(phyloseq)               # phyloseq：處理微生物組學數據
library(tidyverse)              # tidyverse：資料處理與視覺化
library(magrittr)               # magrittr：提供 %<>% 管道運算符
library(openxlsx)               # openxlsx：處理 Excel 存檔
library(MicrobiotaProcess)      # MicrobiotaProcess：導入 QIIME2 資料
library(treeio)                 # treeio、ape：處理系統發生樹
library(ape)            
library(stringr)                # stringr：處理字串（屬級名稱清理）
library(coin)                   # coin：Wilcoxon 檢定
library(qiime2R)                # qiime2R：讀取 .qza 格式的 QIIME2 結果
```


# Wilcoxon in asv
## prepare data

```{r}
setwd("D:/yuyu/master/113-1/TA/113-1 Bioinformatics/Final") # 設定工作目錄
otu <- "phyloseq/rel-feature-table.qza"
tax <- "phyloseq/taxonomy.qza"
metadata <- read.table("phyloseq/weis_metadata_1210.tsv", sep = "\t", header = TRUE)
asv <- otu %>% MicrobiotaProcess::read_qza(.) %>% .$otutab %>% data.frame
tax %<>% MicrobiotaProcess::read_qza(.) %>% .$taxatab  %>% data.frame

## 屬（Genus）與門（Phylum）標籤整理
Genus <-  str_split(tax$Genus, "g__", simplify = TRUE)[,2]
Genus[Genus  == ""] <- "Unidentified"
Phylum <- str_split(tax$Phylum, "p__", simplify = TRUE)[,2]
Phylum[Phylum  == ""] <- "Unidentified"
Genus %<>% str_replace_all(c("-" = "_", "\\[" = "", "\\]" = ""))
Phylum %<>% str_replace_all(c("-" = "_", "\\[" = "", "\\]" = ""))
```
- 設定工作目錄並讀取 OTU 表（相對豐度）、分類表（taxonomy）、元數據（metadata）
- read_qza()：從 .qza 格式讀取資料
- otutab：ASV 表
- taxatab：分類表
- 去除 QIIME2 給定的分類前綴 (g__, p__)
- 無法分類的標記為 Unidentified
- 清理特殊字元（- 轉為 _，去除 [、]）


##　 過濾低豐富度 ASV
```{r}
rownames(asv) <- paste0("asv", 1:nrow(asv))
Group <- ifelse(metadata$Entacapone == "yes", 0L, 1L)
few <- apply(asv, 1, function(x) mean(x > 0) < 0.1)
Genus <- Genus[!few]
Phylum <- Phylum[!few]
asv_new <- asv[!few,]
asv_new <- asv_new[,colnames(asv_new) %in% metadata$ID]
```
－ few：篩選 至少在 10% 樣本中出現的 ASV
－ 去除低豐度 ASV
- 只保留 metadata 中有對應 ID 的樣本


## 計算 Fold Change
```{r}
new_data <- cbind(Group, t(asv_new)) %>% as.data.frame
new_data <- new_data[!is.na(new_data$Group),]
new_data$Group <- as.factor(new_data$Group)
k <- aggregate(new_data[2:ncol(new_data)], list(new_data$Group), mean)
tmp <- apply(k, 2, function(x) all(x==0))
k <- k[, !tmp]
new_data <- new_data[, !tmp]
Genus <- Genus[!tmp[-1]]
Phylum <- Phylum[!tmp[-1]]
k <- k[,-1]
```
- aggregate()：計算每組的 ASV 平均豐度
- 移除所有為 0 的 ASV
  
```{r}
result <- data.frame()
for (i in 1:ncol(k)){
  FoldChange <- (k[2,i] + .Machine$double.eps)/(k[1,i] + .Machine$double.eps)
  result[i,1] <- FoldChange
}
```
- Fold Change 計算：(Treatment + ε) / (Control + ε)
- 加入機器最小值 .Machine$double.eps 避免除以 0






# wilcoxon test or T test

## wilcoxon test (二選一)
####　 Wilcoxon 檢定，用於比較兩組的 ASV 豐度分佈是否有統計顯著性
```{r}
# for (i in 1:ncol(k)){
#   formula <- reformulate("Group", response =  colnames(new_data)[i+1])
#   p <- coin::wilcox_test(formula,data=new_data) %>% pvalue
#   result[i,2] <- p
# }
```

或

## T test (二選一)
```{r}
for (i in 1:ncol(k)){
  formula <- reformulate("Group", response =  colnames(new_data)[i+1])
  equal_var <- bartlett.test(formula, data=new_data)$p.value
  if(equal_var > 0.05){
      p <- t.test(formula, data=new_data, var.equal=TRUE)$p.value

   }else{
      p <- t.test(formula, data=new_data, var.equal=FALSE)$p.value
   }
  result[i,2] <- p
}
```

- Bartlett’s test：檢查變異數是否相等
- 若變異數相等 (p > 0.05)：使用 pooled variance T 檢定
- 若變異數不相等：使用 Welch’s T 檢定

```{r}
for (i in 1:ncol(k)){
  criteria <- ifelse(k[1,i] == 0 | k[2,i] == 0, "1", "0")
  result[i,3] <- criteria
}

result$Phylum <- Phylum
result$Genus <- Genus
colnames(result) <- c("FoldChange","pvalue","Criteria","Phylum","Genus")
result <- result[order(result$pvalue),]

```







# volcano plot
## colored by phylum
```{r}
result$sig_p <- result$Phylum %>%
  ifelse(is.na(.), "Unidentified", .) %>%
  ifelse(result$pvalue < 0.05, ., NA)
result$sig_p %<>% factor(levels = {
  result$sig_p %>% .[. != "Unidentified"] %>% factor %>% levels %>% c(., "Unidentified")
})

result %>%
  ggplot() +
  geom_point(aes(x = log2(FoldChange), y = -log10(pvalue)), color = "gray", size = 1) +
  geom_point(data = result[result$pvalue < 0.05,], aes(x = log2(FoldChange), y = -log10(pvalue), color = sig_p), size = 1) +
  labs(x = "log2 fold change", y = expression(-log10*paste(" ","p-value")), color="Phylum") +
  scale_y_continuous(breaks = seq(from = 0, to = 3, by = 0.5),limits = c(0,3)) +
  scale_x_continuous(breaks = seq(from = -45, to = 45, by = 15),limits = c(-50,50)) +
  geom_hline(yintercept = -log10(0.05),color='red',linetype = "dashed")+
  geom_vline(xintercept = max(log2(result[result$Criteria==0,1])),color='red')+
  geom_vline(xintercept = min(log2(result[result$Criteria==0,1])),color='red')+
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "black",fill=NA),
        legend.position = "right",
        axis.title = element_text(size = rel(2)),
        axis.text = element_text(size = rel(1.5)),
        legend.text = element_text(size=12),
        legend.title  = element_text(size=13),
        legend.margin = margin(17, 2, -25, 2),
        plot.margin = margin(0.6, 0.6 , 0.6, 0.3, "cm"),
        axis.title.x = element_text(vjust = -0.8),
        axis.title.y = element_text(vjust = 1.3))
```
- 火山圖設定：
  - x 軸：log2(Fold Change)
  - y 軸：-log10(p-value)
  - 紅線：顯著性閾值 (p = 0.05)
  - 顏色區分不同分類單元


## colored by genus
```{r}
result$sig_g <- result$Genus %>%
  ifelse(is.na(.), "Unidentified", .) %>%
  ifelse(result$pvalue < 0.05, ., NA)
result$sig_g %<>% factor(levels = {
  result$sig_g %>% .[. != "Unidentified"] %>% factor %>% levels %>% c(., "Unidentified")
})
result %>%
  ggplot() +
  geom_point(aes(x = log2(FoldChange), y = -log10(pvalue)), color = "gray", size = 1) +
  geom_point(data = result[result$pvalue < 0.05,], aes(x = log2(FoldChange), y = -log10(pvalue), color = sig_g), size = 1) +
  labs(x = "log2 fold change", y = expression(-log10*paste(" ","p-value")), color="Genus") +
  scale_y_continuous(breaks = seq(from = 0, to = 3, by = 0.5),limits = c(0,3)) +
  scale_x_continuous(breaks = seq(from = -45, to = 45, by = 15),limits = c(-50,50)) +
  geom_hline(yintercept = -log10(0.05),color='red',linetype = "dashed")+
  geom_vline(xintercept = max(log2(result[result$Criteria==0,1])),color='red')+
  geom_vline(xintercept = min(log2(result[result$Criteria==0,1])),color='red')+
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "black",fill=NA),
        legend.position = "right",
        axis.title = element_text(size = rel(2)),
        axis.text = element_text(size = rel(1.5)),
        legend.text = element_text(size=12),
        legend.title  = element_text(size=13),
        legend.margin = margin(17, 2, -25, 2),
        plot.margin = margin(0.6, 0.6 , 0.6, 0.3, "cm"),
        axis.title.x = element_text(vjust = -0.8),
        axis.title.y = element_text(vjust = 1.3))
```





# Wilcoxon in Genus level(屬級（Genus Level）Wilcoxon 檢定)
類似於 ASV 層級的分析，這部分將 以屬級（Genus）為單位 進行 Wilcoxon 檢定，並輸出 Excel 結果。
## prepare data
```{r}
setwd("D:/yuyu/master/113-1/TA/113-1 Bioinformatics/Final")
otug <- "phyloseq/rel-genus-table.qza"
metadata <- read.table("phyloseq/weis_metadata_1210.tsv", sep = "\t", header = TRUE)

asv <- otug %>% MicrobiotaProcess::read_qza(.) %>% .$otutab %>% data.frame
Genus <- str_split(rownames(asv), ";g__", simplify = TRUE)[,2] %>% paste0("g__", .)
Genus[Genus  == "g__"] <- paste0("Unidentified_", 1:sum(Genus == "g__"))
for(name in Genus[duplicated(Genus)] %>% unique) Genus[Genus == name] <- paste0(name, 1:sum(Genus == name))
Genus %<>% str_replace_all(c("-" = "_", "\\[" = "", "\\]" = ""))
rownames(asv) <- Genus
Group <- ifelse(metadata$Entacapone == "yes", 0L, 1L)
few <- apply(asv, 1, function(x) mean(x > 0) < 0.1)
Genus <- Genus[!few]
asv_new <- asv[!few,]
asv_new <- asv_new[,colnames(asv_new) %in% metadata$ID]


new_data <- cbind(Group, t(asv_new)) %>% as.data.frame
new_data <- new_data[!is.na(new_data$Group),]
new_data$Group <- as.factor(new_data$Group)
k <- aggregate(new_data[2:ncol(new_data)], list(new_data$Group), mean)
k <- k[, !(apply(k, 2, function(x) all(x==0)))]
k <- k[,-1]


result <- data.frame()
for (i in 1:ncol(k)){
  FoldChange <- (k[2,i] + .Machine$double.eps)/(k[1,i] + .Machine$double.eps)
  result[i,1] <- FoldChange
}
```



## wilcoxon test
```{r}
for (i in 1:ncol(k)){
  formula <- as.formula(paste(colnames(new_data)[i+1],"~","Group"))
  p <- coin::wilcox_test(formula,data=new_data) %>% pvalue
  result[i,2]<-p
}

for (i in 1:ncol(k)){
  criteria <- ifelse(k[1,i] == 0 | k[2,i] == 0, "1", "0")
  result[i,3] <- criteria
}

result$Genus <- Genus
colnames(result) <- c("FoldChange","pvalue","Criteria","Genus")
result <- result[order(result$pvalue),]
```

```{r}
# save
# wb <- createWorkbook()
# addWorksheet(wb, "All")
# writeData(wb, "All", result %>% select(FoldChange, pvalue, Genus))
# saveWorkbook(wb, file = "Wilcoxon_Genus.xlsx", overwrite = TRUE)
```


