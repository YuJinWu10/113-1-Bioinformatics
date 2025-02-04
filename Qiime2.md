為了研究服用Entacapone是否與腸道菌群結構的變化相關，我們進行多樣性分析，這邊我們做Alpha diversity和Beta Diversity

QIIME2 微生物組學分析，包括 Alpha 多樣性、Beta 多樣性分析、PERMANOVA 檢定

# Install the required packages
```{r}
#if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#BiocManager::install("MicrobiotaProcess")

#if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#BiocManager::install("phyloseq")
```


# Library the required packages
```{r}
library(MicrobiotaProcess)      # MicrobiotaProcess：處理 QIIME2 轉換的 phyloseq 物件
library(phyloseq)               # phyloseq：微生物組學的標準數據處理
library(tidyverse)              # tidyverse：數據處理（ggplot2、dplyr、tidyr 等）
library(treeio)                 # treeio & ape：處理微生物系統發生樹（phylogenetic tree）
library(ggh4x)
library(ape)
library(coin)                   # coin：Wilcoxon 檢定
library(vegan)                  # vegan：PERMANOVA（adonis2）檢定
library(qiime2R)                # qiime2R：載入 QIIME2 產生的 .qza 文件
```



# Prepare data
## 讀取 QIIME2 資料 & 導入 QIIME2 數據為 phyloseq 物件
```{r}
setwd("D:/yuyu/master/113-1/TA/113-1 Bioinformatics/Final") #　設定工作目錄
otu <- "phyloseq/feature-table.qza"
rep <- "phyloseq/rep-seqs.qza"
tree <- treeio::read.treeqza("phyloseq/rooted-tree.qza")
tree <- ape::multi2di(tree)
tax <- "phyloseq/taxonomy.qza"
sample <- "phyloseq/weis_metadata_1210.tsv"

## 導入 QIIME2 數據為 phyloseq 物件
ps <- MicrobiotaProcess::import_qiime2(otuqza = otu,
                                       taxaqza = tax,
                                       refseqqza = rep,
                                       mapfilename = sample,
                                       treeqza = tree)

```
- 讀取 .qza 格式的 OTU 表、代表性序列、分類表、系統發生樹
- multi2di()：將樹轉換為二叉樹（避免多重分叉）
- import_qiime2() 轉換 QIIME2 格式為 phyloseq 物件
- ps 是後續分析的主體


# alpha diversity
- alpha diversity是指在單一樣本或區域內的物種豐富度和多樣性。它衡量的是單一樣本內物種的多樣性，常見的指標有Observed, Chao1, Shannon, Simpson
  - Observed：實際觀察到的物種數量
  - Chao1：基於物種豐富度的估計，考慮到未被觀察到的稀有物種。
  - Shannon：基於物種的數量和均勻度計算多樣性，值越大表示多樣性越高。
  - Simpson：另一種多樣性指標，值越接近1，表示多樣性越高。
    
## alpha diversity plot: Observed, Chao1, Shannon, Simpson
```{r}
# 計算 Alpha 多樣性
alphaobj <- get_alphaindex(ps) %>% as.data.frame
head(as.data.frame(alphaobj))
colnames(alphaobj)[which(colnames(alphaobj) == "Observe")] <- "Observed"

## 繪製 Alpha 多樣性箱型圖
alphaobj %>%
  pivot_longer(c("Observed", "Chao1", "Shannon", "Simpson"), names_to = "metric", values_to = "Alpha Diversity Measure") %>%
  mutate(metric = factor(metric, levels = c("Observed", "Chao1", "Shannon", "Simpson"))) %>%
  ggplot() +
    geom_boxplot(aes(y = `Alpha Diversity Measure`, fill = Entacapone), outlier.shape = 21, width = 5) +
  facet_grid2(.~metric, scales = "free_y", independent = "y") +
  theme_minimal() +
  theme(axis.text=element_text(size=10), 
        axis.title=element_text(size=10)
  )+
  scale_x_discrete(labels = element_blank(), name = "") +
  scale_fill_manual("Entacapone", values=c("#E69F00", "#0072B2")) 
```
- 生成了Alpha diversity測量的數據框alphaobj，並將「Observe」變量重新命名為「Observed」。
- 將多個Alpha diversity指標（Observed, Chao1, Shannon, Simpson）轉換為長格式，並生成盒狀圖比較不同樣本群組（Entacapone）下的多樣性。
- fill='你想要比對的變數名稱'
- get_alphaindex(ps)：計算 Alpha 多樣性指標（Observed、Chao1、Shannon、Simpson）
- pivot_longer()：將指標轉為長格式
- geom_boxplot()：繪製箱型圖
- facet_grid2()：按指標分類
- scale_fill_manual()：設定 Entacapone 組別的顏色



## exact wilcoxon test for alpha diversity (Wilcoxon 檢定)
```{r}
metric <- c("Observed", "Chao1", "Shannon", "Simpson")
alpha_diversity <- data.frame()
for(i in seq_along(metric)) {
  alpha_diversity[1,i] <- coin::wilcox_test(reformulate("factor(Entacapone)", metric[[i]]), data = alphaobj) %>% pvalue()
}

colnames(alpha_diversity) <- metric
rownames(alpha_diversity) <- "p.value"
alpha_diversity
```
 - wilcox_test()：檢測不同 Entacapone 群組間的 Alpha 多樣性是否有統計差異
 - 透過wilcoxon test進行每個Alpha diversity指標的Wilcoxon秩和檢驗，並返回p-value，這部分用來比較兩個群組在每個指標上的顯著性差異。
 - Alpha diversity四個指標（Observed, Chao1, Shannon, Simpson）的P-value都大於0.05，也就是說Alpha diversity四個指標（Observed, Chao1, Shannon, Simpson）在服用Entacapone和未有服用Entacapone的患者之間沒有發現統計上顯著差異（pObserved = 0.6223937；pChao1 = 0.4689388；pShannon = 0.5429662；pSimpson = 0.3106355） 。
 - 不管是它的豐富度還是均勻度

# beta diversity
- Beta Diversity是描述不同樣本或群落之間的物種差異，表示不同群落之間的多樣性差異。常見的指標有canberra, unweighted Unifraca, weighted Unifrac
  -  Canberra 距離：衡量兩個樣本中物種豐富度的差異。
  -  Unifrac（Unweighted 和 Weighted）：基於系統發生樹（phylogenetic tree）計算的樣本之間的差異。Unweighted Unifrac距離只考慮是否存在物種，而 Weighted Unifrac距離會考慮物種的豐富度。
(生成三個不同距離度量的PCoA圖（主坐標分析圖），這些圖用來可視化群組之間的Beta多樣性差異，每個PCoA圖都顯示了樣本在兩個主坐標軸上的分佈，並且透過橢圓形來描繪群組的分佈範圍)

- Alpha Diversity 和 Beta Diversity 的比較：
  - Alpha Diversity 專注於單一樣本內的多樣性，通常用來衡量群落內的物種豐富度和均勻性。如果你想知道某個樣本內有多少不同的物種，以及這些物種分佈是否均勻，Alpha Diversity 是適合的指標。
  - Beta Diversity 則是用來比較不同樣本之間的物種組成差異。它更多的是比較兩個或多個樣本群落之間的多樣性變化，從而反映它們的相似性或差異。

## beta diversity plot: Canberra PCoA
```{r}
pcoa_canberra <- get_pcoa(obj = ps, 
                          distmethod = "canberra", 
                          method = "hellinger")
beta_canberra <- cbind(ps@sam_data, pcoa_canberra@pca$vectors)

## 繪圖
ggplot(beta_canberra, aes(x = Axis.1, y = Axis.2, color = Entacapone)) +
  geom_point() +
  stat_ellipse(type = "t") +
  scale_colour_manual("Entacapone", values=c("#E69F00", "#0072B2")) +
  labs(title = "Canberra PCoA",
       x = paste0("Axis.1 [", (pcoa_canberra@pca$values["Relative_eig"][1,1] * 100) %>% round(1), "%]"),
       y = paste0("Axis.2 [", (pcoa_canberra@pca$values["Relative_eig"][2,1] * 100) %>% round(1), "%]")) +
  theme_classic()
```
- get_pcoa()：計算主坐標分析（PCoA）
- cbind()：合併 PCoA 結果與 metadata

## beta diversity plot: Unweighted UniFrac PCoA
```{r}
pcoa_unweighted_unifrac <- get_pcoa(obj = ps,
                                    distmethod = "Unweighted-UniFrac",
                                    method = "hellinger")
beta_unweighted_unifrac <- cbind(ps@sam_data, pcoa_unweighted_unifrac@pca$vectors)
ggplot(beta_unweighted_unifrac, aes(x = Axis.1, y = Axis.2, color = Entacapone)) +
  geom_point() +
  stat_ellipse(type = "t") +
  scale_colour_manual("Entacapone", values=c("#E69F00", "#0072B2")) +
  labs(title = "Unweighted Unifrac PCoA",
       x = paste0("Axis.1 [", (pcoa_unweighted_unifrac@pca$values["Relative_eig"][1,1] * 100) %>% round(1), "%]"),
       y = paste0("Axis.2 [", (pcoa_unweighted_unifrac@pca$values["Relative_eig"][2,1] * 100) %>% round(1), "%]")) +
  theme_classic()
```
- Unweighted UniFrac：考慮分支長度但不計入相對豐度

## beta diversity plot: Weighted UniFrac PCoA
```{r}
pcoa_weighted_unifrac <- get_pcoa(obj = ps, 
                                  distmethod = "weighted-UniFrac", 
                                  method = "hellinger")
beta_weighted_unifrac <- cbind(ps@sam_data, pcoa_weighted_unifrac@pca$vectors)
ggplot(beta_weighted_unifrac, aes(x = Axis.1, y = Axis.2, color = Entacapone)) +
  geom_point() +
  stat_ellipse(type = "t") +
  scale_colour_manual("Entacapone", values=c("#E69F00", "#0072B2")) +
  labs(title = "Weighted Unifrac PCoA",
       x = paste0("Axis.1 [", (pcoa_weighted_unifrac@pca$values["Relative_eig"][1,1] * 100) %>% round(1), "%]"),
       y = paste0("Axis.2 [", (pcoa_weighted_unifrac@pca$values["Relative_eig"][2,1] * 100) %>% round(1), "%]")) +
  theme_classic()
```
- Weighted UniFrac：計入分支長度與相對豐度


## PERMANOVA test for beta diversity(Calculate p-values) [PERMANOVA 分析]
- PERMANOVA : Permutational analysis of variance
- 使用adonis2進行PERMANOVA（基於距離矩陣的方差分析）來比較不同群組之間的Beta diversity顯著性，並生成p-value。
```{r}
metric <- c("canberra", "unweighted-UniFrac", "weighted-UniFrac")
beta_diversity <- data.frame()
for(i in seq_along(metric)) {
  distme <- get_dist(ps, distmethod = metric[[i]], method="hellinger")
  sampleda <- data.frame(sample_data(ps), check.names=FALSE)
  sampleda <- sampleda[match(colnames(as.matrix(distme)),rownames(sampleda)),,drop=FALSE]
  sampleda$Entacapone <- factor(sampleda$Entacapone)
  set.seed(1234)
  beta_diversity[1,i] <- adonis2(distme ~ Entacapone, data = sampleda, permutations = 9999)$`Pr(>F)`[1]
}
colnames(beta_diversity) <- metric
rownames(beta_diversity) <- "p.value"
beta_diversity
```
-Beta Diversity 的三個指標中的兩個（canberra距離和加權Unifrac距離）顯示服用Entacapone和未有服用Entacapone的患者之間的群落結構發生顯著變化（pCanberra = 0.0371 ； pUnweighted-unifrac = 0.1859；pweighted-unifrac = 0.0048)。


- Alpha Diversity 強調的是單個樣本內的多樣性，可以用於測量每個樣本群落的物種豐富度和均勻性。
- Beta Diversity 則比較不同樣本之間的物種差異，測量樣本之間的群落變異性。
