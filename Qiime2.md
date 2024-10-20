# Install the required packages
```{r}
#if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#BiocManager::install("MicrobiotaProcess")

#if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#BiocManager::install("phyloseq")
```


# Library the required packages
```{r}
library(MicrobiotaProcess)
library(phyloseq)
library(tidyverse)
library(treeio)
library(ggh4x)
library(ape)
library(coin)
library(vegan)
```



# Prepare data
```{r}
setwd("D:/yuyu/master/113-1/TA/113-1 Bioinformatics/Final")
otu <- "phyloseq/feature-table.qza"
rep <- "phyloseq/rep-seqs.qza"
tree <- treeio::read.treeqza("phyloseq/rooted-tree.qza")
tree <- ape::multi2di(tree)
tax <- "phyloseq/taxonomy.qza"
sample <- "phyloseq/weis_metadata_1210.tsv"
ps <- MicrobiotaProcess::import_qiime2(otuqza = otu,
                                       taxaqza = tax,
                                       refseqqza = rep,
                                       mapfilename = sample,
                                       treeqza = tree)

```



# alpha diversity
- alpha diversity是指在單一樣本或區域內的物種豐富度和多樣性。它衡量的是單一樣本內物種的多樣性，常見的指標有Observed, Chao1, Shannon, Simpson
  - Observed：實際觀察到的物種數量
  - Chao1：基於物種豐富度的估計，考慮到未被觀察到的稀有物種。
  - Shannon：基於物種的數量和均勻度計算多樣性，值越大表示多樣性越高。
  - Simpson：另一種多樣性指標，值越接近1，表示多樣性越高。
    
## alpha diversity plot: Observed, Chao1, Shannon, Simpson
```{r}
alphaobj <- get_alphaindex(ps) %>% as.data.frame
head(as.data.frame(alphaobj))
colnames(alphaobj)[which(colnames(alphaobj) == "Observe")] <- "Observed"
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


## exact wilcoxon test for alpha diversity
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
 - 透過wilcox test進行每個Alpha diversity指標的Wilcoxon秩和檢驗，並返回p-value，這部分用來比較兩個群組在每個指標上的顯著性差異。

 - 結論：P-value都大於0.05，也就是說對我們總體來說有服用Entacapone跟沒有服用Entacapone的帕金氏斯症患者中沒有太大的差異性，不管是它的豐富度還是

# beta diversity
- Beta Diversity是描述不同樣本或群落之間的物種差異，表示不同群落之間的多樣性差異。常見的指標有canberra, unweighted Unifraca, weighted Unifrac
  -  Canberra 距離：基於兩個樣本中物種豐富度的距離衡量。
  -  Unifrac（Unweighted 和 Weighted）：基於系統發生樹（phylogenetic tree）計算的樣本之間的差異。Unweighted Unifrac距離只考慮是否存在物種，而 Weighted Unifrac距離會考慮物種的豐富度。
(生成三個不同距離度量的PCoA圖（主坐標分析圖），這些圖用來可視化群組之間的Beta多樣性差異，每個PCoA圖都顯示了樣本在兩個主坐標軸上的分佈，並且透過橢圓形來描繪群組的分佈範圍)

- Alpha Diversity 和 Beta Diversity 的比較：
  - Alpha Diversity 專注於單一樣本內的多樣性，通常用來衡量群落內的物種豐富度和均勻性。如果你想知道某個樣本內有多少不同的物種，以及這些物種分佈是否均勻，Alpha Diversity 是適合的指標。
  - Beta Diversity 則是用來比較不同樣本之間的物種組成差異。它更多的是比較兩個或多個樣本群落之間的多樣性變化，從而反映它們的相似性或差異。

## beta diversity plot: canberra
```{r}
pcoa_canberra <- get_pcoa(obj = ps, 
                          distmethod = "canberra", 
                          method = "hellinger")
beta_canberra <- cbind(ps@sam_data, pcoa_canberra@pca$vectors)
ggplot(beta_canberra, aes(x = Axis.1, y = Axis.2, color = Entacapone)) +
  geom_point() +
  stat_ellipse(type = "t") +
  scale_colour_manual("Entacapone", values=c("#E69F00", "#0072B2")) +
  labs(title = "Canberra PCoA",
       x = paste0("Axis.1 [", (pcoa_canberra@pca$values["Relative_eig"][1,1] * 100) %>% round(1), "%]"),
       y = paste0("Axis.2 [", (pcoa_canberra@pca$values["Relative_eig"][2,1] * 100) %>% round(1), "%]")) +
  theme_classic()
```

## beta diversity plot: unweighted Unifraca
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


## beta diversity plot: weighted Unifrac
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


## PERMANOVA test for beta diversity
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

- Alpha Diversity 強調的是單個樣本內的多樣性，可以用於測量每個樣本群落的物種豐富度和均勻性。
- Beta Diversity 則比較不同樣本之間的物種差異，測量樣本之間的群落變異性。
