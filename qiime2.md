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


# beta diversity
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
