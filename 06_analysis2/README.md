# Analyze data

The scripts below were used to analyze the data. Given that this meta-analysis requires high computational resources, we decided to release the processed data together with the code to perform the analyses. All dataset are available inside the folder `02_processed_data`.

### Load libraries

```R
library("tidyverse")
library("metafor")
library("multcomp")
library("lme4")
library("car")
library("emmeans")
library("lmerTest")
```

### Diversity data

Make sure to load correct the datasets. Examples below refer to the 16S dataset.

```R
df <- read.table("02_processed_data/16S/diversity_16S.txt", header = T, sep = "\t")

calc.div.metr <- function(dataset, metric){
  div.df <- dataset %>%
    group_by(BioProject, group) %>%
    summarize(div = mean(get(metric)),
              sd = sd(get(metric)),
              n = length(get(metric))) %>%
    filter(div > 0) %>%
    pivot_wider(names_from = group, values_from = c(div, sd, n)) %>% filter(!is.na(div_Unburned))
  
  df <- metafor::escalc(measure = "SMD", 
                        vtype = "UB", 
                        data = div.df,
                        m1i = div_Burned,
                        m2i = div_Unburned, 
                        sd1i = sd_Burned,
                        sd2i = sd_Unburned,
                        n1i = n_Burned,
                        n2i = n_Unburned)
  names(df)[names(df) == 'yi']  <- "effsize"
  names(df)[names(df) == 'vi'] <- "variance"
  
  model.1 <- rma(method="REML", yi = df$effsize, vi = df$variance)
  
  b <- data.frame(metric = metric,
                  EffSize = model.1$b,
                  se = model.1$se,
                  QM = model.1$QM,
                  pval = model.1$pval)
  
  return(b)
}

list.metrics <- c("observed", "chao1", "diversity_shannon", "dominance_dbp", "evenness_pielou", "rarity_log_modulo_skewness", "phylogenetic_div", "MNTD")
res <- sapply(list.metrics, calc.div.metr, dataset = df, simplify = FALSE, USE.NAMES = TRUE)
res.df <- do.call(rbind, res) %>% mutate(color = ifelse(pval < 0.05, "darkblue", "darkgray"))

res.df <- res.df %>% mutate(metric = dplyr::recode(metric,
                                            `phylogenetic_div` = "phylogenetic diversity",
                                            `rarity_log_modulo_skewness` = "rarity",
                                            `evenness_pielou` = "evenness",
                                            `dominance_dbp` = "dominance",
                                            `diversity_shannon` = "shannon diversity"))

res.df$metric <- factor(res.df$metric , levels=c("observed",
                                                 "shannon diversity",
                                                 "chao1",
                                                 "phylogenetic diversity",
                                                 "dominance",
                                                 "evenness",
                                                 "rarity",
                                                 "MNTD"))

plot <- ggplot(res.df, aes(metric, EffSize, color = color)) +
  theme_bw(base_size = 15) +
  geom_point() +
  geom_errorbar(aes(ymin=EffSize-se*1.96, ymax=EffSize+se*1.96), width=.2, position=position_dodge(0.05)) +
  theme(axis.line.x = element_line(size = 0.5),
        axis.ticks.x = element_line(size = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  scale_color_identity() +
  xlab("") +
  ylab("Effect size (±95% CI)") +
  ylim(-1.5, 1.5)+
  scale_x_discrete(limits=rev) +
  coord_flip()
plot
```

This code will outputs the following graph:

<img width="993" alt="image" src="https://github.com/user-attachments/assets/74dcd3b0-098b-43d5-bc70-729f702eafb9" />

In the case of ITS data, the code will output the following graph:

<img width="993" alt="image" src="https://github.com/user-attachments/assets/ee09be7f-fd8d-4800-a7b8-e60bc8c6c2cb" />

### iCAMP

Make sure to load correct the datasets. Examples below refer to the 16S dataset.

```R
data.icamp <- read.table("02_processed_data/16S/icamp_16S.txt", header = T, sep = "\t") %>% as_tibble()

df2 <- df %>% dplyr::select(c("SampleID", "group"))
colnames(data.icamp) <- c( "BioProject", "id", "heterog. sel.", "homog. sel.", "dispersal lim.", "homog. dispersal", "drift + others")
df.icamp <- merge(df2, data.icamp, by.x = "SampleID", by.y = "id") %>% as_tibble()

list.metrics <- colnames(df.icamp)[4:8]
res <- sapply(list.metrics, calc.div.metr, dataset = df.icamp, simplify = FALSE, USE.NAMES = TRUE)
res.df <- do.call(rbind, res) %>% mutate(color = ifelse(pval < 0.05, "darkblue", "darkgray"))

res.df$metric <- factor(res.df$metric , levels=c("heterog. sel.",
                                                 "homog. sel.",
                                                 "homog. dispersal",
                                                 "dispersal lim.",
                                                 "drift + others"))

plot <- ggplot(res.df, aes(metric, EffSize, color = color)) +
  theme_bw(base_size = 15) +
  geom_point() +
  geom_errorbar(aes(ymin=EffSize-se*1.96, ymax=EffSize+se*1.96), width=.2, position=position_dodge(0.05)) +
  theme(axis.line.x = element_line(size = 0.5),
        axis.ticks.x = element_line(size = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  scale_color_identity() +
  xlab("") +
  ylab("Effect size (±95% CI)") +
  scale_x_discrete(limits=rev) +
  coord_flip()
plot
```

This code will outputs the following graph:

<img width="993" alt="image" src="https://github.com/user-attachments/assets/419a4583-2109-4a32-98c5-7d8fb7cdb32c" />

In the case of ITS data, the code will output the following graph:

<img width="993" alt="image" src="https://github.com/user-attachments/assets/8c30ae91-b0c3-47d5-a121-58d2a6cee5d0" />

### EcolUtils

Make sure to load correct the datasets. Examples below refer to the 16S dataset.

```R
data.ecolutils <- read.table("02_processed_data/16S/ecolutils_16S.txt", header = T, sep = "\t") %>% separate(group, c("bp", "group")) %>% as_tibble()

eu.df <- data.ecolutils %>% group_by(bp, group, sign) %>% count() %>% ungroup() %>%
  group_by(bp, group) %>% mutate(tot = sum(n), prop = n/tot) %>% filter(sign != "NON SIGNIFICANT") %>%
  janitor::remove_empty("rows")

eu.df <- eu.df %>% mutate(sign = dplyr::recode(sign,
                                                   `GENERALIST` = "generalist",
                                                   `SPECIALIST` = "specialist"))

plot <- ggplot(eu.df, aes(x = sign, y = prop, fill = group)) +
  theme_bw(base_size = 15) +
  geom_boxplot(outlier.colour="black", notch=F, outlier.shape=NA) +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        axis.title.y = element_text(color="black"),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "none") +
  ylab("proportion of total ASVs") +
  scale_fill_manual(name = "Legend", values=c("#e41a1c", "#377eb8", "black"),
                    labels = c("Burned", "Unburned"),
                    breaks = c("Burned", "Unburned"))
plot
```

This code will outputs the following graph:

<img width="754" alt="image" src="https://github.com/user-attachments/assets/c752520a-f4d3-44e1-8cd7-b848e47a8690" />

In the case of ITS data, the code will output the following graph:

<img width="754" alt="image" src="https://github.com/user-attachments/assets/05680644-070e-4bda-872b-3a74332f167b" />

This code can be used to test for differences between groups:

```R
model <- lmer(prop ~ group * sign * (1|bp), data = eu.df)
Anova(model)
ph <- emmeans(model, "group", by = "sign")
pairs(ph)
```
