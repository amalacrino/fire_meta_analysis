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

```

This code will outputs the following graph:


In the case of ITS data, the code will output the following graph:

### EcolUtils

Make sure to load correct the datasets. Examples below refer to the 16S dataset.

```R

```

This code will outputs the following graph:


In the case of ITS data, the code will output the following graph:
