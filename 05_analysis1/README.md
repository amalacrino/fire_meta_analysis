# Create datasets for analysis

The scripts below create the datasets that can be used for further analysis. The code assumes that you have loaded the correct `phyloseq` object (named `ps`) in your R environment.

### Load libraries

```R
library("tidyverse")
library("phyloseq")
library("microbiome")
library("picante")
library("iCAMP")
```

### Data cleanup

```R
ps <- subset_taxa(ps, Class !="Chloroplast") #run this for 16S
ps <- subset_taxa(ps, Order !="Chloroplast") #run this for 16S
ps <- subset_taxa(ps, Family !="Mitochondria") #run this for 16S
ps <- filter_taxa(ps, function (x) {sum(x > 0) > 1}, prune=TRUE)
ps <- prune_samples(sample_sums(ps) >= 1000, ps)
```

### Calculate diversity metrics

```R
div <- microbiome::alpha(ps)
otus <- as.data.frame(t(otu_table(ps)))
tree <- phy_tree(ps)
div.pd <- pd(otus, tree, include.root = FALSE)
div.2 <- cbind(sample_data(ps), div)
div.2 <- cbind(div.2, div.pd)

mntd.calc <- function(x){
  print(paste0("############ BIOPROJECT ", x, " STARTING at ", Sys.time()))
  ks <- sample_data(ps)[["BioProject"]] %in% x
  ps <- prune_samples(samples = ks, ps)
  ps <- filter_taxa(ps, function (x) {sum(x > 0) > 1}, prune=TRUE)
  
  comm <- as.data.frame(t(otu_table(ps)))
  phy <- phy_tree(ps)
  phy.dist <- cophenetic(phy)
  mntd.df <- ses.mntd(comm, phy.dist, null.model = "richness", abundance.weighted = T, runs = 999)
  
  print(paste0("-------------- BIOPROJECT ", x, " COMPLETED at ", Sys.time()))
  return(mntd.df)
}

mntd.res <- sapply(list.gm, mntd.calc, simplify = FALSE, USE.NAMES = TRUE)
mntd.res.df <- do.call(rbind, mntd.res)
div.2 <- cbind(div.2, mntd.res.df)
write.table(div.2, "diversity.txt", quote = F, sep = "\t")
```

### iCAMP

```R
icamp.calc <- function(x){
  print(paste0("############ BIOPROJECT ", x, " STARTING at ", Sys.time()))
  ks <- sample_data(ps)[["BioProject"]] %in% x
  ps <- prune_samples(samples = ks, ps)
  ps <- filter_taxa(ps, function (x) {sum(x > 0) > 1}, prune=TRUE)
  preabs <- as.data.frame(t(otu_table(ps)))
  phy <- phy_tree(ps)
  sss <- icamp.big(preabs, phy, nworker = 64)
  df <- sss$detail$processes$CbMPDiCBraya %>%
    select(-sample2) %>%
    group_by(sample1) %>%
    summarise(across(everything(), mean))
  df$BioProject <- paste0(x)
  write.table(df, paste0("icamp_out/", x, ".txt"), quote = F, sep = "\t")
  unlink("path.rda")
  unlink("pd.bin")
  unlink("pd.desc")
  unlink("pd.taxon.name.csv")
  unlink("iCAMP.process.CbMPDiCBraya.csv")
  unlink("iCAMP.iCAMP.Confidence.detail.rda")
  print(paste0("-------------- BIOPROJECT ", x, " COMPLETED at ", Sys.time()))
  return(df)
  }

icamp.res <- sapply(list.gm, icamp.calc, simplify = FALSE, USE.NAMES = TRUE)
icamp.res.df <- do.call(rbind, icamp.res)
write.table(icamp.res.df, "icamp.txt", quote = F, sep = "\t")
```

### EcolUtils

```R
sample_data(ps)$group2 <- paste0(ks.tmp$BioProject,"-",ks.tmp$Burned_Unburned)
list.gm <- unique(sample_data(ps)$group2)

sncm.calc <- sapply(list.gm,  
                    function(x){
                      print(paste0("############ BIOPROJECT ", x, " STARTING at ", Sys.time()))
                      ks <- sample_data(ps)[["group2"]] %in% x
                      ps <- prune_samples(samples = ks, ps)
                      ps <- filter_taxa(ps, function (x) {sum(x > 0) > 1}, prune=TRUE)
                      otut <- t(as.data.frame(otu_table(ps)@.Data))
                      spp.out <- spec.gen(otut, niche.width.method = "levins", perm.method = "quasiswap",
                                          n = 10, probs = c(0.025, 0.975))
                      spp.out <- tibble::rownames_to_column(spp.out, "ASV")
                      spp.out$group <- paste0(x)
                      print(paste0("-------------- BIOPROJECT ", x, " COMPLETED at ", Sys.time()))
                      Sys.time()
                      return(spp.out)},
                    simplify = FALSE,USE.NAMES = F)

sncm.results <- do.call(rbind, sncm.calc)
write.table(sncm.results, "ecolutils.txt", quote = F, sep = "\t")
```
