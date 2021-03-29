Data preparation
================

This document describes steps to process the gigantic EIGENSTRAT data
combining the Reich lab ancient DNA genotypes, SGDP genotypes and
Vindija and Altai Neanderthal genomes into a smaller, more manageable
EIGENSTRAT data for this course.

If you participated in this course, you don’t need to concern yourselves
with this. I used this to generated data for the *admixr* demonstration
and also to have something to share with you, in case you would like to
test *admixr* on your own.

## Download full EIGENSTRAT data used in [Petr *et al.* (2019)](https://www.pnas.org/content/116/5/1639)

This data is a combination of about 500 samples (ancient and
present-day) used in a study by [Fu *et al.*
(2016)](https://www.nature.com/articles/nature17993), as well as Vindija
and Altai Neanderthal genomes I analyzed in my paper on selection
against Neanderthal introgression.

``` bash
scp bionc04:/mnt/expressions/mp/nea-over-time/data/eigenstrat/bigyri_ho/all.{snp,ind,geno} .
```

## Prepare a table of sample names, ages, and population assignment

``` r
library(tidyverse)

sgdp <-
  read_tsv(url("http://simonsfoundation.s3.amazonaws.com/share/SCDA/datasets/10_24_2014_SGDP_metainformation_update.txt")) %>%
  select(Panel, name=SGDP_ID, Region, Country, Latitude, Longitude) %>%
  filter(complete.cases(.)) %>%
  filter(Panel == "C") %>%
  mutate(name=str_replace(name, "-", "_")) %>%
  select(-Panel) %>% 
  select(-Country, pop=Region) %>%
  mutate(age=0) %>%
  group_by(name, age, pop) %>%
  ungroup %>%
  select(-Latitude, -Longitude)

emhs <-
  read_delim(url("https://raw.githubusercontent.com/bodkan/nea-over-time/master/data/emh_ages.txt"), delim=" ", col_names=c("name", "age")) %>%
  mutate(pop="EMH") %>%
  filter(name != "Oase1") %>%
  mutate(name = case_when(
    name == "UstIshim" ~ "new_UstIshim",
    name == "Loschbour" ~ "new_Loschbour",
    TRUE ~ name
  ))

emh_keep <-
  readRDS(url("https://github.com/bodkan/nea-over-time/raw/master/data/rds/nea_estimates.rds")) %>%
  filter(stat == "direct_f4", snp_count >= 200000, pop == "EMH", sites == "all", C == "Yoruba") %>%
   .$X

modern_keep <- c(
  "S_Papuan_1",
  "S_Adygei_1", "S_Orcadian_1", "S_Spanish_1", "S_Greek_2", "S_French_1", "S_Russian_1", "S_Icelandic_2", "S_Estonian_1",   "S_Czech_2", "S_Tuscan_1",
  "S_Han_1", "S_Japanese_1", "S_Korean_2", "S_Dai_3", "S_Burmese_2", "S_Tujia_1", "S_Uygur_1", "S_Hezhen_2",
  "S_Dinka_1", "S_Mbuti_1", "S_Yoruba_1"
)

samples <- bind_rows(emhs, sgdp) %>%
  filter(name %in% c(emh_keep, modern_keep))

write_tsv(samples, "samples.tsv")
```

## Subset the gigantic EIGENSTRAT data to a managable size

``` r
library(admixr)

all_snps <- eigenstrat("all")

geno <- read_geno(all_snps)
ind <- read_ind(all_snps)

ind_emh <- filter(ind, label %in% emh_keep)
ind_sgdp <- filter(ind, id %in% str_replace(modern_keep, "_(\\d+)$", "-\\1"))
ind_other <- filter(ind, label %in% c("new_Vindija", "new_Altai", "Chimp"))

new_ind <- bind_rows(ind_emh, ind_sgdp, ind_other)
new_geno <- geno[, new_ind$id] %>%
  rename(Vindija = new_Vindija,
         Altai = new_Altai)

new_ind <- mutate(
  new_ind,
  id = str_replace(id, "^new_", ""),
  label = str_replace(label, "^new_", "")
) %>%
  mutate(id = str_replace("^S_", "", id)) %>%
  mutate(id = str_replace("-\\d+$", "", id))

write_geno(new_geno, "subset.geno")
write_ind(new_ind, "subset.ind")
file.copy(all_snps$snp, "subset.snp")
```

## Clean up

``` bash
rm all.{geno,ind,snp}
```
