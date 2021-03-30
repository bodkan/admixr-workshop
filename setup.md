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
mkdir data/
scp bionc04:/mnt/expressions/mp/nea-over-time/data/eigenstrat/bigyri_ho/all.{snp,ind,geno} data/
```

## Prepare a table of sample names, ages, and population assignment

``` r
library(tidyverse)
```

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.0 ──

    ## ✓ ggplot2 3.3.3     ✓ purrr   0.3.4
    ## ✓ tibble  3.1.0     ✓ dplyr   1.0.2
    ## ✓ tidyr   1.1.2     ✓ stringr 1.4.0
    ## ✓ readr   1.4.0     ✓ forcats 0.5.0

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

``` r
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
```

    ## Warning: Missing column names filled in: 'X12' [12], 'X13' [13], 'X14' [14],
    ## 'X15' [15], 'X16' [16], 'X17' [17], 'X18' [18], 'X19' [19], 'X20' [20],
    ## 'X21' [21], 'X22' [22], 'X23' [23], 'X24' [24]

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   .default = col_logical(),
    ##   Panel = col_character(),
    ##   SGDP_ID = col_character(),
    ##   Population_ID = col_character(),
    ##   Region = col_character(),
    ##   Country = col_character(),
    ##   Contributor = col_character(),
    ##   Gender = col_character(),
    ##   Latitude = col_double(),
    ##   Longitude = col_double(),
    ##   Coverage = col_double(),
    ##   HetRateAuto = col_double()
    ## )
    ## ℹ Use `spec()` for the full column specifications.

``` r
emhs <-
  read_delim(url("https://raw.githubusercontent.com/bodkan/nea-over-time/master/data/emh_ages.txt"), delim=" ", col_names=c("name", "age")) %>%
  mutate(pop="EMH") %>%
  filter(name != "Oase1")
```

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   name = col_character(),
    ##   age = col_double()
    ## )

``` r
emh_keep <-
  readRDS(url("https://github.com/bodkan/nea-over-time/raw/master/data/rds/nea_estimates.rds")) %>%
  filter(stat == "direct_f4", snp_count >= 200000, pop == "EMH", sites == "all", C == "Yoruba") %>%
   .$X

modern_keep <- c(
  "S_Papuan_1",
  "S_Adygei_1", "S_Orcadian_1", "S_Spanish_1", "S_Greek_2", "S_French_1", "S_Russian_1", "S_Icelandic_2", "S_Estonian_1",   "S_Czech_2", "S_Tuscan_1",
  "S_Han_1", "S_Japanese_1", "S_Korean_2", "S_Dai_3", "S_Burmese_2", "S_Tujia_1", "S_Uygur_1", "S_Hezhen_2",
  "S_Dinka_1", "S_Mbuti_1", "S_Yoruba_1", "S_Saharawi_1"
)

# save a combined table of EMH and present-day samples with clean names
samples <- bind_rows(emhs, sgdp) %>%
  filter(name %in% c(emh_keep, modern_keep)) %>%
  mutate(name = str_replace_all(name, "^S_", "")) %>%
  mutate(name = str_replace_all(name, "_\\d+$", ""))

write_tsv(samples, "data/samples.tsv")
```

## Subset the gigantic EIGENSTRAT data to a managable size

``` r
library(admixr)

all_snps <- eigenstrat("data/all")

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
  mutate(id = str_replace(id, "^S_", "")) %>%
  mutate(id = str_replace(id, "-\\d+$", ""))

# call random alleles in diploid individuals
for (s in colnames(new_geno)) {
  # get the positions of het sites in a given individual
  het_pos <- !is.na(new_geno[[s]]) & new_geno[[s]] == 1

  if (!any(het_pos, na.rm = T)) next

  # randomly sample hom REF/ALT alleles at het positions
  new_geno[[s]] <- replace(new_geno[[s]],
                           het_pos,
                           2 * rbinom(sum(het_pos, na.rm = T), 1, 0.5)) %>% as.integer
}

write_geno(new_geno, "data/snps.geno")
write_ind(new_ind, "data/snps.ind")
file.copy(all_snps$snp, "data/snps.snp")
```

    ## [1] TRUE

## Clean up

``` bash
rm data/all.{geno,ind,snp}
```
