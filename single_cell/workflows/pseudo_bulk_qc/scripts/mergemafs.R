#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)
args <- commandArgs(TRUE)

input = args[1]
output = args[2]
maf = data.table::fread(input)

filtmaf <- filter(maf, str_detect(Consequence, "frameshift|stop") | IMPACT == "HIGH") %>%
    group_by_at(vars(-contains("depth"), -contains("count"))) %>%
    summarise(t_depth = sum(t_depth),
           t_ref_count = sum(t_ref_count),
           t_alt_count = sum(t_alt_count),
           n_depth = sum(n_depth),
           n_ref_count = sum(n_ref_count),
           n_alt_count = sum(n_alt_count),
           nlibrary = n()
         ) %>%
    ungroup() %>%
    mutate(tVAF = t_alt_count / t_depth, nVAF = n_alt_count / n_depth) %>%
    dplyr::select(id, Hugo_Symbol, Chromosome, Start_Position,
      Reference_Allele, Variant_Type, Tumor_Seq_Allele1,
      Tumor_Seq_Allele2, Consequence, IMPACT, tVAF, nVAF, nlibrary) %>%
    dplyr::arrange(id, Chromosome, Start_Position)

write_delim(filtmaf, output, delim = "\t")
