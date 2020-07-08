library(tidyverse)

print("Read in files:")
print("Files:")
print(snakemake@input)

maf <- data.frame()
for (f in snakemake@input){
  print(paste0("Sample: ",str_split(basename(f), "[.]")[[1]][1]))
  maftemp <- read_delim(f, delim = "\t", guess_max = 10^7, col_types = list(PHENO = "c", SOMATIC = "c")) %>%
    mutate(id = str_split(basename(f), "[.]")[[1]][1]) %>%
    filter(t_alt_count > 2)
  maf <- bind_rows(maf, maftemp)
}

write_delim(maf, snakemake@output[[1]], delim = "\t")

print("Filter for high impact variants")

filtmaf <- filter(maf, str_detect(Consequence, "frameshift|stop") | IMPACT == "HIGH") %>%
    mutate(tVAF = t_alt_count / t_depth, nVAF = n_alt_count / n_depth) %>%
    dplyr::select(id, Hugo_Symbol, Chromosome, Start_Position,
      Reference_Allele, Variant_Type, Tumor_Seq_Allele1,
      Tumor_Seq_Allele2, Consequence, IMPACT, tVAF, nVAF) %>%
    dplyr::arrange(id, Chromosome, Start_Position)

write_delim(filtmaf, snakemake@output[[2]], delim = "\t")
