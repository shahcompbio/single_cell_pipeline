library(tidyverse)

print("Read in files:")
print("Files:")
print(snakemake@input)

snvs <- data.frame()
for (f in snakemake@input){
  print(paste0("Sample: ",str_split(basename(f), "[.]")[[1]][1]))
  snvstemp <- read_csv(f, guess_max = 10^7, col_types = list(chrom = "c")) %>%
    mutate(id = str_split(basename(f), "[.]")[[1]][1]) %>%
    filter(alt_counts > 1)
  snvs <- bind_rows(snvs, snvstemp)
}

filtsnvs <- snvs %>%
    group_by_at(vars(-contains("counts"), -num_cells)) %>%
    summarise(alt_counts = sum(alt_counts),
           ref_counts = sum(ref_counts),
           total_counts = sum(total_counts),
           num_cells = sum(num_cells),
           nlibrary = n()
         ) %>%
    ungroup() %>%
    mutate(tVAF = alt_counts / total_counts) %>%
    dplyr::select(chrom,coord,ref,alt,gene_name,effect,effect_impact,is_cosmic,
      amino_acid_change,num_cells,alt_counts,ref_counts,total_counts,
      id, tVAF, nlibrary) %>%
    dplyr::arrange(id, chrom, coord)

write_delim(filtsnvs, snakemake@output[[1]], delim = "\t")
