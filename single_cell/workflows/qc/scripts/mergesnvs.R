library(tidyverse)

args <- commandArgs(TRUE)
input = data.table::fread(args[1])
print(input)
output = args[2]


filtsnvs <- input %>%
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

write_delim(filtsnvs, output, delim = "\t")
