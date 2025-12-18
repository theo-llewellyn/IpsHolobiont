#0.1 threshold
contamdf.prev <- isContaminant(ps, method="prevalence", neg="is_NC", normalize = TRUE)
table(contamdf.prev$contaminant)
hist(contamdf.prev$p, breaks = 50)
#add taxonomy
contamdf.prev %>% rownames_to_column(var = 'OTU_ID') %>%
  left_join(read_tsv('RESULTS/OTUs/OTUstaxonomy.txt')[,c('OTU_ID','species')]) -> contam_taxa

#0.2 threshold
contamdf.prev02 <- isContaminant(ps, method="prevalence", neg="is_NC", threshold=0.2)
table(contamdf.prev02$contaminant)
#check which taxa
contamdf.prev02 %>% rownames_to_column(var = 'OTU_ID') %>%
  left_join(read_tsv('RESULTS/OTUs/OTUstaxonomy.txt')[,c('OTU_ID','species')]) -> contam_taxa0.2

# OTU table
OTU_table1 <- OTU_table[-1,]
OTU_table1[-which(contamdf.prev02$contaminant),] %>%
  fwrite(., paste0(otu_output_path, "otu_table_minusNC.txt"), sep = "\t")

# Write the filtered fasta file
filtered_fasta_minusNC <- filtered_fasta[c(OTU_table1[-which(contamdf.prev02$contaminant),]$OTU_ID)]
# Fasta file
writeXStringSet(filtered_fasta_minusNC, paste0(otu_output_path, "sequences_minusNC.fasta"))
#taxonomy table
read_tsv('RESULTS/OTUs/OTUstaxonomy.txt') %>%
  filter(OTU_ID %in% OTU_table1[-which(contamdf.prev02$contaminant),]$OTU_ID) %>% 
  arrange(desc(abundance)) %>%
  fwrite(., paste0(otu_output_path, "taxonomy_minusNC.txt"), sep = "\t")
