#PIME
library(tidyverse)
library(vegan)
library(phyloseq)
library(pime)


setwd("RESULTS/")

#read in dissimilarity matrices
OTU_table <- read_tsv('OTUs/OTUsotu_table_minusNC.txt') %>% 
  column_to_rownames(var = "OTU_ID")

# Calculate relative abundance
rel_abund <- sweep(OTU_table, 2, colSums(OTU_table), FUN = "/")
# Apply 0.1% relative abundance threshold
OTU_table[rel_abund <= 0.001] <- 0
# Remove OTUs that are now 0 across all samples
nonzero_rows <- rowSums(OTU_table) > 0
OTU_table <- OTU_table[nonzero_rows, ]

#transpose
OTU_table <- OTU_table %>%
  t() %>%
  as.data.frame()

#read in sample data
sample_data <- read_csv('../sample_metadata.csv') %>%
  dplyr::select(c(`Sequence-ID`,species,country)) %>% 
  na.omit() %>%
  filter(`Sequence-ID` %in% rownames(OTU_table)) %>%
  arrange(`Sequence-ID`)
#change all non I. typographus to UK native
for(i in 1:nrow(sample_data)){
  if(sample_data$species[i] != 'Ips typographus'){
    sample_data$species[i] <- 'UKnative'
  }
}


#find out which samples dont have metadata for and remove them
missing_metadata <- setdiff(rownames(OTU_table), sample_data$`Sequence-ID`)
OTU_table_reduced <- OTU_table[ ! rownames(OTU_table) %in% missing_metadata, ] %>%
  mutate_if(is.numeric, ~1 * (. > 0))


#link data into phyloseq object
taxonomy <- read_tsv('OTUs/OTUstaxonomy_minusNC.txt') %>%
  select(-c(ReferenceID, rank, score, cutoff, abundance)) %>%
  column_to_rownames('OTU_ID')

ps <- phyloseq(otu_table(OTU_table_reduced, taxa_are_rows=FALSE), 
               tax_table(as.matrix(taxonomy)))
sample_metadata <- phyloseq::sample_data(sample_data %>% column_to_rownames('Sequence-ID'))
physeq <- merge_phyloseq(ps,sample_metadata)



###########################
# Run PIME
#estimate OOB error rate
source('../../ANALYSIS_SCRIPTS/pime_no_prev_filter.R')
pime.oob.error(physeq, "species")
# 0

#split tables by species
per_variable_obj <- pime.split.by.variable(physeq, "species")
per_variable_obj


#calculates prevalence for taxa at 0% prev. filtering
prevalences <- pime.prevalence_noprev_filter(per_variable_obj)
prevalences

#get best prevalence value
set.seed(8)
best.prev <- pime.best.prevalence_noprev_filter(prevalences, "species")

#add column showing whether that OTU contributes to defining Ips or Native
imp_table <- imp$Prevalence
imp_table$indicator_species <- NULL

for(i in 1:nrow(imp_table)){
  if(imp_table[i,'Ips.typographus']>imp_table[i,'UKnative']){
    imp_table[i,'indicator_species'] = 'UKnative'
  } else{
    imp_table[i,'indicator_species'] = 'Ips.typographus'
  }
}

write_csv(imp_table, 'PIME_Results_IpstypographusvsUKnative.csv')
head(imp)


#likelihood of introducing bias
#randomises sample labels into arbitrary groups and calculates OOB. Used to check whether diff. in original groups due to chance
randomized <- pime.error.prediction(physeq, "species", bootstrap = 100, parallel = FALSE, max.prev = 10)
randomized$Plot
randomized$`Results table`

#replicate without randomisation to show how consistent results are
replicated.oob.error <- pime.oob.replicate(prevalences, "species", bootstrap = 100, parallel = TRUE)
replicated.oob.error$Plot
replicated.oob.error$`Results table`

#show random and real together
png('FIGURES/PIME_prediction_error.png', res = 400, height = 2000, width = 2000)
ggplot(data = randomized$`Results table`, 
       mapping = aes(x='5', y = `Prevalence5%`)) +
  geom_boxplot(aes(col = "Randomised")) +
  geom_boxplot(data = replicated.oob.error$`Results table`, aes(col = 'Original')) +
  ylim(0,1) +
  xlab('Prevalence %') +
  ylab('OOB error') +
  theme_bw() +
  coord_fixed() +
  scale_color_manual(name = "Dataset",
                     values = c("Randomised" = "black", "Original" = "red"))
dev.off()
