library(ggplot2)
library(plyr)
library(tidyverse)
library(cowplot)

setwd('')

#read in EPA taxonomy
EPA_table <- read_tsv('RESULTS/EPA/Microascales/per_query.tsv') %>%
  bind_rows(read_tsv('RESULTS/EPA/Ophiostomatales/per_query.tsv'))
#split taxopath into separate columns
EPA_table$name <- gsub('_.*','',EPA_table$name) 
#bind with dnabarcoder
full_table <- EPA_table %>%
  cbind(str_split_fixed(EPA_table$taxopath, ";", 8)) %>%
  #join with the dnabarcoder taxonomy
  left_join(read_tsv('RESULTS/OTUs/OTUstaxonomy_minusNC.txt'), by = c("name" = "OTU_ID"))

#edit synonymy manually
write_csv(full_table,'RESULTS/EPA_vs_dnabarcoder.csv')
#read back in
full_table <- read_csv('RESULTS/EPA_vs_dnabarcoder.csv')

#Add column for whether dnabarcoder and T-BAS agree
full_table[,c('phylum.match','class.match','order.match','family.match','genus.match','species.match')] <- NA


#dnabarcoder vs TBAS
for (i in 1:nrow(full_table)) {
  #if blast species matches TBAS species then add Y(es) to match column else add N(o)
  if (is.na(grepl(full_table$species[i], full_table$`8`[i], fixed = TRUE))) {
    full_table$species.match[i] <- NA
  } else if (grepl(full_table$species[i], full_table$`8`[i], fixed = TRUE)) {
    full_table$species.match[i] <- 'Y'
  } else {
    full_table$species.match[i] <- 'N'
  }
  #genus level
  if (is.na(grepl(full_table$genus[i], full_table$`7`[i], fixed = TRUE))) {
    full_table$genus.match[i] <- NA
  } else if (grepl(full_table$genus[i], full_table$`7`[i], fixed = TRUE)) {
    full_table$genus.match[i] <- 'Y'
  } else {
    full_table$genus.match[i] <- 'N'
  }
  #family
  if (is.na(grepl(full_table$family[i], full_table$`6`[i], fixed = TRUE))) {
    full_table$family.match[i] <- NA
  } else if (grepl(full_table$family[i], full_table$`6`[i], fixed = TRUE)) {
    full_table$family.match[i] <- 'Y'
  } else {
    full_table$family.match[i] <- 'N'
  }
  #order
  if (is.na(grepl(full_table$order[i], full_table$`5`[i], fixed = TRUE))) {
    full_table$order.match[i] <- NA
  } else if (grepl(full_table$order[i], full_table$`5`[i], fixed = TRUE)) {
    full_table$order.match[i] <- 'Y'
  } else {
    full_table$order.match[i] <- 'N'
  }
  #class
  if (is.na(grepl(full_table$class[i], full_table$`4`[i], fixed = TRUE))) {
    full_table$class.match[i] <- NA
  } else if (grepl(full_table$class[i], full_table$`4`[i], fixed = TRUE)) {
    full_table$class.match[i] <- 'Y'
  } else {
    full_table$class.match[i] <- 'N'
  }
  #phylum
  if (is.na(grepl(full_table$phylum[i], full_table$`3`[i], fixed = TRUE))) {
    full_table$phylum.match[i] <- NA
  } else if (grepl(full_table$phylum[i], full_table$`3`[i], fixed = TRUE)) {
    full_table$phylum.match[i] <- 'Y'
  } else {
    full_table$phylum.match[i] <- 'N'
  }
  
}


#count how OTUs ID agree at each taxon rank
dnabarcoder_vs_EPA <- c()
for(col in colnames(full_table)[c(27:32)]){
  #count how many agree Y and divide the total number where they could make a comparison for that rnak
  dnabarcoder_vs_EPA <- c(dnabarcoder_vs_EPA,length(grep('Y',full_table[[col]]))/sum(!is.na(full_table[[col]])))
}

dnabarcoder_vs_EPA_table <- cbind(tibble(agreement = dnabarcoder_vs_EPA),rank = c('phylum','class','order','family','genus','species'))

(agree_plot <- ggplot(data = dnabarcoder_vs_EPA_table, aes(x = rank, y = agreement*100)) +
    geom_point() +
    geom_line(aes(group=1)) +
    scale_x_discrete(limits = c('phylum','class','order','family','genus','species')) +
    ylab('EPA vs dnabarcoder ID agreement (%)') +
    xlab('Taxon rank') +
    theme_minimal() +
    theme(aspect.ratio=1) +
    ylim(c(0,100))
)


#count how many IDs at each rank for both tools
dnabarcoder_vs_EPA_number <- c()
for(col in c(3:8,'phylum','class','order','family','genus','species')){
  #count how many agree Y and divide the total number where they could make a comparison for that rank
  dnabarcoder_vs_EPA_number <- c(dnabarcoder_vs_EPA_number,sum(!is.na(full_table[[col]])))
}

dnabarcoder_vs_EPA_number_table <- cbind(tibble(number = dnabarcoder_vs_EPA_number),rank = c('phylum','class','order','family','genus','species','phylum','class','order','family','genus','species'), tool = c(rep('EPA',6),rep('dnabarcoder',6)))

(id_plot <- ggplot(data = dnabarcoder_vs_EPA_number_table, aes(x = rank, y = (number/126)*100, col = tool, group = tool)) +
    geom_point() +
    geom_line() +
    scale_x_discrete(limits = c('phylum','class','order','family','genus','species')) +
    ylab('OTUs identified (%)') +
    xlab('Taxon rank') +
    theme_minimal() +
    theme(aspect.ratio=1) +
    ylim(c(0,100)) +
    labs(col = "Tool"))

png("FIGURES/EPA_vs_dnabarcoder.png",  res = 600, width = 250, height = 100, units = 'mm')
plot_grid(id_plot,agree_plot, nrow = 1, align = 'hv', labels = c('a','b'))
dev.off()
