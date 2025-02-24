
#Load packages
suppressPackageStartupMessages(library(rmarkdown))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(vegan))
suppressPackageStartupMessages(library(phyloseq))
suppressPackageStartupMessages(library(ggplot2))

setwd("/Users/jasper/Documents/GradSchool/Y3/EHICollaboration/vulture_metagenomics")
#Load Data
#The goal will be to sort data based on metadata information (Remove jejunum and cecum samples, sort by species etc)
metadata_table = read_tsv("/Users/jasper/Documents/GradSchool/Y3/EHICollaboration/vulture_metagenomics/data/DMB0117_metadata.tsv")
count_table = read_tsv("/Users/jasper/Documents/GradSchool/Y3/EHICollaboration/vulture_metagenomics/data/DMB0117_counts.tsv")
mags_table = read_tsv("/Users/jasper/Documents/GradSchool/Y3/EHICollaboration/vulture_metagenomics/data/DMB0117_mag_info.tsv")
coverage_table = read_tsv("/Users/jasper/Documents/GradSchool/Y3/EHICollaboration/vulture_metagenomics/data/DMB0117_coverage.tsv")
#kegg_table = read_tsv("/Users/jasper/Documents/GradSchool/Y3/EHICollaboration/vulture_metagenomics/data/DMB0117_merged_kegg.tsv")
tree_table = read_tsv("/Users/jasper/Documents/GradSchool/Y3/EHICollaboration/vulture_metagenomics/data/DMB0117.tree")

#Convert to dataframe
count_data = as.data.frame(count_table)
mags_data = as.data.frame(mags_table)
coverage_data = as.data.frame(coverage_table)

#Create matrix for metadata sorted by species
metadata_sort = metadata_table[order(as.matrix(metadata_table[,"species"])),]
#Add state data to metadata
metadata_sort$state = substr(word(metadata_sort$region),1,nchar(word(metadata_sort$region))-1)
#Create matrix of count data to be sorted like metadata_sort
count_sort = count_data[,match( t(as_tibble(metadata_sort[,1])), colnames(count_data) ) ]
#Above removes the bins, so they must be added back in
rownames(count_sort) = count_data[,1]
#Transpose data
count_sort = t(count_sort)
#Same for coverage data
coverage_sort = coverage_data[,match( t(as_tibble(metadata_sort[,1])), colnames(coverage_data))]
coverage_sort = t(coverage_sort)
#Same for MAGs
mags_sort = mags_data[match( colnames(count_sort), mags_data$genome),]


#Export filtered data
write.csv(x = metadata_sort, file = "metadata_sort.tsv", sep = "\t")
write.csv(x = count_sort, file = "count_sort.tsv", sep = "\t")
write.csv(x = coverage_sort, file = "coverage_sort.tsv", sep = "\t")
write.table(x = mags_sort, file = "mags_sort.tsv", sep = "\t")

#Order bins
mags_ordered = mags_table[match(colnames(count_sort), mags_table$genome),]


#Tidy data
Samples=rownames(count_sort)
count_mags = count_sort
colnames(count_mags) = mags_ordered$family

tidied = pivot_longer(data = cbind(Samples, as.data.frame(count_sort)), cols = -Samples, names_to = "family", values_to = "count")
#tidied = as.data.frame(tidied) %>% mutate(phylum = str_replace_all(phylum, "d__Bacteria;p__", "")) %>% mutate(phylum = str_replace(phylum, "d__Bacteria;__", "UnclassifiedBacteria")) %>% mutate(phylum = str_replace_all(phylum, "d__Archaea;p__", "")) %>% mutate(phylum = str_replace(phylum, "d__Archaea;__", "UnclassifiedArchaea")) %>% mutate(Samples = str_replace_all(Samples, "1766-..-", ""))
tidied$family = rep(mags_ordered$family,times=nrow(tidied)/nrow(mags_ordered))


#Find relative abundances
rel.count = tidied %>% group_by(Samples) %>% mutate(rel_abund = count/sum(count))
rel.count %>% ggplot(aes(x = Samples, y = rel_abund, fill = family), xaxt="n") + geom_col() + 
  scale_fill_discrete(name=NULL) + 
  labs(x="Sample", y= "Relative Abundance (%)") + 
  theme_classic() + theme(legend.key.size = unit(5, "pt")) + geom_vline(xintercept = 50.5, "solid", "black", 1.5) + theme(axis.text.x = element_blank())
#Left of line is Cathartes aura
#Right of line is Coragyps atratus

#PCoA Bray-Curtis
count_dist <- vegdist(count_sort, method="bray")
dispersion <- betadisper(count_dist, group = metadata_sort$state_sex)
permutest(dispersion)
#Above shows close to difference between species 0.063
#High difference between states (not considering species) 0.001
#Sex not significant 0.788
#Sample Type high significance 0.008
plot(dispersion, hull=FALSE, ellipse=TRUE) ##sd ellipse

#Combine metadata columns
metadata_sort$species_sample = paste(metadata_sort$species, metadata_sort$sample_type, sep=", ")
#p=0.093
metadata_sort$species_state = paste(metadata_sort$species, metadata_sort$state, sep=", ")
#p=0.002
metadata_sort$species_sex = paste(metadata_sort$species, metadata_sort$sex, sep=", ")
#p=0.215
metadata_sort$sample_state = paste(metadata_sort$sample_type, metadata_sort$state, sep=", ")
#p=0.12
metadata_sort$sample_sex = paste(metadata_sort$sample_type, metadata_sort$sex, sep=", ")
#p=0.031
metadata_sort$state_sex = paste(metadata_sort$state, metadata_sort$sex, sep=", ")
#p=0.011





