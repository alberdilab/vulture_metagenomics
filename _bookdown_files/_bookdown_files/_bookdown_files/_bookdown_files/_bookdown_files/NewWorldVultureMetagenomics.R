# Base
library(R.utils)
library(knitr)
library(tidyverse)
library(devtools)
library(tinytable)
library(rairtable)

# For tree handling
library(ape)
library(phyloseq)
library(phytools)

# For plotting
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(ggnewscale)
library(gridExtra)
library(ggtreeExtra)
library(ggtree)
library(ggh4x)

# For statistics
library(spaa)
library(vegan)
library(Rtsne)
library(geiger)
library(hilldiv)
library(distill)
library(ANCOMBC)
library(lme4)

library(distillery)
library(hilldiv2)
library(distillR)

#Load Data

#Metadata
sample_metadata <- read_tsv("data/DMB0117_metadata.tsv.gz") %>%
  rename(sample=1)
#Read Counts
read_counts <- read_tsv("data/DMB0117_counts.tsv.gz") %>%
  rename(genome=1)
#Genome base Hits
genome_coverage <- read_tsv("data/DMB0117_coverage.tsv.gz") %>%
  rename(genome=1)
#Genome Taxonomy
genome_metadata <- read_tsv("data/DMB0117_mag_info.tsv.gz") %>%
  rename(length=mag_size)

#Genome Tree
genome_tree <- read_tree("data/DMB0117.tree")
genome_tree$tip.label <- str_replace_all(genome_tree$tip.label,"'", "") #remove single quotes in MAG names
genome_tree <- keep.tip(genome_tree, tip=genome_metadata$genome) # keep only MAG tips

#Genome Annotations (Requires Personal Access Token)
#Need key, until then cannot run anything involving genome_annotations.tsv.xz
#airtable("MAGs", "appWbHBNLE6iAsMRV") %>% #get base ID from Airtable browser URL
# read_airtable(., fields = c("ID","mag_name","number_genes","anno_url"), id_to_col = TRUE) %>% #get 3 columns from MAGs table
#  filter(mag_name %in% paste0(genome_metadata$genome,".fa")) %>% #filter by MAG name
#  filter(number_genes > 0) %>% #genes need to exist
#  select(anno_url) %>% #list MAG annotation urls
#  pull() %>%
#  read_tsv() %>% #load all tables
#  rename(gene=1, genome=2, contig=3) %>% #rename first 3 columns
#  write_tsv(file="data/genome_annotations.tsv.xz") #write to overall compressed file

#genome_annotations <- read_tsv("data/genome_annotations.tsv.xz") %>%
#  rename(gene=1, genome=2, contig=3)

min_coverage=0.3
read_counts_filt <- genome_coverage %>%
  mutate(across(where(is.numeric), ~ ifelse(. > min_coverage, 1, 0))) %>%
  mutate(across(-1, ~ . * read_counts[[cur_column()]])) 


genome_counts_filt %>%
  mutate_at(vars(-genome),~./sum(.)) %>% #apply TSS nornalisation
  pivot_longer(-genome, names_to = "sample", values_to = "count") %>% #reduce to minimum number of columns
  left_join(., genome_metadata, by = join_by(genome == genome)) %>% #append genome metadata
  left_join(., sample_metadata, by = join_by(sample == sample)) %>% #append sample metadata
  filter(!is.na(count)) %>%
  ggplot(aes(y=count,x=sample, fill=phylum, group=phylum)) + #grouping enables keeping the same sorting of taxonomic units
  geom_bar(stat="identity", colour="white", linewidth=0.1) + #plot stacked bars with white borders
  scale_fill_manual(values=phylum_colors) +
  labs(x = "Relative abundance", y ="Samples") +
  facet_nested(. ~ species.y + sample_type ,  scales="free", space="free") + #facet per day and treatment
  scale_y_continuous(expand = c(0.001, 0.001)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.x = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
        legend.position = "none",
        strip.background.x=element_rect(color = NA, fill= "#f4f4f4"))




