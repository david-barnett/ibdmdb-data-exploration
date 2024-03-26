library(tidyverse)
library(phyloseq)
library(here)

data_dir <- here("data")

metaphlan3 <- file.path(data_dir, "input", "taxonomic_profiles_3.tsv") %>%
  read_tsv(show_col_types = FALSE) %>%
  rename_with(\(x) str_remove(x, "_profile")) %>%
  rename(Taxon = `Feature\\Sample`) %>%
  filter(Taxon == "UNKNOWN" | grepl("s__", Taxon)) %>%
  separate_wider_delim(
    cols = Taxon, delim = "|",
    names = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
    too_few = "align_end"
  ) %>%
  mutate(
    Taxon = Species,
    across(Kingdom:Genus, \(x) ifelse(is.na(x), yes = "UNKNOWN", no = x))
  )

profiles <- metaphlan3 %>%
  dplyr::select(!Kingdom:Species) %>%
  column_to_rownames("Taxon") %>%
  as.matrix() %>%
  t() %>%
  as.data.frame()

taxonomy <- metaphlan3 %>%
  dplyr::select(Kingdom:Species, Taxon) %>%
  column_to_rownames("Taxon")

metadata <- file.path(data_dir, "input", "hmp2_metadata_2018-08-20.csv") %>%
  read_csv(name_repair = "universal_quiet", show_col_types = FALSE) %>%
  filter(
    data_type %in% c("metagenomics", "metabolomics"),
    External.ID %in% rownames(profiles)
  ) %>%
  select(!where(\(x) all(is.na(x)))) %>% # not all NA
  select(!where(\(x) is.character(x) && length(unique(x)) == 1)) # uninformative

metagenomics_metadata <- metadata %>%
  filter(data_type == "metagenomics") %>%
  mutate(sample_ID = External.ID) %>%
  column_to_rownames("External.ID")

ps <- phyloseq(
  otu_table(profiles, taxa_are_rows = FALSE),
  tax_table(as.matrix(taxonomy)),
  sample_data(metagenomics_metadata)
)

write_rds(ps, file = file.path(data_dir, "phyloseq.rds"))
