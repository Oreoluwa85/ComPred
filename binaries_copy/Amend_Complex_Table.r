#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)
directory <- getwd()
dat <- as.data.frame(read.csv(file ="~/Documents/Complex_Portal/Complex_Portal_dbfiles/Complex_Portal_Table_raw.tsv", sep = '\t', header = TRUE, fill = TRUE))
col <- dat$'Identifiers..and.stoichiometry..of.molecules.in.complex'
dat_changed <- as.data.frame(dat %>%
  mutate(UniprotID = strsplit(as.character(col), "\\|")) %>%
  unnest(UniprotID) %>%
  filter(UniprotID != "") %>%
  select(1:4,UniprotID, 6:19))
write.csv(dat_changed, paste(directory,"/Complex_Portal_Table.tsv",sep=""), row.names = FALSE)
