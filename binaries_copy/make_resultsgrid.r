#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)
library(dplyr)
library(dbplyr)
library(ggplot2)
library(RSQLite)
library(tidyr)
library(tidyverse)
library(phytools)
library(ape)
library(treeio)
library(ggtree)
library(svglite)
sqlite <- dbDriver("SQLite")
mydb <- dbConnect(sqlite, '~/Documents/Complex_Portal/ComplexPortalSearch2.db')
dbListTables(mydb)
results <- dbSendQuery(mydb, "SELECT `taxid`, `#Complex ac` FROM Complex_portal_components
                         INNER JOIN Complex_portal_Results
                         ON Complex_portal_components.'UniprotID' = Complex_portal_Results.'query'
                         ORDER BY taxid, '#Complex ac'")
#results <- dbBind(results, data.frame()
data = dbFetch(results)
df <- data %>%
  add_column(add_column = 100)
resultsgrid <- as.data.frame(df %>% 
                               pivot_wider(
                                 names_from = `#Complex ac`, 
                                 values_from = `add_column`,
                                 values_fn = min
                               ))
write.csv(resultsgrid, "~/Documents/Complex_Portal/resultsgrid_test.csv", row.names = FALSE)
