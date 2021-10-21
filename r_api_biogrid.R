# Load libraries
library(httr)
library(jsonlite)
library(tidyverse)
library(UpSetR)
library(enrichR)

# Set up working directory
setwd("~/GitHub/lqts-ppi")

# Set up urls and access key ####
base_url <- 'https://webservice.thebiogrid.org'

request_inter_url <- paste0(base_url, '/interactions/')

request_org_url <- paste0(base_url, '/organisms/')

ACCESS_KEY <- '1d175e1ea79a720986c1e6228f01d738'


# Get list of supported organisms in BioGRiD ####
orgList <- GET(request_org_url, 
           query = list(accesskey = ACCESS_KEY, format = 'json'))

# Convert json object to R object
orgRList <- fromJSON(rawToChar(orgList$content))

# Create a data frame from json-to-R converted object
orgDF <- tibble::enframe(unlist(orgRList), 
                         name = 'taxId', value = 'species')

# Add gene(s) of interest ####
geneList <- c('ANK2', 'KCNQ1', 'KCNH2', 'SCN5A', 'KCNE1', 'KCNE2', 'KCNJ2', 'CACNA1C', 'CAV3', 'SCN4B', 'AKAP9', 'SNTA1', 'ALG10')

species_id <- c(orgDF$taxId[orgDF$species == 'Homo sapiens'],
                orgDF$taxId[orgDF$species == 'Rattus norvegicus'],
                orgDF$taxId[orgDF$species == 'Mus musculus'])

evidenceList = c("POSITIVE GENETIC", "PHENOTYPIC ENHANCEMENT")

# Setting up parameters ####
params <- list(
  accesskey = ACCESS_KEY,
  format = 'json', # Return results in TAB2 format
  geneList = paste(geneList, collapse = '|'), # Must be | separated
  searchNames = 'true', # Search against official names
  includeInteractors = 'true', # Set to true to get any interaction involving EITHER gene, set to false to get interactions between genes
  includeInteractorInteractions = 'false',
  taxId = paste(species_id, collapse = '|'),
  evidenceList = paste(evidenceList, collapse = '|'), # evidence parameters to exclude or include, see below
  includeEvidence = 'false' # # If false "evidenceList" is evidence to exclude, if true "evidenceList" is evidence to show
)

intList <- GET(request_inter_url,
               query = params)

intRList <- fromJSON(rawToChar(intList$content))

intDF <- tibble::enframe(unlist(intRList))

# for loop to get a tidy data frame of interactions ####
df_aggregate <- c()

for (interaction in 1:length(intRList)) {
  df_interaction <- tibble::enframe(intRList[[interaction]])
  df_aggregate <- rbind(pivot_wider(df_interaction,
                                    names_from = name,
                                    values_from = value),
                        df_aggregate)
  # df_interaction <- pivot_wider(df_interaction,
  #                               names_from = name,
  #                               values_from = value)
  # df_aggregate[nrow(df_aggregate) + 1,] <- df_interaction
}


# Extracting inner list elements as data frames ####
# df <- intRList[[1]]
# 
# test <- tibble::enframe(unlist(intRList[[1]]))
# 
# df <- tibble::enframe(unlist(df))
# 
# df_wide <- pivot_wider(df, names_from = name, values_from = value)
# 
# df2 <- intRList[[2]]
# 
# df2 <- tibble::enframe(unlist(df2))
# 
# df2_wide <- pivot_wider(df2, names_from = name, values_from = value)
# 
# df12_wide <- rbind(df_wide, df2_wide)


# Making QUANTITATION numeric ####
df_aggregate$QUANTITATION <- as.numeric(df_aggregate$QUANTITATION)

# df_aggregate$OFFICIAL_SYMBOL_A <- as.character(df_aggregate$OFFICIAL_SYMBOL_A)
# 
# df_aggregate$OFFICIAL_SYMBOL_B <- as.character(df_aggregate$OFFICIAL_SYMBOL_B)

# Creating a list of biogrid interactions ####
biogrid <- list()

for (gene in geneList) {
  biogrid[length(biogrid) + 1] <- list(na.omit(
    case_when(
      toupper(df_aggregate$OFFICIAL_SYMBOL_A) == gene ~ toupper(df_aggregate$OFFICIAL_SYMBOL_B),
      toupper(df_aggregate$OFFICIAL_SYMBOL_B) == gene ~ toupper(df_aggregate$OFFICIAL_SYMBOL_A)
      )
  )
  )
}

# Name lists ####
names(biogrid) <- geneList

# biogrid <- lapply(biogrid, function(biogrid) biogrid[!is.na(biogrid)])

# Keep a siongle record of each interaction #####
biogrid <- lapply(biogrid, function(biogrid) unique(biogrid))

# Checking if gene is within interactor list and removing from the list ####
for (gene in geneList) {
  if (gene %in% biogrid[[gene]]){
    biogrid[[gene]] <- biogrid[[gene]][-which(biogrid[[gene]] == gene)]
  }
}

# Check if the above chunk was successful
# for (gene in geneList) {
#   print(gene %in% biogrid[[gene]])
# }

# Create UpSet plot of common interactors
upset_plot <- upset(fromList(biogrid), nsets = 13, order.by = "freq", nintersects = NA)

# Enrichr analysis with enrichR ####

setEnrichrSite("Enrichr") # Human genes
websiteLive <- TRUE

dbs <- listEnrichrDbs()

if (is.null(dbs)) websiteLive <- FALSE
if (websiteLive) head(dbs)

dbs <- c("GO_Molecular_Function_2021", "GO_Cellular_Component_2021", "GO_Biological_Process_2021")
if (websiteLive) {
  ank2_go_enriched <- enrichr(biogrid[["ANK2"]], dbs)
}

dbs_kegg <- c("KEGG_2021_Human")
if (websiteLive) {
  ank2_enriched_kegg <- enrichr(biogrid[["ANK2"]], dbs_kegg)
}

if (websiteLive) enriched[["GO_Molecular_Function_2015"]]

ank2_go_mol_plot <- if (websiteLive) plotEnrich(ank2_go_enriched[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")

ank2_go_com_plot <- if (websiteLive) plotEnrich(ank2_go_enriched[[2]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")

ank2_go_pro_plot <- if (websiteLive) plotEnrich(ank2_go_enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")

# KEGG
ank2_kegg_plot <- if (websiteLive) plotEnrich(ank2_enriched_kegg[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")

# KCNQ1 ####
if (websiteLive) {
  kcnq1_go_enriched <- enrichr(biogrid[["KCNQ1"]], dbs)
}

kcnq1_go_mol_plot <- if (websiteLive) plotEnrich(kcnq1_go_enriched[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")

kcnq1_go_com_plot <- if (websiteLive) plotEnrich(kcnq1_go_enriched[[2]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")

kcnq1_go_pro_plot <- if (websiteLive) plotEnrich(kcnq1_go_enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")

if (websiteLive) {
  kcnq1_enriched_kegg <- enrichr(biogrid[["KCNQ1"]], dbs_kegg)
}

kcnq1_kegg_plot <- if (websiteLive) plotEnrich(kcnq1_enriched_kegg[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")

cowplot::plot_grid(ank2_go_com_plot, kcnq1_go_com_plot)

cowplot::plot_grid(ank2_go_mol_plot, kcnq1_go_mol_plot)

cowplot::plot_grid(ank2_go_pro_plot, kcnq1_go_pro_plot)
