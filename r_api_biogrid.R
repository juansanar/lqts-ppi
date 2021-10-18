# Load libraries
library(httr)
library(jsonlite)
library(tidyverse)
library(magrittr)
library(UpSetR)

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

# Create UpSet plot
upset(fromList(biogrid), nsets = 13, order.by = "freq", nintersects = NA)

# Finding common interactors ####

