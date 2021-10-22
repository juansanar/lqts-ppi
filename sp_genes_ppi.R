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
geneList_sp <- c('MYH7', 'MYBPC3', 'TNNT2', 'TNNI3', 'MYL2', 'MYL3', 'ACTC1', 'TPM1', 'TNNC1', 'MYH6', 'CSRP3', 'DES', 'TCAP', 'PDLIM3', 'PLN', 'LDB3', 'LMNA', 'VCL', 'RBM20', 'TTN')

species_id <- c(orgDF$taxId[orgDF$species == 'Homo sapiens'],
                orgDF$taxId[orgDF$species == 'Rattus norvegicus'],
                orgDF$taxId[orgDF$species == 'Mus musculus'])

evidenceList = c("POSITIVE GENETIC", "PHENOTYPIC ENHANCEMENT")

# Setting up parameters ####
params <- list(
  accesskey = ACCESS_KEY,
  format = 'json', # Return results in TAB2 format
  geneList = paste(geneList_sp, collapse = '|'), # Must be | separated
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
df_aggregate_sp <- c()

for (interaction in 1:length(intRList)) {
  df_interaction <- tibble::enframe(intRList[[interaction]])
  df_aggregate_sp <- rbind(pivot_wider(df_interaction,
                                    names_from = name,
                                    values_from = value),
                        df_aggregate_sp)
  # df_interaction <- pivot_wider(df_interaction,
  #                               names_from = name,
  #                               values_from = value)
  # df_aggregate_sp[nrow(df_aggregate_sp) + 1,] <- df_interaction
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
df_aggregate_sp$QUANTITATION <- as.numeric(df_aggregate_sp$QUANTITATION)

# df_aggregate_sp$OFFICIAL_SYMBOL_A <- as.character(df_aggregate_sp$OFFICIAL_SYMBOL_A)
# 
# df_aggregate_sp$OFFICIAL_SYMBOL_B <- as.character(df_aggregate_sp$OFFICIAL_SYMBOL_B)

# Creating a list of biogrid interactions ####
biogrid_sp <- list()

for (gene in geneList_sp) {
  biogrid_sp[length(biogrid_sp) + 1] <- list(na.omit(
    case_when(
      toupper(df_aggregate_sp$OFFICIAL_SYMBOL_A) == gene ~ toupper(df_aggregate_sp$OFFICIAL_SYMBOL_B),
      toupper(df_aggregate_sp$OFFICIAL_SYMBOL_B) == gene ~ toupper(df_aggregate_sp$OFFICIAL_SYMBOL_A)
    )
  )
  )
}

# Name lists ####
names(biogrid_sp) <- geneList_sp

# biogrid <- lapply(biogrid, function(biogrid) biogrid[!is.na(biogrid)])

# Keep a siongle record of each interaction #####
biogrid_sp <- lapply(biogrid_sp, function(biogrid_sp) unique(biogrid_sp))

# Checking if gene is within interactor list and removing from the list ####
for (gene in geneList_sp) {
  if (gene %in% biogrid_sp[[gene]]){
    biogrid_sp[[gene]] <- biogrid_sp[[gene]][-which(biogrid_sp[[gene]] == gene)]
  }
}

# Check if the above chunk was successful
# for (gene in geneList_sp) {
#   print(gene %in% biogrid[[gene]])
# }

# Create UpSet plot of common interactors
upset_plot <- upset(fromList(biogrid_sp), nsets = 13, order.by = "freq", nintersects = NA)