xxx <-read.tsv("Users/40085784/Desktop/SRR873610hungate_alignment.primaryalignment.bam.single_assembly_stats.txt")
getwd()



if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("KEGGREST")

library(KEGGREST)
library(dplyr)


library(KEGGREST)

# Load necessary library
library(KEGGREST)

# Function to convert UniProt IDs to KEGG IDs and then to KO numbers
convert_uniprot_to_ko <- function(uniprot_ids) {
  ko_numbers <- list()
  
  for (uniprot_id in uniprot_ids) {
    # Convert UniProt ID to KEGG ID
    kegg_id_result <- keggConv("genes", paste("uniprot:", uniprot_id, sep=""))
    if (length(kegg_id_result) > 0) {
      kegg_id <- names(kegg_id_result)
      
      # Retrieve KO number using KEGG ID
      ko_result <- keggLink("ko", kegg_id)
      if (length(ko_result) > 0) {
        ko_number <- names(ko_result)
        ko_numbers[[uniprot_id]] <- ko_number
      } else {
        ko_numbers[[uniprot_id]] <- NA
      }
    } else {
      ko_numbers[[uniprot_id]] <- NA
    }
  }
  
  return(ko_numbers)
}

# Example usage
uniprot_ids <- c("P12345", "Q9Y2X3")
ko_numbers <- convert_uniprot_to_ko(uniprot_ids)
print(ko_numbers)


# Assuming data is already a dataframe
colnames(data) <- c("V1", "V2", "V3", "V4", "V5")

# Display the dataframe to verify
print(data)


colnames(xgenes_with_uniprot)  <-c("X","Y","Z","UNIPROTID")


# Load necessary library
library(KEGGREST)

# Function to convert UniProt IDs to KEGG IDs and then to KO numbers
convert_uniprot_to_ko <- function(uniprot_ids) {
  ko_numbers <- list()
  
  for (uniprot_id in uniprot_ids) {
    try({
      # Convert UniProt ID to KEGG ID
      kegg_id_result <- keggConv("genes", paste("uniprot:", uniprot_id, sep=""))
      if (length(kegg_id_result) > 0) {
        kegg_id <- names(kegg_id_result)[1]
        
        # Retrieve KO number using KEGG ID
        ko_result <- keggLink("ko", kegg_id)
        if (length(ko_result) > 0) {
          ko_number <- names(ko_result)
          ko_numbers[[uniprot_id]] <- ko_number
        } else {
          ko_numbers[[uniprot_id]] <- NA
        }
      } else {
        ko_numbers[[uniprot_id]] <- NA
      }
    }, silent = TRUE)  # continue on error
  }
  
  return(ko_numbers)
}

# Example usage
uniprot_ids <- c("P12345", "Q9Y2X3")
ko_numbers <- convert_uniprot_to_ko(uniprot_ids)
print(ko_numbers)










# Load necessary library
library(KEGGREST)

# Function to convert UniProt IDs to KO numbers
convert_uniprot_to_ko <- function(uniprot_ids) {
  ko_numbers <- sapply(uniprot_ids, function(id) {
    tryCatch({
      # Convert UniProt ID to KEGG ID
      kegg_id_result <- keggConv("genes", paste("uniprot:", id, sep=""))
      if (length(kegg_id_result) > 0) {
        kegg_id <- names(kegg_id_result)[1]
        
        # Retrieve KO number using KEGG ID
        ko_result <- keggLink("ko", kegg_id)
        if (length(ko_result) > 0) {
          ko_number <- sub("ko:", "", names(ko_result)[1])
          return(ko_number)
        } else {
          return(NA)
        }
      } else {
        return(NA)
      }
    }, error = function(e) {
      return(NA)
    })
  })
  
  return(ko_numbers)
}

# Extract UniProt IDs from the dataframe
uniprot_ids <- xgenes_with_uniprot$UNIPROTID

# Convert UniProt IDs to KO numbers
ko_numbers <- convert_uniprot_to_ko(uniprot_ids)

# Add KO numbers to the dataframe
data$KO <- ko_numbers

# Display the updated dataframe
print(data)




save.image("my_workspace.RData")