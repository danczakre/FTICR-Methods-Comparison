### Testing compound name conversion using the KEGG REST API

library(readxl)
library(rvest)
library(dplyr)

output_prefix = "Output"

# Set working directory
setwd("~/Documents/East River/New Comparative Analyses/NonSPE_Neg/")

# Load in molecular file from ftmsRanalysis
mol = read.csv("Processed_NonSPE_Neg_NoP_Mol.csv")
mol = mol[!is.na(mol$MolForm),]
mol = mol[!duplicated(mol$MolForm),]

# Loop through all of the compounds and use the REST API
mf.table = NULL # Empty object to store information

for(i in 1:nrow(mol)){
  mf = mol$MolForm[i] # Set current KO number
  page = paste0("http://rest.kegg.jp/find/compound/", mf, "/formula") # Find entry on the API

  if(length(grep(" ", page)) > 0){
    page = gsub(" ", "%20", page)
  } # Need to fix spaces in page name
  
  page = read_html(page) # Load the HTML page
  
  if(length(page) == 1){
    page = data.frame(CPD = NA, KEGG_MF = NA, Detected_MF = mf)
  } else {
    page = page %>% html_node("body") %>% html_children() %>% html_text()
    page = gsub("'", "", page)
    page = read.table(text = page, sep = "\t", header = F)
    
    if((length(which(tolower(page$V2) == tolower(mf))) == 0)){
      page = data.frame(CPD = NA, KEGG_MF = NA, Detected_MF = mf)
    } else {
      mf.loc = which(tolower(page$V2) == tolower(mf)) # Finding exact location of compound
      page = data.frame(CPD = page$V1[mf.loc], KEGG_MF = page$V2[mf.loc], Detected_MF = mf) # REporting the CPD of the exact match
    } # This loop captures misidentified compounds (i.e., something with an identical partial name)
}
  
  mf.table = rbind(mf.table, page)
  
  print(paste("Just finished with", mf, "which is entry", which(mol$MolForm %in% mf)))
  
  rm(df, page, split.page, mf.loc, j)
} # Loop through each KO; need the number variant because it occasionally times out

# Cleaning up duplicates and missing values
mf.table = mf.table[!is.na(mf.table$CPD),]

# Writing results out
write.table(mf.table, paste0(output_prefix, "_KEGG_CPD_FTICR_Matching.txt"), sep = "\t", quote = F, row.names = F)
