### Testing compound name conversion using the KEGG REST API

library(readxl)
library(rvest)
library(dplyr)

# Set working directory
setwd("~/Documents/MSC-2 Chitin Degradation/Transcript Mapping/")

# Load in annotations
annotations = read.csv("../annotations_alt.csv")

# Loop through all of the compounds and use the REST API
brite.table = NULL # Empty object to store information

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

# Working with the CPD values
incidence = data.frame(Detected_MF = mol$MolForm, Incidence = rowSums(data)/ncol(data))
merged.cpd = mf.table[,c("Detected_MF", "CPD")] %>% 
  left_join(incidence, by = "Detected_MF")
col = colorRampPalette(c("dodgerblue2", "gray90", "firebrick2"))
col = data.frame(Incidence = unique(merged.cpd$Incidence[order(merged.cpd$Incidence)]), 
                 Color = col(length(unique(merged.cpd$Incidence))))
merged.cpd = merged.cpd %>% left_join(col, by = "Incidence")
merged.cpd = merged.cpd[,c("CPD", "Color")]
write.table(merged.cpd, 
            paste0(dataset, "_CPD_with_Incidence_Colors.txt"), 
            quote = F, row.names = F, sep = "\t")
