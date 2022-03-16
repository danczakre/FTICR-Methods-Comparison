### Testing compound name conversion using the KEGG REST API

library(rvest)
library(dplyr)

dataset = "NonSPE_Pos"

# Load data
setwd("~/Documents/East River/New Comparative Analyses/NonSPE_Neg/")
comp.tab = read.table(list.files(pattern = "FTICR_Matching"), sep = "\t", comment.char = "", header = 1)

# Look up compound information on KEGG
path.tab = NULL

for(curr.cpd in comp.tab$CPD){
  page = paste0("http://rest.kegg.jp/get/", curr.cpd)
  page = read_html(page) # Load the HTML page
  page = page %>% html_node("body") %>% html_children() %>% html_text()
  page = str_extract_all(page, ".*")
  page = page[[1]][-which(page[[1]] %in% "")]
  
  if(length(grep("PATHWAY", page)) == 0){
    page = data.frame(CPD_ID = curr.cpd, Map_ID = "Unavailable", Pathway_Name = "Unavailable")
    path.tab = rbind(path.tab, page)
  } else {
    if(length(grep("ENZYME", page)) == 0){
      if(length(grep("MODULE", page)) == 1){
        page = page[grep("PATHWAY", page):(grep("MODULE", page)-1)]
      } else if(length(grep("BRITE", page)) == 1){
        page = page[grep("PATHWAY", page):(grep("BRITE", page)-1)]
      } else {
        page = page[grep("PATHWAY", page):(grep("DBLINKS", page)-1)]
      } # If both BRITE and ENZYME are missing, DBLINKS is next
    } else if(length(grep("MODULE", page)) == 1){
      page = page[grep("PATHWAY", page):(grep("MODULE", page)-1)]
    } else {
      page = page[grep("PATHWAY", page):(grep("ENZYME", page)-1)]
    } # This if-loop is necessary if the entry doesn't have an associated enzyme
    
    page = gsub("PATHWAY", "", page)
    page = gsub("^.*map", "map", page)
    page = data.frame(Map_ID = str_extract(page, "map[0-9]{5}"), Pathway_Name = gsub("map[0-9]{5}  ", "", page))
    page = data.frame(CPD_ID = curr.cpd, page)
    
    path.tab = rbind(path.tab, page)
  } # This if-loop controls whether the pathway exists at all
  
  print(paste0("Finished processing compound ", curr.cpd, " which is entry ", which(comp.tab$CPD %in% curr.cpd),
              " (", round((which(comp.tab$CPD %in% curr.cpd)/length(comp.tab$CPD))*100, digits = 2), "%)"))
  
  rm(page)
}

rm(curr.cpd)

# Writing results
write.table(path.tab, 
            paste0(dataset, "_Pathways_from_Compounds.txt"), 
            quote = F, row.names = F, sep = "\t")
