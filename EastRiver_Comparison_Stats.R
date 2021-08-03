#### East River ESI comparative analyses statistics

library(vegan)
library(reshape2)
library(ggplot2)
library(ggthemes)
library(stringr)
library(dplyr)

# Paired stats switch
pair = F

# Load in data
setwd("~/Documents/East River/New Comparative Analyses/")
sneg.data = read.csv("SPE_Neg/Processed_SPE_Neg_NoP_Data.csv", row.names = 1)
sneg.mol = read.csv("SPE_Neg/Processed_SPE_Neg_NoP_Mol.csv", row.names = 1)

spos.data = read.csv("SPE_Pos/Processed_SPE_Pos_NoP_Data.csv", row.names = 1)
spos.mol = read.csv("SPE_Pos/Processed_SPE_Pos_NoP_Mol.csv", row.names = 1)

nneg.data = read.csv("NonSPE_Neg/Processed_NonSPE_Neg_NoP_Data.csv", row.names = 1)
nneg.mol = read.csv("NonSPE_Neg/Processed_NonSPE_Neg_NoP_Mol.csv", row.names = 1)

npos.data = read.csv("NonSPE_Pos/Processed_NonSPE_Pos_NoP_Data.csv", row.names = 1)
npos.mol = read.csv("NonSPE_Pos/Processed_NonSPE_Pos_NoP_Mol.csv", row.names = 1)


# ####################### #
#### Clean-up and stuff ####
# ####################### #

# Selecting peaks with molecular formulas
sneg.data = sneg.data[!is.na(sneg.mol$MolForm),]
sneg.mol = sneg.mol[!is.na(sneg.mol$MolForm),]

spos.data = spos.data[!is.na(spos.mol$MolForm),]
spos.mol = spos.mol[!is.na(spos.mol$MolForm),]

nneg.data = nneg.data[!is.na(nneg.mol$MolForm),]
nneg.mol = nneg.mol[!is.na(nneg.mol$MolForm),]

npos.data = npos.data[!is.na(npos.mol$MolForm),]
npos.mol = npos.mol[!is.na(npos.mol$MolForm),]

# Removing duplicate molforms
sneg.data = sneg.data[!duplicated(sneg.mol$MolForm),]
sneg.mol = sneg.mol[!duplicated(sneg.mol$MolForm),]

spos.data = spos.data[!duplicated(spos.mol$MolForm),]
spos.mol = spos.mol[!duplicated(spos.mol$MolForm),]

nneg.data = nneg.data[!duplicated(nneg.mol$MolForm),]
nneg.mol = nneg.mol[!duplicated(nneg.mol$MolForm),]

npos.data = npos.data[!duplicated(npos.mol$MolForm),]
npos.mol = npos.mol[!duplicated(npos.mol$MolForm),]

# Reordering data
sneg.data = sneg.data[,order(str_extract(colnames(sneg.data), "_[0-9]{3}_"))]
spos.data = spos.data[,order(str_extract(colnames(spos.data), "_[0-9]{3}_"))]
nneg.data = nneg.data[,order(str_extract(colnames(nneg.data), "_[0-9]{3}_"))]
npos.data = npos.data[,order(str_extract(colnames(npos.data), "_[0-9]{3}_"))]

# Merging data into lists
data = list(sneg.data, spos.data, nneg.data, npos.data)
names(data) = c("SPE_Neg", "SPE_Pos", "Non_Neg", "Non_Pos")

mol = list(sneg.mol, spos.mol, nneg.mol, npos.mol)
names(mol) = c("SPE_Neg", "SPE_Pos", "Non_Neg", "Non_Pos")

# Readding MolForm and combining data into one object
sneg.data = cbind(MolForm = sneg.mol$MolForm, sneg.data)
spos.data = cbind(MolForm = spos.mol$MolForm, spos.data)
nneg.data = cbind(MolForm = nneg.mol$MolForm, nneg.data)
npos.data = cbind(MolForm = npos.mol$MolForm, npos.data)

mol_data = sneg.data %>% full_join(spos.data, by = "MolForm") %>% full_join(nneg.data, by = "MolForm") %>%
  full_join(npos.data, by = "MolForm")
row.names(mol_data) = mol_data$MolForm; mol_data = mol_data[,-1]
mol_data[is.na(mol_data)] = 0
mol_data[mol_data > 0] = 1

rm(list = ls(pattern = "\\.data|\\.mol"))


# ####################### #
#### Calculating stuff #### 
# ####################### #

### Characteristics stats
char = NULL

for(a in 1:length(data)){
  # Select data
  temp.data = data[[a]]
  temp.mol = mol[[a]]
  
  # Empty object
  temp.char = data.frame(NOSC = rep(NA, ncol(temp.data)), DBE = NA, AI_Mod = NA, NtoC = NA,
                         row.names = colnames(temp.data))
  
  # Averaging by sample
  for(i in 1:nrow(temp.char)){
    temp = temp.mol[which(temp.data[,i] > 0),]
    
    temp.char$NOSC[i] = mean(temp$NOSC)
    temp.char$DBE[i] = mean(temp$DBE)
    temp.char$AI_Mod[i] = mean(temp$AI_Mod)
    temp.char$NtoC[i] = mean(temp$NtoC_ratio)

    rm(temp)
  }
  
  # Adding information to overall dataset
  temp.char = melt(as.matrix(temp.char))
  temp.char$Type = names(data)[a]
  char = rbind(char, temp.char)
  
  rm(temp.char, temp.data, temp.mol)
}

rm(a)

# Plotting results
print(
  ggplot(char, aes(x = Type, y = value))+
    geom_boxplot(aes(color = Type)) + geom_jitter(aes(color = Type))+
    facet_wrap(.~Var2, ncol = 2, scales = "free_y")+
    scale_color_calc() + xlab(NULL) + theme_bw() + 
    theme(legend.position = "none", text = element_text(size = 14),
          axis.text.x = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black"),
          panel.border = element_rect(size = 1, colour = "black"),
          panel.background = element_blank(),
          panel.grid = element_blank())
)

# Calculating stats
char.stats = NULL # Empty object for stat results
type = unique(char$Type) # Methods
char.list = unique(char$Var2) # Elemental compositions

for(curr.char in char.list){
  temp = char[which(char$Var2 %in% curr.char),]
  
  for(i in 1:(length(type)-1)){
    for(j in (i+1):length(type)){
      temp.type = temp[which(temp$Type %in% c(type[i], type[j])),]
      
      mwu = wilcox.test(value~Type, data = temp.type, paired = T)
      temp.stat = data.frame(Comparison = paste0(curr.char, " ", type[i], "-", type[j]),
                             MWU = mwu$statistic, p.value = mwu$p.value)
      char.stats = rbind(char.stats, temp.stat)
    }
  }
  rm(temp, temp.type, temp.stat, mwu)
}

rm(curr.char, i, j)

# p-value adjustment
char.stats$p.value = p.adjust(char.stats$p.value, method = "fdr")


### Elemental composition analyses
el.comp = NULL

# Loop through samples to get relative abundance
for(a in 1:length(data)){
  # Select data
  temp.data = data[[a]]
  temp.mol = mol[[a]]
  
  # Temporary elemental composition object
  temp.comp = matrix(data = 0, nrow = ncol(temp.data), ncol = length(unique(temp.mol$El_comp)),
                     dimnames = list(colnames(temp.data), unique(temp.mol$El_comp)))
  
  for(i in 1:nrow(temp.comp)){
    temp = temp.mol[which(temp.data[,i] > 0),] # Mol data for a given sample
    
    for(j in 1:ncol(temp.comp)){
      temp.comp[i,j] = length(which(temp$El_comp %in% colnames(temp.comp)[j]))
    }
  }
  
  # Converting to relative abundance
  temp.comp = t(apply(temp.comp, 1, function(x) (x/sum(x))*100))
  
  # Adding information to overall dataset
  temp.comp = melt(temp.comp)
  temp.comp$Type = names(data)[a]
  el.comp = rbind(el.comp, temp.comp)
  
  rm(temp.comp, temp.data, temp.mol, temp)
}

rm(a)

# Plot the results
print(
  ggplot(el.comp, aes(x = Type, y = value))+
    geom_boxplot(aes(color = Type)) + geom_jitter(aes(color = Type))+
    facet_wrap(.~Var2, ncol = 2, scales = "free_y")+
    scale_color_calc() + xlab(NULL) + ylab("Rel. Abund. (%)") + theme_bw() + 
    theme(legend.position = "none", text = element_text(size = 14),
          axis.text.x = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black"),
          panel.border = element_rect(size = 1, colour = "black"),
          panel.background = element_blank(),
          panel.grid = element_blank())
)

# Calculating stats
el.stats = NULL # Empty object for stat results
type = unique(el.comp$Type) # Methods
elem = unique(el.comp$Var2) # Elemental compositions

for(curr.elem in elem){
  temp = el.comp[which(el.comp$Var2 %in% curr.elem),]
  
  for(i in 1:(length(type)-1)){
    for(j in (i+1):length(type)){
      temp.type = temp[which(temp$Type %in% c(type[i], type[j])),]
      
      mwu = wilcox.test(value~Type, data = temp.type, paired = pair)
      temp.stat = data.frame(Comparison = paste0(curr.elem, " ", type[i], "-", type[j]),
                             MWU = mwu$statistic, p.value = mwu$p.value)
      el.stats = rbind(el.stats, temp.stat)
    }
  }
  rm(temp, temp.type, temp.stat, mwu)
}

rm(curr.elem, i, j)

# p-value adjustment
el.stats$p.value = p.adjust(el.stats$p.value, method = "fdr")


### Compoud class analyses
comp.class = NULL

# Loop through samples to get relative abundance
for(a in 1:length(data)){
  # Select data
  temp.data = data[[a]]
  temp.mol = mol[[a]]
  
  # Temporary elemental composition object
  temp.comp = matrix(data = 0, nrow = ncol(temp.data), ncol = length(unique(temp.mol$Class)),
                     dimnames = list(colnames(temp.data), unique(temp.mol$Class)))
  
  for(i in 1:nrow(temp.comp)){
    temp = temp.mol[which(temp.data[,i] > 0),] # Mol data for a given sample
    
    for(j in 1:ncol(temp.comp)){
      temp.comp[i,j] = length(which(temp$Class %in% colnames(temp.comp)[j]))
    }
  }
  
  rm(i,j)
  
  # Converting to relative abundance
  temp.comp = t(apply(temp.comp, 1, function(x) (x/sum(x))*100))
  
  # Adding information to overall dataset
  temp.comp = melt(temp.comp)
  temp.comp$Type = names(data)[a]
  comp.class = rbind(comp.class, temp.comp)
  
  rm(temp.comp, temp.data, temp.mol, temp)
}

rm(a)

# Plot the results
print(
  ggplot(comp.class, aes(x = Type, y = value))+
    geom_boxplot(aes(color = Type)) + geom_jitter(aes(color = Type))+
    facet_wrap(.~Var2, ncol = 2, scales = "free_y")+
    scale_color_calc() + xlab(NULL) + ylab("Rel. Abund. (%)") + theme_bw() + 
    theme(legend.position = "none", text = element_text(size = 14),
          axis.text.x = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black"),
          panel.border = element_rect(size = 1, colour = "black"),
          panel.background = element_blank(),
          panel.grid = element_blank())
)

# Calculating stats
comp.stats = NULL # Empty object for stat results
type = unique(comp.class$Type) # Methods
comp = unique(comp.class$Var2) # Elemental compositions

for(curr.comp in comp){
  temp = comp.class[which(comp.class$Var2 %in% curr.comp),]
  
  for(i in 1:(length(type)-1)){
    for(j in (i+1):length(type)){
      temp.type = temp[which(temp$Type %in% c(type[i], type[j])),]
      
      mwu = wilcox.test(value~Type, data = temp.type, paired = pair)
      temp.stat = data.frame(Comparison = paste0(curr.comp, " ", type[i], "-", type[j]),
                             MWU = mwu$statistic, p.value = mwu$p.value)
      comp.stats = rbind(comp.stats, temp.stat)
    }
  }
  rm(temp, temp.type, temp.stat, mwu)
}

rm(curr.comp, i, j)

# p-value adjustment
comp.stats$p.value = p.adjust(comp.stats$p.value, method = "fdr")


### Permanova
# Distance
dist = vegdist(t(mol_data), method = "jaccard")

# Make factors
factors = data.frame(SPE_Mode = str_extract(row.names(as.matrix(dist)), "^SPE|^NonSPE"),
                     ESI_Mode = str_extract(row.names(as.matrix(dist)), "Pos|Neg"))
factors$SPE_ESI = paste0(factors$SPE_Mode, "_", factors$ESI_Mode)

# Overall stats
method.perm = adonis(dist~factors$SPE_ESI)

# Pairwise stats
pair.perm = NULL
uniq.type = unique(factors$SPE_ESI)

for(i in 1:(length(uniq.type)-1)){
  for(j in (i+1):length(uniq.type)){
    w = which(factors$SPE_ESI %in% c(uniq.type[i], uniq.type[j]))
    temp = as.matrix(dist)[w,w]
    temp.perm = adonis(as.dist(temp)~factors$SPE_ESI[w])
    
    pair.perm = rbind(pair.perm,
                      data.frame(Comparison = paste(uniq.type[i], "-", uniq.type[j]),
                                 Pseudo_F = temp.perm$aov.tab$F.Model[1],
                                 p_value = temp.perm$aov.tab$`Pr(>F)`[1]))
  }
}

# Plotting NMDS
nms = metaMDS(dist)
nms = as.data.frame(scores(nms))
nms$SPE_ESI = paste0(factors$SPE_Mode, "_", factors$ESI_Mode)

print(
  ggplot(nms, aes(x = NMDS1, y = NMDS2))+
    geom_point(aes(color = SPE_ESI), size = 4)+
    scale_color_calc() + theme_bw() + 
    theme(text = element_text(size = 14),
          axis.text.x = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black"),
          panel.border = element_rect(size = 1, colour = "black"),
          panel.background = element_blank(),
          panel.grid = element_blank())
)
