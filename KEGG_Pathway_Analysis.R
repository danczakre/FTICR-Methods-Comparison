### Processing East River pathway mapping data

library(dplyr)
library(reshape2)
library(ggplot2)
library(ggpubr)

setwd("~/Documents/East River/New Comparative Analyses/")

# Load in pathway data
path = cbind(Method = "SPE (+)", read.table("SPE_Pos/SPE_Pos_Pathways_from_Compounds.txt", sep = "\t", header = 1))
path = rbind(path,
             cbind(Method = "SPE (-)", read.table("SPE_Neg/SPE_Neg_Pathways_from_Compounds.txt", sep = "\t", header = 1))
)
path = rbind(path,
             cbind(Method = "NonSPE (+)", read.table("NonSPE_Pos/NonSPE_Pos_Pathways_from_Compounds.txt", sep = "\t", header = 1))
)
path = rbind(path,
              cbind(Method = "NonSPE (-)", read.table("NonSPE_Neg/NonSPE_Neg_Pathways_from_Compounds.txt", sep = "\t", header = 1))
)

# Removing CPDs without pathway information
path = path[-which(path$Pathway_Name %in% "Unavailable"),]

# Add in count data for casting data
path$Count = 1

# Casting and rearranging data
path = dcast(path, Method~Pathway_Name, value.var = "Count", fun.aggregate = sum)
row.names(path) = path$Method; path = path[,-1]
rel.path = melt(apply(path, 1, function(x) (x/sum(x))*100))
path = melt(as.matrix(path))

# Plotting data
plot1 = path %>% group_by(Var1) %>% top_n(10, value) %>%
  ggplot(aes(x = reorder(Var2, -value), y = value))+
  geom_bar(stat = "identity")+
  facet_wrap(.~Var1, ncol = 1)+
  xlab(NULL) + ylab("Compound Count")+
  theme_bw() + theme(axis.text.x = element_text(color = "black", angle = 90, hjust = 1, vjust = 0.5),
                     axis.text = element_text(color = "black", size = 12),
                     axis.title = element_text(color = "black", size = 14),
                     axis.ticks = element_line(color = "black"),
                     panel.background = element_blank(),
                     panel.grid = element_blank(),
                     panel.border = element_rect(color = "black"),
                     strip.text = element_text(color = "black", size = 12))

plot2 = rel.path %>% group_by(Var2) %>% top_n(10, value) %>%
  ggplot(aes(x = reorder(Var1, -value), y = value))+
  geom_bar(stat = "identity")+
  facet_wrap(.~Var2, ncol = 1)+
  xlab(NULL) + ylab("Relative Abundance (%)")+
  theme_bw() + theme(axis.text.x = element_text(color = "black", angle = 90, hjust = 1, vjust = 0.5),
                     axis.text = element_text(color = "black", size = 12),
                     axis.title = element_text(color = "black", size = 14),
                     axis.ticks = element_line(color = "black"),
                     panel.background = element_blank(),
                     panel.grid = element_blank(),
                     panel.border = element_rect(color = "black"),
                     strip.text = element_text(color = "black", size = 12))

ggarrange(plot1, plot2, labels = c("a)", "b)"), 
          font.label = list(size = 18, face = "plain"))
