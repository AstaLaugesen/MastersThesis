#Speciale plots
#####relative abundances####
rm(list=ls())

#loading libraries
library(tidyverse)
library(phyloseq)
library(readxl)
library(ggplot2)


#reading in files
setwd("~/data")
otu_mat<- t(read.table("defensesystems/defsystems_abundance.txt", header=TRUE, sep="\t", row.names = 1))
tax_mat<- read_excel("defensesystems/deftaxonomy.xlsx", sheet = "Data")
sampleIDs <-read.table("clinical/samples.abc.tab", header=FALSE, sep="\t")

#defining row names in the taxonomy table
tax_mat <- data.frame(tax_mat, row.names = 1)
#removing the patient no. that wasn't sequenced (no sample number)
sampleIDs <- sampleIDs[-c(403),]
#adding the column with the S[...] identifiers and adding column titles
#(Done because samples can't start with a number)
sampleIDs$V0 <- paste("S", sampleIDs$V1, sep="")
colnames(sampleIDs) <- c("SampleName", "ABCno", "NewNames")
sample_df <- sampleIDs

#defining rows in sample tables
row.names(sample_df) <- sample_df$NewNames
#removing the new names since they're row names now
sample_df <- sample_df %>% select(-NewNames)


#And now putting it all together:

#making OTU and tax table into matrix
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)

#Transform into phyloseq object
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
SAMPLES = sample_data(sample_df)

DefSys <- phyloseq(OTU,TAX, SAMPLES)

library(devtools)
#install_github("JStokholm/rabuplot")
library(rabuplot)

#Overview:
rabuplot(DefSys, predictor="ABCno", type="Type", main="Relative Abundances of Defense Systems in All Samples", bar_chart_stacked = TRUE) + 
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
  xlab("All 662 Samples") +
  ylab("Relative Abundance")+
  labs(title = "Relative Abundances of Defense Systems in Samples", tag = "A")
rabuplot(DefSys, type="Type", bar_chart_stacked = TRUE) + 
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
  xlab("") +
  ylab("Average Relative Abundance")+
  labs(title = "Average Relative Abundances of Defense Systems", tag = "B")
ggsave("MASTplots/overview.pdf", width = 20, height =10 )


#####count/percentages plot####
#loading libraries
library(tidyverse)
library(plyr)
library(dplyr)
library(ggplot2)

#load in data
setwd("~/data")
table_clusters<- read.table("defensesystems/ALL_defsys_contigIn200kCluster_unique.txt", header=FALSE, sep="\t")
table <- table_clusters
#set the actual header as a header
colnames(table) <- table_clusters[c(145433),]
#remove that row now that its a header
table <- table[-c(145433),]
#make the 0's and 1's into actual numbers, not characters
table <- transform(table, contig_in_200kclust = as.numeric(contig_in_200kclust))


table %>% group_by(contig_in_200kclust) %>% dplyr::summarise(N=n())
#table %>% group_by(type) %>% summarize(N=n()))

#summarize and see how many 1's there are for each type
summaryTab <- table %>%
  group_by(type) %>%
  dplyr::summarize(sum(contig_in_200kclust))
#find the total amount of that type
rows_type <- ddply(table, c("type"), summarise, nrows = length(type))
#join the tables
summaryTab <- summaryTab %>% inner_join(rows_type,by="type")
colnames(summaryTab) <- c("type", "in_cluster", "total_n_rows")

#calculate percentages
summaryTab <- summaryTab %>% 
  mutate(percentage_in = in_cluster/total_n_rows,
         percentage_out = 1-percentage_in)
#pivot to long format
long_summaryTab<- summaryTab %>%
  select(-in_cluster, -total_n_rows) %>% 
  pivot_longer(!type, names_to = "condition", values_to = "percentage")
long_summaryTab <- long_summaryTab %>% inner_join(summaryTab,by="type") %>% 
  select(-in_cluster, -percentage_in, -percentage_out)


# plot --------------------------------------------------------------------

#install.packages("ggthemes")
library(ggplot2)
library(ggthemes)

brks <- c(0, 0.25, 0.5, 0.75, 1)

#RM is most abundant. Better ways to show this
ggplot(data=summaryTab,aes(x=reorder(type, -in_cluster),y=in_cluster)) +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust = 1))

#make plot:
#unsorted:
#long_summaryTab<- long_summaryTab %>% fct_reorder(condition, n, .desc = TRUE)
long_summaryTab$condition <- factor(long_summaryTab$condition, levels = c("percentage_out", "percentage_in"))
ggplot(data=long_summaryTab,aes(x=type, y=percentage,fill=condition)) +
  geom_bar(stat="identity") +
  scale_y_continuous(breaks = brks, labels = scales::percent(brks))+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust = 1))+
  ggtitle("Percentage of defense systems in clusters over 200k across all children")+
  scale_fill_manual(name="In cluster",labels=c("Not In","In"), values=c("#D3494E","#23C552"))+
  ylab("Percentage") +
  xlab("Types of Defense Systems")+
  geom_text(aes(y=0.07,label=paste0("n=",total_n_rows)),angle = 90, size=2.5)
#ggsave("MASTplots/DefenseSystemsInClustOver200kbp.png", width = 11, height =6)

#sorted:
ggplot(data=long_summaryTab,aes(x=reorder(type, -total_n_rows), y=percentage,fill=condition)) +
  geom_bar(stat="identity") +
  scale_y_continuous(breaks = brks, labels = scales::percent(brks))+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust = 1))+
  ggtitle("Percentage of defense systems in clusters over 200k across all children")+
  scale_fill_manual(name="In cluster",labels=c("Not In","In"), values=c("#D3494E","#23C552"))+
  ylab("Percentage") +
  xlab("Types of Defense Systems")+
  geom_text(aes(y=0.07,label=paste0("n=",total_n_rows)),angle = 90, size=2.5)
ggsave("MASTplots/DefenseSystemsInClustOver200kbpSIZESORT.png", width = 11.5, height =6)




