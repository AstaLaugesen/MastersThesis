rm(list=ls())

#loading libraries
library(tidyverse)
library(ggplot2)

#load in data
setwd("~/data")
Abi2<- read.table("defensesystems/fractions/BactDefsysFraction_Abi2.txt", header=TRUE, sep="\t", row.names = "Sample")
AbiEii<-read.table("defensesystems/fractions/BactDefsysFraction_AbiEii.txt", header=TRUE, sep="\t", row.names = "Sample")
AbiH<-read.table("defensesystems/fractions/BactDefsysFraction_AbiH.txt", header=TRUE, sep="\t", row.names = "Sample")
AVAST<-read.table("defensesystems/fractions/BactDefsysFraction_AVAST.txt", header=TRUE, sep="\t", row.names = "Sample")
BREX<-read.table("defensesystems/fractions/BactDefsysFraction_BREX.txt", header=TRUE, sep="\t", row.names = "Sample")
Cas<-read.table("defensesystems/fractions/BactDefsysFraction_Cas.txt", header=TRUE, sep="\t", row.names = "Sample")
CasClass1<-read.table("defensesystems/fractions/BactDefsysFraction_CasClass1.txt", header=TRUE, sep="\t", row.names = "Sample")
CasClass2<-read.table("defensesystems/fractions/BactDefsysFraction_CasClass2.txt", header=TRUE, sep="\t", row.names = "Sample")
CBASS<-read.table("defensesystems/fractions/BactDefsysFraction_CBASS.txt", header=TRUE, sep="\t", row.names = "Sample")
Gabija<-read.table("defensesystems/fractions/BactDefsysFraction_Gabija.txt", header=TRUE, sep="\t", row.names = "Sample")
Hachiman<-read.table("defensesystems/fractions/BactDefsysFraction_Hachiman.txt", header=TRUE, sep="\t", row.names = "Sample")
Kiwa<-read.table("defensesystems/fractions/BactDefsysFraction_Kiwa.txt", header=TRUE, sep="\t", row.names = "Sample")
LamassuFam<-read.table("defensesystems/fractions/BactDefsysFraction_LamassuFam.txt", header=TRUE, sep="\t", row.names = "Sample")
Mokosh<-read.table("defensesystems/fractions/BactDefsysFraction_Mokosh.txt", header=TRUE, sep="\t", row.names = "Sample")
Nhi<-read.table("defensesystems/fractions/BactDefsysFraction_Nhi.txt", header=TRUE, sep="\t", row.names = "Sample")
Retron<-read.table("defensesystems/fractions/BactDefsysFraction_Retron.txt", header=TRUE, sep="\t", row.names = "Sample")
RM<-read.table("defensesystems/fractions/BactDefsysFraction_RM.txt", header=TRUE, sep="\t", row.names = "Sample")
RstPARIS<-read.table("defensesystems/fractions/BactDefsysFraction_RstPARIS.txt", header=TRUE, sep="\t", row.names = "Sample")
SEFIR<-read.table("defensesystems/fractions/BactDefsysFraction_SEFIR.txt", header=TRUE, sep="\t", row.names = "Sample")
Septu<-read.table("defensesystems/fractions/BactDefsysFraction_Septu.txt", header=TRUE, sep="\t", row.names = "Sample")
Shedu<-read.table("defensesystems/fractions/BactDefsysFraction_Shedu.txt", header=TRUE, sep="\t", row.names = "Sample")
ShosTA<-read.table("defensesystems/fractions/BactDefsysFraction_ShosTA.txt", header=TRUE, sep="\t", row.names = "Sample")
SoFIC<-read.table("defensesystems/fractions/BactDefsysFraction_SoFIC.txt", header=TRUE, sep="\t", row.names = "Sample")
Thoeris<-read.table("defensesystems/fractions/BactDefsysFraction_Thoeris.txt", header=TRUE, sep="\t", row.names = "Sample")
Wadjet<-read.table("defensesystems/fractions/BactDefsysFraction_Wadjet.txt", header=TRUE, sep="\t", row.names = "Sample")

#test if binary-like distribution is seen
RM_test <- RM %>% rownames_to_column() %>% pivot_longer(-rowname, names_to = "bacteria", values_to = "value")
ggplot(RM_test, aes(value)) + 
  geom_histogram(binwidth = 0.1)+ 
  xlim(c(-0.05,2))+
  ggtitle("Occurence of ratio data",
          subtitle="How often do species of bacterial genus have a defense system")


#make it into 0 and 1
Abi2_abs <- Abi2 %>% 
  mutate_all(
    funs(case_when(
      . < 0.15  ~ 0,
      . >= 0.15 ~ 1))) 

AbiH_abs <- AbiH %>% 
  mutate_all(
    funs(case_when(
      . < 0.15  ~ 0,
      . >= 0.15 ~ 1)))

AbiEii_abs <- AbiEii %>% 
  mutate_all(
    funs(case_when(
      . < 0.15  ~ 0,
      . >= 0.15 ~ 1)))
AVAST_abs <- AVAST %>% 
  mutate_all(
    funs(case_when(
      . < 0.15  ~ 0,
      . >= 0.15 ~ 1)))
BREX_abs  <- BREX %>% 
  mutate_all(
    funs(case_when(
      . < 0.15  ~ 0,
      . >= 0.15 ~ 1)))
Cas_abs  <- Cas %>% 
  mutate_all(
    funs(case_when(
      . < 0.15  ~ 0,
      . >= 0.15 ~ 1)))
CasClass1_abs  <- CasClass1 %>% 
  mutate_all(
    funs(case_when(
      . < 0.15  ~ 0,
      . >= 0.15 ~ 1)))
CasClass2_abs  <- CasClass2 %>% 
  mutate_all(
    funs(case_when(
      . < 0.15  ~ 0,
      . >= 0.15 ~ 1)))
CBASS_abs <- CBASS %>% 
  mutate_all(
    funs(case_when(
      . < 0.15  ~ 0,
      . >= 0.15 ~ 1)))
Gabija_abs <- Gabija %>% 
  mutate_all(
    funs(case_when(
      . < 0.15  ~ 0,
      . >= 0.15 ~ 1)))
Hachiman_abs <- Hachiman %>% 
  mutate_all(
    funs(case_when(
      . < 0.15  ~ 0,
      . >= 0.15 ~ 1)))
Kiwa_abs  <- Kiwa %>% 
  mutate_all(
    funs(case_when(
      . < 0.15  ~ 0,
      . >= 0.15 ~ 1)))
LamassuFam_abs <- LamassuFam %>% 
  mutate_all(
    funs(case_when(
      . < 0.15  ~ 0,
      . >= 0.15 ~ 1)))
Mokosh_abs <- Mokosh %>% 
  mutate_all(
    funs(case_when(
      . < 0.15  ~ 0,
      . >= 0.15 ~ 1)))
Nhi_abs <- Nhi %>% 
  mutate_all(
    funs(case_when(
      . < 0.15  ~ 0,
      . >= 0.15 ~ 1)))
Retron_abs<- Retron %>% 
  mutate_all(
    funs(case_when(
      . < 0.15  ~ 0,
      . >= 0.15 ~ 1)))
RM_abs <- RM %>% 
  mutate_all(
    funs(case_when(
      . < 0.15  ~ 0,
      . >= 0.15 ~ 1)))
RstPARIS_abs <- RstPARIS %>% 
  mutate_all(
    funs(case_when(
      . < 0.15  ~ 0,
      . >= 0.15 ~ 1)))
SEFIR_abs <- SEFIR %>% 
  mutate_all(
    funs(case_when(
      . < 0.15  ~ 0,
      . >= 0.15 ~ 1)))
Septu_abs <- Septu %>% 
  mutate_all(
    funs(case_when(
      . < 0.15  ~ 0,
      . >= 0.15 ~ 1)))
Shedu_abs <- Shedu %>% 
  mutate_all(
    funs(case_when(
      . < 0.15  ~ 0,
      . >= 0.15 ~ 1)))
ShosTA_abs <- ShosTA %>% 
  mutate_all(
    funs(case_when(
      . < 0.15  ~ 0,
      . >= 0.15 ~ 1)))
SoFIC_abs <- SoFIC %>% 
  mutate_all(
    funs(case_when(
      . < 0.15  ~ 0,
      . >= 0.15 ~ 1)))
Thoeris_abs <- Thoeris %>% 
  mutate_all(
    funs(case_when(
      . < 0.15  ~ 0,
      . >= 0.15 ~ 1)))
Wadjet_abs <- Wadjet %>% 
  mutate_all(
    funs(case_when(
      . < 0.15  ~ 0,
      . >= 0.15 ~ 1)))

#make into a tensor, and indicate row names etc
all_tables <- array(list(Abi2_abs, AbiEii_abs, AbiH_abs, AVAST_abs, BREX_abs, Cas_abs, 
                         CasClass1_abs,CasClass2_abs, CBASS_abs, Gabija_abs, Hachiman_abs, 
                         Kiwa_abs, LamassuFam_abs, Mokosh_abs, Nhi_abs, Retron_abs, RM_abs, 
                         RstPARIS_abs, SEFIR_abs, Septu_abs, Shedu_abs, ShosTA_abs, 
                         SoFIC_abs, Thoeris_abs, Wadjet_abs) %>% unlist, 
                    dim=c(662, 50, 25), 
                    dimnames = list(sample = rownames(Abi2),  bacteria = colnames(Abi2), 
                                    defsys = c("Abi2", "AbiEii", "AbiH", "AVAST", "BREX", "Cas", 
                                               "CasClass1", "CasClass2", "CBASS", "Gabija", 
                                               "Hachiman", "Kiwa", "LamassuFam", "Mokosh", "Nhi",
                                               "Retron", "RM","RstPARIS", "SEFIR", "Septu", 
                                               "Shedu", "ShosTA", "SoFIC", "Thoeris", "Wadjet")))
#we can view it with 
all_tables[,,1:2]

##### Parafac ####
library(multiway)
library(tidyverse)
#install.packages("devtools")
#devtools::install_github("thomasp85/patchwork")
library(patchwork)

#no centering or scaling, not needed here when its binary. Just run the models
pfac <- parafac(all_tables, nfac = 3, nstart = 100)
pfac2 <- parafac(all_tables, nfac = 2, nstart = 100)
#run with constraints
pfac2const <- parafac(all_tables, nfac = 2, nstart = 100, const = c("nonneg","nonneg", "uncons"))

# Extract three modes
#if the original plot, run this:
children <- pfac$A %>% data.frame %>% bind_cols(children = dimnames(all_tables)[[1]])
bact <- pfac$B %>% data.frame %>% bind_cols(bact = dimnames(all_tables)[[2]])
DefSys <- pfac$C %>% data.frame %>% bind_cols(DefSys = dimnames(all_tables)[[3]])
#if the optimized plot, run this:
#constrained
children <- pfac2const$A %>% data.frame %>% bind_cols(children = dimnames(all_tables)[[1]])
bact <- pfac2const$B %>% data.frame %>% bind_cols(bact = dimnames(all_tables)[[2]])
DefSys <- pfac2const$C %>% data.frame %>% bind_cols(DefSys = dimnames(all_tables)[[3]])



bactplot <- bact %>% pivot_longer(-bact, names_to = "comp", values_to = "value") %>%
  filter(abs(value) > 0.75) %>% 
  ggplot(aes(bact, value, fill = bact)) +
  geom_col() +
  facet_wrap(~ comp) +
  ylab("Loading")+
  xlab("Bacteria")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.4))+
  theme(legend.position="none")
bactplot

defsysplot <- DefSys %>% pivot_longer(-DefSys, names_to = "comp", values_to = "value") %>%
  #filter(abs(value) > 1.5) %>% 
  ggplot(aes(DefSys, value, fill=DefSys)) +
  geom_col() +
  facet_wrap(~ comp) +
  xlab("Defense Systems")+
  ylab("Loading")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.4))+
  theme(legend.position="none")
defsysplot

#sort the data of the children (neighbor joining):
tree <- ape::nj(children %>% column_to_rownames(var = "children") %>% dist)
#plot(tree, "u")
library(ape)
write.tree(tree, file = "tree.nwk")
#save and do following unix command:
#cat tree.nwk | sed 's/1809/\n1809/g' | sed 's/:.\+//' > tree.order.txt

namesorder <- scan("tree.order.txt", what = character())

childplotsort <- children %>%
  pivot_longer(-children, names_to = "comp", values_to = "value") %>%
  mutate(children = factor(children, levels = namesorder)) %>% 
  ggplot(aes(children, value)) +
  geom_col()+
  facet_wrap(~ comp) +
  ylab("Loading")+
  xlab("Children")+
  theme(axis.text.x = element_blank()) +
  geom_vline(xintercept=seq(1.5, 662-0.5, 60), linewidth=0.2, colour="blue")

#for displaying optimized one
childplotsort /
  bactplot/
  defsysplot +
  plot_annotation(title = 'Parafac Model Over Defense Systems in Top 50 Most Abundant Bacteria', 
                  subtitle = 'Bacteria with all absolute loadings of <0.75 omitted')
# for displaying big one:
childplotsort /
  bactplot/
  defsysplot +
  plot_annotation(title = 'Parafac Model Over Defense Systems in Top 50 Most Abundant Bacteria', 
                  subtitle = 'All 50 bacteria included, non-negative contraints')





##### run the tests on what is the optimal model ####

#set number of components we want to test
maxcomp <- 10
#run parafac for those
pfac_Rsq <- sapply(1:maxcomp, function(i){
  parafac(all_tables, nfac = i, nstart = 25)$Rsq
})

## PCA alternative
liste <- list()
for(i in 1:dim(all_tables)[3]){
  liste[[i]] <- all_tables[,,i]
  names(liste)[i] <- dimnames(all_tables)[[3]][i]
}
wide_alltables <- bind_rows(liste %>% lapply(data.frame) %>% lapply(rownames_to_column), .id = "id") %>% 
  pivot_longer(Bacteroides:Faecalibacillus, names_to = "bacteria", values_to = "value") %>% 
  unite("name", c(id, bacteria)) %>% 
  pivot_wider()
defsys_pca <- prcomp(wide_alltables %>% column_to_rownames, center=FALSE)

#calculate R^2 to compare
varexpdata <- data.frame(comp = 1:length(defsys_pca$sdev), sd = defsys_pca$sdev) %>% 
  mutate(var = sd^2,
         cumvar = cumsum(var),
         cumpct_pca = cumvar/max(cumvar)) %>% 
  slice(1:maxcomp) %>% 
  bind_cols(parafac_Rsq = pfac_Rsq)

#make plot
varexpdata %>% 
  pivot_longer(cumpct_pca:parafac_Rsq) %>% 
  mutate(value=round(value, 2))%>% 
  ggplot(aes(comp, value, color = name, label=value)) +
  geom_point(size = 2) +
  geom_line() +
  ylim(c(0.04,0.96)) + 
  #for not clogging the plot, only display parafac values
  geom_text(nudge_y = -0.05, size=3, aes(label=ifelse(name=="parafac_Rsq", value,"")))+
  ylab("Variance Explained")+
  xlab("Number of Components")+
  ggtitle("Variance Explained (R^2) - Parafac VS PCA",
          subtitle = "For Top 50 most abundant bacteria")+
  scale_colour_discrete(name="Method",
                        labels=c("PCA", "Parafac"))+
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9,10))






#####Rabuplot####
#load in virome data
library(phyloseq)
load("~/data/vOTUs.RData", verbose= TRUE)
load("~/data/VAMB/220729_VAMB_final_phyloseq_with_tree.RData", verbose = TRUE)

#agglomerate VAMB data to genus level
final.physeq
VAMB <- tax_glom(final.physeq, taxrank="Genus")

#vOTU data
vOTUs

#check what genera are available to compare
tax_table(vOTUs) %>% 
  as.data.frame %>% 
  select(hostGenus) %>% 
  distinct(., hostGenus) %>% 
  arrange(., hostGenus)



###BIFIDOBACTERIUM###
#subset viruses that infect Bifidobacterium
vOTUs.Bifidobacterium_unglom <- subset_taxa(vOTUs, hostGenus == "Bifidobacterium")
vOTUs.Bifidobacterium <- tax_glom(vOTUs.Bifidobacterium_unglom, taxrank="family")


#cluster 125 is Bifidobacterium, lets take it out
clustname <- tax_table(VAMB) %>% 
  as.data.frame() %>% 
  filter(str_detect(Genus, 'Bifidobacterium')) %>% 
  rownames()

#samples in VAMB: 662
Bifidobacterium.VAMBabundance <- otu_table(VAMB) %>% 
  as.data.frame %>% 
  select(paste(clustname)) %>% 
  rownames_to_column() %>% 
  arrange(rowname) %>% 
  #change rownames so they match virome data
  mutate(rowname = gsub("[.]", "-", rowname)) %>% 
  mutate(rowname = gsub("X", "S", rowname))

#Samples in virome: 648
tableotu <- otu_table(vOTUs.Bifidobacterium) %>% 
  as.data.frame %>% 
  select(sort(names(.)))

#choose only samples that are in both tables:
Bifidobacterium.VAMBabundance <- Bifidobacterium.VAMBabundance[Bifidobacterium.VAMBabundance$rowname %in%
                                                                 colnames(otu_table(vOTUs.Bifidobacterium)),]
Bifidobacterium.VAMBvector <- Bifidobacterium.VAMBabundance[[paste(clustname)]]

#and then divide otu table of vOTUs with the vector
otu_table(vOTUs.Bifidobacterium) <- otu_table(sweep(tableotu,2,Bifidobacterium.VAMBvector,'/'), taxa_are_rows=TRUE)
##ok done dividing, lets go on

#add row names to table
Bifidobacterium_table <- rownames_to_column(all_tables[,,1] %>% as.data.frame) %>% 
  select(rowname)
#add each column from each defense system to table
for(i in 1:length(defsyslist)){
  column <- all_tables[,"Bifidobacterium",defsyslist[i]] %>% as.data.frame()
  colnames(column) <- c(defsyslist[i])
  #print(column)
  Bifidobacterium_table[defsyslist[i]] <- column
}

sample_data(vOTUs.Bifidobacterium) <- left_join(get_variable(vOTUs.Bifidobacterium), 
                                                Bifidobacterium_table, by = c("SampleName" = "rowname")) %>% 
  column_to_rownames(., var = "NewNames")


library(rabuplot)

Bifidobacterium.abi2plot <- rabuplot(vOTUs.Bifidobacterium, "Abi2", type = "family", p_adjust=TRUE, p_adjust_method = "fdr")
Bifidobacterium.CasClass1plot <- rabuplot(vOTUs.Bifidobacterium, "CasClass1", type = "family", p_adjust=TRUE, p_adjust_method = "fdr")
Bifidobacterium.RMplot <- rabuplot(vOTUs.Bifidobacterium, "RM", type = "family", p_adjust=TRUE, p_adjust_method = "fdr")
Bifidobacterium.wadjetplot <- rabuplot(vOTUs.Bifidobacterium, "Wadjet", type = "family", p_adjust=TRUE, p_adjust_method = "fdr")


(Bifidobacterium.abi2plot+theme(axis.text.y = element_text(size = 13)) |
    Bifidobacterium.CasClass1plot+theme(axis.text.y = element_text(size = 13)))  /
  (Bifidobacterium.RMplot+theme(axis.text.y = element_text(size = 13)) | 
     Bifidobacterium.wadjetplot+theme(axis.text.y = element_text(size = 13)))+
  plot_annotation(title = 'Differential Analysis for Abundances of Virus Targeting Bifidobacterium - Defense Systems Absence/Presence (no log transform)',
                  subtitle = "Abundances adjusted for bacteria abundance, p-values adjusted with FDR method", 
                  theme = theme(plot.title = element_text(size = 14.3)))
#ggsave("MASTplots/fractions_bifido_NOTtransformedAndAgglomerated.jpeg", width=13.5, height=7)



#####CLR TRANSFORM####


#extract family names and OTU names, for replacing OTU names with real family names
glom_vir_tax <- tax_table(vOTUs) %>% data.frame()
glom_vir_tax <- rownames_to_column(glom_vir_tax)
OTUnameID <- data.frame(glom_vir_tax$family, glom_vir_tax$rowname) %>%
  rename(family=glom_vir_tax.family) %>% 
  rename(OTUname=glom_vir_tax.rowname)

#####Bifidobacterium####

#subset and agglomerate
vOTUs.Bifidobacterium_unglom <- subset_taxa(vOTUs, hostGenus == "Bifidobacterium")
vOTUs.Bifidobacterium <- tax_glom(vOTUs.Bifidobacterium_unglom, taxrank="family")

#otu names to column
vOTUs.Bifidobacterium2 <- otu_table(vOTUs.Bifidobacterium) %>%
  data.frame() %>% 
  rownames_to_column()
#join ID and OTU names
Bifidobacterium_newnames <- left_join(vOTUs.Bifidobacterium2, OTUnameID, by = c("rowname" = "OTUname"))
Bifidobacterium_newnames <- Bifidobacterium_newnames %>% 
  column_to_rownames(var="family") %>% select(-rowname)

#CLR transform
clr_Bifidobacterium <- Bifidobacterium_newnames %>% 
  #find minimum where we have replaced 0'es with NA to not just find that 0 is minimum
  + (min(replace(Bifidobacterium_newnames, Bifidobacterium_newnames<=0, NA), na.rm = TRUE)/2) %>% 
  clr() %>% 
  t() %>% 
  data.frame() %>% 
  rownames_to_column()

#add row names to table
Bifidobacterium_table <- rownames_to_column(all_tables[,,1] %>% as.data.frame) %>% 
  select(rowname)
#add each column from each defense system to table
for(i in 1:length(defsyslist)){
  column <- all_tables[,"Bifidobacterium",defsyslist[i]] %>% as.data.frame()
  colnames(column) <- c(defsyslist[i])
  #print(column)
  Bifidobacterium_table[defsyslist[i]] <- column
}

Bifidobacterium_table <- Bifidobacterium_table %>%  
  mutate(rowname = gsub("18097", "S18097", rowname)) %>% 
  mutate(rowname = gsub("-", ".", rowname))
Bifidobacterium_table

#take out rownames and add again (since row names disappear when joining)
clr_Bifidobacterium_joined <- left_join(clr_Bifidobacterium, Bifidobacterium_table, by = c("rowname" = "rowname"))

#change rownames
rownames(clr_Bifidobacterium_joined)<- clr_Bifidobacterium_joined[,"rowname"]

###### for RM #####
#select RM
bactRM <- clr_Bifidobacterium_joined %>%
  select(ends_with('viridae'), RM)
#pivot longer for plotting
bactRM_long <- bactRM %>% pivot_longer(-RM, names_to = "Virus", values_to = "clr_abund")

otulist<- colnames(bactRM %>% select(-RM))
pvalues <- rep(0, length(otulist))
donetable=data.frame(otulist, pvalues)
for(i in 1:length(otulist)){
  wilcoxtest<- subset(bactRM_long , Virus == otulist[i]) %>% 
    wilcox_test(clr_abund ~ RM) %>%
    add_significance()
  donetable$pvalues[i] = wilcoxtest$p
}
#add adjusted pvalues to dataset
p.adj <- p.adjust(donetable$pvalues)
donetableadj <- data.frame(donetable, p.adj)
bactRM_long_pvalues <- left_join(bactRM_long, donetableadj, by = c("Virus" = "otulist"))

library(ggpubr)
# Custom formatting function
format_pval <- function(pval){
  pval <- scales::pvalue(pval, accuracy= 0.0001, add_p = FALSE)
  gsub(pattern = "(=|>)", replacement = " ", x = pval)
}

#plot
#select significant ones from donetableadj
bifidobacteriumRMplot <-subset(bactRM_long_pvalues , p.adj < 0.05) %>% 
  ggplot(., aes(x=Virus, y=clr_abund, fill=as.factor(RM))) +
  geom_violin(alpha=0.5)+
  coord_flip()+
  #add median
  stat_summary(fun = median, geom = "point",
               position = position_dodge(0.9),
               color="#301934")+
  #add p-value
  geom_text(aes(y=8.5,label=paste0("p=",format_pval(pvalues),"\n","p.adjust=",format_pval(p.adj))),check_overlap = TRUE, size=2.6)+
  scale_fill_discrete(name = "RM", breaks=c('1', '0'), labels = c("Has system", "Does not have"))+
  ggtitle("CLR Transformed Relative Abundances")+
  theme(axis.text.y = element_text(size = 13),
        axis.title.y = element_blank())+
  scale_y_continuous(limits = c(NA,11))


###### for Abi2 #####
#select Abi2
bactAbi2 <- clr_Bifidobacterium_joined %>%
  select(ends_with('viridae'), Abi2)
#pivot longer for plotting
bactAbi2_long <- bactAbi2 %>% pivot_longer(-Abi2, names_to = "Virus", values_to = "clr_abund")

otulist<- colnames(bactAbi2 %>% select(-Abi2))
pvalues <- rep(0, length(otulist))
donetable=data.frame(otulist, pvalues)
for(i in 1:length(otulist)){
  wilcoxtest<- subset(bactAbi2_long , Virus == otulist[i]) %>% 
    wilcox_test(clr_abund ~ Abi2) %>%
    add_significance()
  donetable$pvalues[i] = wilcoxtest$p
}

p.adj <- p.adjust(donetable$pvalues)
donetableadj <- data.frame(donetable, p.adj)
bactAbi2_long_pvalues <- left_join(bactAbi2_long, donetableadj, by = c("Virus" = "otulist"))


#plot
#select significant ones from donetableadj
bifidobacteriumAbi2plot <- subset(bactAbi2_long_pvalues , p.adj < 0.05) %>% 
  ggplot(., aes(x=Virus, y=clr_abund, fill=as.factor(Abi2))) +
  geom_violin(alpha=0.5)+
  coord_flip()+
  stat_summary(fun = median, geom = "point",
               position = position_dodge(0.9),
               color="#301934")+
  geom_text(aes(y=9.5,label=paste0("p=",format_pval(pvalues),"\n","p.adjust=",format_pval(p.adj))),check_overlap = TRUE, size=2.6)+
  scale_fill_discrete(name = "Abi2", breaks=c('1', '0'), labels = c("Has system", "Does not have"))+
  ggtitle("CLR Transformed Relative Abundances")+
  theme(axis.text.y = element_text(size = 13),
        axis.title.y = element_blank())+
  scale_y_continuous(limits = c(NA,12))


###### for CasClass1 #####
#select CasClass1
bactCasClass1 <- clr_Bifidobacterium_joined %>%
  select(ends_with('viridae'), CasClass1)
#pivot longer for plotting
bactCasClass1_long <- bactCasClass1 %>% pivot_longer(-CasClass1, names_to = "Virus", values_to = "clr_abund")

otulist<- colnames(bactCasClass1 %>% select(-CasClass1))
pvalues <- rep(0, length(otulist))
donetable=data.frame(otulist, pvalues)
for(i in 1:length(otulist)){
  wilcoxtest<- subset(bactCasClass1_long , Virus == otulist[i]) %>% 
    wilcox_test(clr_abund ~ CasClass1) %>%
    add_significance()
  donetable$pvalues[i] = wilcoxtest$p
}
#add adjusted p-value
p.adj <- p.adjust(donetable$pvalues)
donetableadj <- data.frame(donetable, p.adj)
bactCasClass1_long_pvalues <- left_join(bactCasClass1_long, donetableadj, by = c("Virus" = "otulist"))

#plot
#select significant ones from donetableadj
bifidobacteriumCasClass1plot <- subset(bactCasClass1_long_pvalues , p.adj < 0.05) %>% 
  ggplot(., aes(x=Virus, y=clr_abund, fill=as.factor(CasClass1))) +
  geom_violin(alpha=0.5)+
  coord_flip()+
  stat_summary(fun = median, geom = "point",
               position = position_dodge(0.9),
               color="#301934")+
  geom_text(aes(y=9,label=paste0("p=",format_pval(pvalues),"\n","p.adjust=",format_pval(p.adj))),check_overlap = TRUE, size=2.6)+
  scale_fill_discrete(name = "Cas Class1",breaks=c('1', '0'), labels = c("Has system", "Does not have"))+
  ggtitle("CLR Transformed Relative Abundances")+
  theme(axis.text.y = element_text(size = 13),
        axis.title.y = element_blank())+
  scale_y_continuous(limits = c(NA,11))


###### for Wadjet #####
#select Wadjet
bactWadjet <- clr_Bifidobacterium_joined %>%
  select(ends_with('viridae'), Wadjet)
#pivot longer for plotting
bactWadjet_long <- bactWadjet %>% pivot_longer(-Wadjet, names_to = "Virus", values_to = "clr_abund")

otulist<- colnames(bactWadjet %>% select(-Wadjet))
pvalues <- rep(0, length(otulist))
donetable=data.frame(otulist, pvalues)
for(i in 1:length(otulist)){
  wilcoxtest<- subset(bactWadjet_long , Virus == otulist[i]) %>% 
    wilcox_test(clr_abund ~ Wadjet) %>%
    add_significance()
  donetable$pvalues[i] = wilcoxtest$p
  #print(paste0(wilcoxtest$p, " i=",otulist[i]))
}
#add adjusted p-value
p.adj <- p.adjust(donetable$pvalues)
donetableadj <- data.frame(donetable, p.adj)
bactWadjet_long_pvalues <- left_join(bactWadjet_long, donetableadj, by = c("Virus" = "otulist"))

#plot
#select significant ones from donetableadj
bifidobacteriumWadjetplot <- subset(bactWadjet_long_pvalues , p.adj < 0.05) %>% 
  ggplot(., aes(x=Virus, y=clr_abund, fill=as.factor(Wadjet))) +
  geom_violin(alpha=0.5)+
  coord_flip()+
  stat_summary(fun = median, geom = "point",
               position = position_dodge(0.9),
               color="#301934")+
  geom_text(aes(y=8,label=paste0("p=",format_pval(pvalues),"\n","p.adjust=",format_pval(p.adj))),check_overlap = TRUE, size=2.6)+
  scale_fill_discrete(name = "Wadjet", breaks=c('1', '0'), labels = c("Has system", "Does not have"),)+
  ggtitle("CLR Transformed Relative Abundances")+
  theme(axis.text.y = element_text(size = 13),
        axis.title.y = element_blank())+
  scale_y_continuous(limits = c(NA,11))



(bifidobacteriumAbi2plot | bifidobacteriumCasClass1plot)/
  (bifidobacteriumRMplot| bifidobacteriumWadjetplot)+
  plot_annotation(title = 'Differential Analysis for CLR Transformed Virus Abundances Targeting Bifidobacterium - Defense System Absence/Presence',
                  subtitle = "p-values adjusted with FDR method")

ggsave("MASTplots/fractions_bifidobacterium_transformedAndAgglomerated.jpeg", width=11.7, height=7)
