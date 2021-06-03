##Code for all analyses in Bell-Dereske and Emery "Contributions of environmental and maternal transmission to the assembly of leaf fungal endophyte communities"
##Note that figures were generated in R, but tweaked and finished in Inkscape


#####Libraries required####
here::i_am("R_code/Bell-Dereske&Evans_Rain_analyses_figures.R")
library(here)
library("phyloseq")
library("ggplot2")
library("vegan")
library(plyr); library(dplyr)
library(reshape2)
library(tidyr)
library(stringr)
library(seqinr)
library(betapart)
library(otuSummary)
library(emmeans)
library(lme4)
library(lmerTest)
#library(ranacapa)
library(cowplot)
library(RColorBrewer)
library(ggh4x)
library(car)
library(indicspecies)
pa=function(x)(ifelse(x>0,1,0))





#####Preprocessing of Community Data####
FunRainLeaf2019.MAP=read.csv(here::here("USEARCHv11","Rain_Leaf_GLBRC_Marsh_metadata.csv"), header = T, row.names = "sampleID_fung")
head(FunRainLeaf2019.MAP)
nrow(FunRainLeaf2019.MAP)
#200


#Read in the fasta file

rep_set.FunRainLeaf2019_ZOTU_full<- read.fasta(here::here("USEARCHv11","rep_set_ITS2_full_FungiRainLeaf2019_zotus.fa"), as.string = TRUE, set.attributes = FALSE)
head(rep_set.FunRainLeaf2019_ZOTU_full)
length(rep_set.FunRainLeaf2019_ZOTU_full)
#18466


#ASV table
FunRainLeaf2019.ZOTU_full=read.delim(here::here("USEARCHv11","ZOTU_table_ITS2_full_FungiRainLeaf2019.txt"),header=T, row.names = 1)
head(FunRainLeaf2019.ZOTU_full)
colnames(FunRainLeaf2019.ZOTU_full)
head(rownames(FunRainLeaf2019.ZOTU_full))

ncol(FunRainLeaf2019.ZOTU_full)
#200
nrow(FunRainLeaf2019.ZOTU_full)
#18465

FunRainLeaf2019_ZOTU_full = otu_table(FunRainLeaf2019.ZOTU_full, taxa_are_rows = TRUE)
FunRainLeaf2019_ZOTU_full[1:10,1:10]

#UNITE8 Taxon table

FunRainLeaf2019.taxa_ZOTU_full_raw_UNITE8 = read.delim(here::here("USEARCHv11","taxonomy_ITS2_full_FungiRainLeaf2019_zotus.sintax"),sep = c("\t"),header = F)
head(FunRainLeaf2019.taxa_ZOTU_full_raw_UNITE8)
nrow(FunRainLeaf2019.taxa_ZOTU_full_raw_UNITE8)
#18466
head(FunRainLeaf2019.taxa_ZOTU_full_raw_UNITE8)
FunRainLeaf2019.taxa_ZOTU_full_raw_80c_UNITE8=FunRainLeaf2019.taxa_ZOTU_full_raw_UNITE8[,c(1,4)]
head(FunRainLeaf2019.taxa_ZOTU_full_raw_80c_UNITE8)

FunRainLeaf2019.taxa_ZOTU_full_raw_80c_UNITE8_sep=FunRainLeaf2019.taxa_ZOTU_full_raw_80c_UNITE8 %>% separate(V4, c("Domain","Phylum","Class","Order","Family","Genus","Species"),sep = ",")
row.names(FunRainLeaf2019.taxa_ZOTU_full_raw_80c_UNITE8_sep)=FunRainLeaf2019.taxa_ZOTU_full_raw_80c_UNITE8_sep$V1
FunRainLeaf2019.taxa_ZOTU_full_raw_80c_UNITE8_sep$V1=NULL
FunRainLeaf2019.taxa_ZOTU_full_raw_80c_UNITE8_sep[is.na(FunRainLeaf2019.taxa_ZOTU_full_raw_80c_UNITE8_sep)] <- "UNKNOWN"
FunRainLeaf2019.taxa_ZOTU_full_raw_80c_UNITE8_sep_mat=as.matrix(FunRainLeaf2019.taxa_ZOTU_full_raw_80c_UNITE8_sep)
head(FunRainLeaf2019.taxa_ZOTU_full_raw_80c_UNITE8_sep_mat)
unique(data.frame(FunRainLeaf2019.taxa_ZOTU_full_raw_80c_UNITE8_sep_mat)$Domain)
TAXA_80C_UNITE8_FunRainLeaf2019_ZOTU_full=tax_table(FunRainLeaf2019.taxa_ZOTU_full_raw_80c_UNITE8_sep_mat)
head(TAXA_80C_UNITE8_FunRainLeaf2019_ZOTU_full)

#Make the phyloseq object

FunRainLeaf2019_ZOTU_full.data=phyloseq(FunRainLeaf2019_ZOTU_full,sample_data(FunRainLeaf2019.MAP),TAXA_80C_UNITE8_FunRainLeaf2019_ZOTU_full)
ntaxa(FunRainLeaf2019_ZOTU_full.data)
#18465
sum(otu_table(FunRainLeaf2019_ZOTU_full.data))
#6988341
max(sample_sums(FunRainLeaf2019_ZOTU_full.data))
#65728
min(sample_sums(FunRainLeaf2019_ZOTU_full.data))
#15


#Remove Non-Fungal ASVs
FunRainLeaf2019_ZOTU_full.fung<-subset_taxa(FunRainLeaf2019_ZOTU_full.data,Domain=="d:Fungi")
ntaxa(FunRainLeaf2019_ZOTU_full.fung)
#15785
sum(otu_table(FunRainLeaf2019_ZOTU_full.fung))
#4916422
max(sample_sums(FunRainLeaf2019_ZOTU_full.fung))
#65077
#paired 57064
#single 48721
#zZOTU_full 59455
min(sample_sums(FunRainLeaf2019_ZOTU_full.fung))
#3
sort(sample_sums(FunRainLeaf2019_ZOTU_full.fung))#Leaf55 and Leaf25 are outliers
#Leaf33        Leaf50        Leaf17  RX026GLBRC26        Leaf64 BlankRainPCR2        Leaf72        Leaf58  RX037GLBRC37  RX065GLBRC65 
#3            13            16            33            66           202           610           675           781          1015 
#Leaf63        Leaf07        Leaf05        Leaf74  RX063GLBRC63        Leaf61        Leaf16        Leaf41        Leaf32        Leaf09 
#1074          1350          1378          1713          1744          1747          1762          1874          2078          2132



#I am going to use CONSTAX for taxonomy since it estimates deeper classifications and catches some crap sequences
FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8.2_all = read.delim(here::here("USEARCHv11",
                                                                          "ZOTU_constax_V2_sintax_fix_classification_all_v4.2.2020_taxa",
                                                                          "combined_taxonomy.txt"),
                                                                 sep = "\t",header = T,fill=T, row.names = 1)
head(FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8.2_all)
nrow(FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8.2_all)
#18295

FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8.2_all[FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8.2_all==""|
                                                          is.na(FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8.2_all)]="Unknown"
head(FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8.2_all)

#There is a formating error where sintax kingdom is a number if there are no lower classifications
#I am going to replace any numbers in Kingdom_SINTAX with the classification from the raw sintax taxonomy

FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8.2_SINTAX = 
  read.delim(here::here("USEARCHv11","ZOTU_constax_V2_sintax_fix_classification_all_v4.2.2020_taxa/otu_taxonomy.SINTAX"),
                                                                     sep = "\t",header = F,fill=T, row.names = 1)

head(FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8.2_SINTAX)
nrow(FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8.2_SINTAX)
#18466
FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8.2_SINTAX_V4=FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8.2_SINTAX
FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8.2_SINTAX_V4[,c("V2","V3")]=NULL
head(FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8.2_SINTAX_V4)


FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8.2_SINTAX_V4$V4=str_replace(FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8.2_SINTAX_V4$V4,"d:","")

FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8.2_all_v2=merge(FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8.2_all,
                                                                 FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8.2_SINTAX_V4,
                                                                 by="row.names")
head(FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8.2_all_v2)
unique(FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8.2_all_v2$Kingdom_SINTAX)
FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8.2_all_v2$Kingdom_SINTAX=ifelse(FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8.2_all_v2$Kingdom_SINTAX=="Viridiplantae"| 
                                                                                 FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8.2_all_v2$Kingdom_SINTAX=="Unknown"|
                                                                                   FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8.2_all_v2$Kingdom_SINTAX=="Fungi"|
                                                                                   FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8.2_all_v2$Kingdom_SINTAX=="Alveolata"|
                                                                                   FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8.2_all_v2$Kingdom_SINTAX=="Heterolobosa"|
                                                                                   FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8.2_all_v2$Kingdom_SINTAX=="Amoebozoa"|
                                                                                   FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8.2_all_v2$Kingdom_SINTAX=="Stramenopila"|
                                                                                   FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8.2_all_v2$Kingdom_SINTAX=="Metazoa",
                                                                                 FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8.2_all_v2$Kingdom_SINTAX,
                                                                                 FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8.2_all_v2$V4)
unique(FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8.2_all_v2$Kingdom_SINTAX)
FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8.2_all_v2$Kingdom_SINTAX=str_replace(FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8.2_all_v2$Kingdom_SINTAX,
                                                                                      "unidentified","Unknown")

unique(FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8.2_all_v2$Kingdom_Consensus)
FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8.2_all_v2$Kingdom_Consensus=ifelse(FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8.2_all_v2$Kingdom_Consensus=="Viridiplantae"| 
                                                                                   FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8.2_all_v2$Kingdom_Consensus=="Unknown"|
                                                                                   FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8.2_all_v2$Kingdom_Consensus=="Fungi"|
                                                                                   FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8.2_all_v2$Kingdom_Consensus=="Alveolata"|
                                                                                   FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8.2_all_v2$Kingdom_Consensus=="Heterolobosa"|
                                                                                   FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8.2_all_v2$Kingdom_Consensus=="Amoebozoa"|
                                                                                   FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8.2_all_v2$Kingdom_Consensus=="Stramenopila"|
                                                                                   FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8.2_all_v2$Kingdom_Consensus=="Metazoa"|
                                                                                     FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8.2_all_v2$Kingdom_Consensus=="Protista"|
                                                                                     FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8.2_all_v2$Kingdom_Consensus=="Rhizaria",
                                                                                 FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8.2_all_v2$Kingdom_Consensus,
                                                                                 FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8.2_all_v2$V4)



unique(FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8.2_all_v2$Kingdom_Consensus)
FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8.2_all_v2$Kingdom_Consensus=str_replace(FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8.2_all_v2$Kingdom_Consensus,
                                                                                      "unidentified","Unknown")

#let's look into the strength of classifications

FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8.2_all_v2$Phyla_unkn=ifelse(FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8.2_all_v2$Phylum_RDP=="Unknown"|
                                                                            FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8.2_all_v2$Phylum_BLAST=="Unknown"|
                                                                          FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8.2_all_v2$Phylum_SINTAX=="Unknown", "TRUE","FALSE")


nrow(subset(FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8.2_all_v2,Phyla_unkn==TRUE))
#6546


#let's look at the OTUs that passed the first round of filtering

ntaxa(FunRainLeaf2019_ZOTU_full.fung)
#15785

head(FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8.2_all)
FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw=merge(data.frame("V2.2_OTUs"=taxa_names(FunRainLeaf2019_ZOTU_full.fung),
                                                                                  "V2.2_Phylum"=data.frame(tax_table(FunRainLeaf2019_ZOTU_full.fung))$Phylum),
                                                                       FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8.2_all_v2,by.x="V2.2_OTUs", by.y = "Row.names", all.x = T)
nrow(FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw)
#15785
head(FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw)
#Let's fix the formaatting so the classifications are the same

FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw$V2.2_Phylum=str_replace(FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw$V2.2_Phylum,"p:","")
FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw$V2.2_Phylum=str_replace(FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw$V2.2_Phylum,"UNKNOWN","Unknown")
FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw[is.na(FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw)] <- "Unknown"
unique(FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw$Kingdom_Consensus)

#let's look into the strength of classifications

FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw$Phyla_unkn=ifelse(FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw$Phylum_RDP=="Unknown"|
                                                                            FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw$Phylum_BLAST=="Unknown"|
                                                                            FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw$Phylum_SINTAX=="Unknown", "TRUE","FALSE")


nrow(subset(FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw,Phyla_unkn==TRUE))
#5076

#How many TAXA have a classification at Phyla across the three measures but differ in there phyla classification

FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw_KNW=subset(FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw,Phyla_unkn==FALSE)
nrow(FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw_KNW)
#10709


FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw_KNW$Phyla_diff=ifelse(FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw_KNW$Phylum_RDP!=
                                                                                FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw_KNW$Phylum_BLAST|
                                                                                FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw_KNW$Phylum_RDP!=
                                                                                FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw_KNW$Phylum_SINTAX|
                                                                                FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw_KNW$Phylum_BLAST!=
                                                                                FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw_KNW$Phylum_SINTAX, "TRUE","FALSE")




nrow(subset(FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw_KNW,Phyla_diff==TRUE))
#13

FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw_KNW_P_dif=subset(FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw_KNW,Phyla_diff==TRUE)



sum(ifelse(FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw_KNW_P_dif$Phylum_RDP==
             FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw_KNW_P_dif$Phylum_BLAST,0,1))
#3
sum(ifelse(FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw_KNW_P_dif$Phylum_SINTAX==
             FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw_KNW_P_dif$Phylum_BLAST,0,1))
#13
sum(ifelse(FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw_KNW_P_dif$Phylum_RDP==
             FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw_KNW_P_dif$Phylum_SINTAX,0,1))
#10

sum(ifelse(FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw_KNW_P_dif$Phylum_SINTAX!=
             FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw_KNW_P_dif$Phylum_BLAST&
             FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw_KNW_P_dif$Phylum_SINTAX==
             FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw_KNW_P_dif$Phylum_RDP,1,0))
#3

unique(with(FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw_KNW_P_dif,interaction(Phylum_SINTAX,Phylum_RDP,Phylum_BLAST)))


data.frame(FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw_KNW_P_dif) %>% group_by(Phylum_SINTAX,Phylum_BLAST,Phylum_RDP)%>%summarise_at(vars(Kingdom_BLAST),~n())

# A tibble: 3 x 4
# Groups:   Phylum_SINTAX, Phylum_BLAST [2]
#Phylum_SINTAX Phylum_BLAST    Phylum_RDP    Kingdom_BLAST
#<chr>         <chr>           <chr>                 <int>
#1 Ascomycota    Basidiomycota Ascomycota                2
#2 Basidiomycota Ascomycota    Ascomycota               10
#3 Basidiomycota Ascomycota    Basidiomycota             1

#How many of these sintax classifications match sintax that classified against v02.02.2019

sum(ifelse(FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw_KNW_P_dif$V2.2_Phylum==
             FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw_KNW_P_dif$Phylum_SINTAX,0,1))

#0
FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw_KNW_P_dif[,c("V2.2_Phylum","Phylum_SINTAX")]
#total number of phyla classification that differ with unknowns included

FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw$Phyla_diff=ifelse(FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw$Phylum_RDP!=
                                                                            FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw$Phylum_BLAST|
                                                                            FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw$Phylum_RDP!=
                                                                            FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw$Phylum_SINTAX|
                                                                            FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw$Phylum_BLAST!=
                                                                            FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw$Phylum_SINTAX, "TRUE","FALSE")
nrow(subset(FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw,Phyla_diff==TRUE))
#3779

FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw_P_dif=subset(FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw,Phyla_diff==TRUE)



sum(ifelse(FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw_P_dif$Phylum_RDP==
             FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw_P_dif$Phylum_BLAST,0,1))
#1276
sum(ifelse(FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw_P_dif$Phylum_SINTAX==
             FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw_P_dif$Phylum_BLAST,0,1))
#3411
sum(ifelse(FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw_P_dif$Phylum_RDP==
             FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw_P_dif$Phylum_SINTAX,0,1))
#2873

sum(ifelse(FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw_P_dif$Phylum_SINTAX!=
             FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw_P_dif$Phylum_BLAST&
             FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw_P_dif$Phylum_SINTAX==
             FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw_P_dif$Phylum_RDP,1,0))


#906

#I am going to use the concensus for all of the times when the three Phyla match
FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw_KNW_NO_P_dif=subset(FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw_KNW,Phyla_diff==FALSE)
nrow(FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw_KNW_NO_P_dif)
head(FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw_KNW_NO_P_dif)
FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_KNW_NO_P_dif=
  FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw_KNW_NO_P_dif[,c("V2.2_OTUs","Kingdom_Consensus","Phylum_Consensus","Class_Consensus","Order_Consensus",
                                                                                   "Family_Consensus","Genus_Consensus","Species_Consensus")]
row.names(FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_KNW_NO_P_dif)=FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_KNW_NO_P_dif$V2.2_OTUs

FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_KNW_NO_P_dif$V2.2_OTUs=NULL
colnames(FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_KNW_NO_P_dif)=c("Kingdom","Phylum","Class","Order","Family","Genus","Species")

head(FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_KNW_NO_P_dif)
#I am going to fix the miss matches by hand

FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw_KNW_P_dif=subset(FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw_KNW,Phyla_diff==TRUE)



FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_KNW_P_dif=FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw_KNW_P_dif[c("V2.2_OTUs","Kingdom_Consensus","Phylum_Consensus","Class_Consensus","Order_Consensus",
                                                                                                "Family_RDP","Genus_RDP","Species_RDP")]

row.names(FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_KNW_P_dif)=FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_KNW_P_dif$V2.2_OTUs

FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_KNW_P_dif$V2.2_OTUs=NULL
colnames(FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_KNW_P_dif)=c("Kingdom","Phylum","Class","Order","Family","Genus","Species")

head(FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_KNW_P_dif)

#Classifications that have at least one unknown from the classifiers

FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw_UNKNW=subset(FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw,Phyla_unkn==TRUE)
nrow(FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw_UNKNW)
#5076
unique(FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw_UNKNW$Kingdom_RDP)
unique(FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw_UNKNW$Kingdom_SINTAX)
unique(FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw_UNKNW$Kingdom_BLAST)
#how many times in Blast Unknown and SINTAX diffs from RDP 
sum(ifelse(FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw_UNKNW$Kingdom_BLAST=="Unknown"&
         FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw_UNKNW$Kingdom_SINTAX!=
         FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw_UNKNW$Kingdom_RDP,1,0))
#145
#how many times in RDP Unknown and SINTAX diffs from BLAST
sum(ifelse(FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw_UNKNW$Kingdom_RDP=="Unknown"&
             FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw_UNKNW$Kingdom_SINTAX!=
             FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw_UNKNW$Kingdom_BLAST,1,0))
#282
#how many times in SINTAX Unknown and RDP diffs from BLAST
sum(ifelse(FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw_UNKNW$Kingdom_SINTAX=="Unknown"&
             FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw_UNKNW$Kingdom_RDP!=
             FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw_UNKNW$Kingdom_BLAST,1,0))
#173

#I am going to use the concensus for all the times that there is one Unknown match since I am not sure how to check this
head(FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw_UNKNW)
FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_UNKNW=
  FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_Taxa_tbl_raw_UNKNW[,c("V2.2_OTUs","Kingdom_Consensus","Phylum_Consensus","Class_Consensus","Order_Consensus",
                                                                                   "Family_Consensus","Genus_Consensus","Species_Consensus")]
row.names(FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_UNKNW)=FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_UNKNW$V2.2_OTUs

FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_UNKNW$V2.2_OTUs=NULL
colnames(FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_UNKNW)=c("Kingdom","Phylum","Class","Order","Family","Genus","Species")

head(FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_UNKNW)

#Combining and formatting for phyloseq
FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_mat=as.matrix(rbind(FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_KNW_NO_P_dif,
                                                          FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_KNW_P_dif,
                                                          FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_all_UNKNW))

head(FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_mat)
unique(data.frame(FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_mat)$Kingdom)
nrow(FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_mat)
#15785
FunRainLeaf2019_ZOTU_full.fung_cons=phyloseq(tax_table(as.matrix(FunRainLeaf2019_ZOTU_full.fung_CONSTAX_UNITE8.2_mat)),otu_table(FunRainLeaf2019_ZOTU_full.fung),
                                                    sample_data(FunRainLeaf2019_ZOTU_full.fung))
ntaxa(FunRainLeaf2019_ZOTU_full.fung_cons)
#15785


#Let's look at the taxa that are not classified to fungi using CONSTAX

sum(taxa_sums(subset_taxa(FunRainLeaf2019_ZOTU_full.fung_cons,Kingdom=="Unknown")))
#49825
unknown_kingdom_names=taxa_names(subset_taxa(FunRainLeaf2019_ZOTU_full.fung_cons,Kingdom=="Unknown"))
length(unknown_kingdom_names)
#253


head(rep_set.FunRainLeaf2019_ZOTU_full)
length(rep_set.FunRainLeaf2019_ZOTU_full)
#18466

rep_set.FunRainLeaf2019_ZOTU_full.fung_cons_UNKNOWN=rep_set.FunRainLeaf2019_ZOTU_full[names(rep_set.FunRainLeaf2019_ZOTU_full) %in% unknown_kingdom_names]
head(rep_set.FunRainLeaf2019_ZOTU_full.fung_cons_UNKNOWN)
length(rep_set.FunRainLeaf2019_ZOTU_full.fung_cons_UNKNOWN)
#253


write.fasta(sequences =rep_set.FunRainLeaf2019_ZOTU_full.fung_cons_UNKNOWN, names = names(rep_set.FunRainLeaf2019_ZOTU_full.fung_cons_UNKNOWN), 
            file.out =here::here("R_file","rep_set.FunRainLeaf2019_ZOTU_full.fung_cons_UNKNOWN_UKNOWN_v4.2.2020.fna"))

#None of these taxa match up well with NCBI Fungi

#Let's remove them 


FunRainLeaf2019_ZOTU_full.fung_cons<-subset_taxa(FunRainLeaf2019_ZOTU_full.fung_cons,Kingdom=="Fungi")
ntaxa(FunRainLeaf2019_ZOTU_full.fung_cons)
#15518

sum(otu_table(FunRainLeaf2019_ZOTU_full.fung_cons))
#4865434

max(sample_sums(FunRainLeaf2019_ZOTU_full.fung_cons))
#59235

min(sample_sums(FunRainLeaf2019_ZOTU_full.fung_cons))
#3

sort(sample_sums(FunRainLeaf2019_ZOTU_full.fung_cons))
#Leaf33        Leaf50        Leaf17  RX026GLBRC26        Leaf64 BlankRainPCR2        Leaf72        Leaf58  RX037GLBRC37  RX065GLBRC65        Leaf63        Leaf07        Leaf05 
#3            13            14            32            65           202           593           670           781          1013          1074          1341          1375 
#Leaf16        Leaf74        Leaf57        Leaf32  RX063GLBRC63        Leaf61        Leaf41        Leaf09        Leaf73        Leaf70        Leaf14        Leaf69        Leaf21 
#1396          1476          1477          1553          1675          1747          1762          2076          2298          2573          2852          2944          2991 
#Leaf56        Leaf45        Leaf13        Leaf34        Leaf30        Leaf19        Leaf25        Leaf66        Leaf22        Leaf71        Leaf27        Leaf10        Leaf03 
#3061          3197          3225          3270          3574          3906          4188          4586          4818          5230          5363          5772          6022 
#Leaf36        Leaf06        Leaf40        Leaf38        Leaf48        Leaf54        Leaf49        Leaf04        Leaf46        Leaf52  RX054Blank54        Leaf31        Leaf18 
#6360          6958          7231          7341          7414          7916          8529          8683          8718          9461          9825         10132         10916



#Remove Malasseziomycetes

FunRainLeaf2019_ZOTU_full.fung_cons<-subset_taxa(FunRainLeaf2019_ZOTU_full.fung_cons,Class!="Malasseziomycetes")
ntaxa(FunRainLeaf2019_ZOTU_full.fung_cons)
#15504
sum(otu_table(FunRainLeaf2019_ZOTU_full.fung_cons))
#4864405
max(sample_sums(FunRainLeaf2019_ZOTU_full.fung_cons))
#59235
min(sample_sums(FunRainLeaf2019_ZOTU_full.fung_cons))
#3
sort(sample_sums(FunRainLeaf2019_ZOTU_full.fung_cons))




#####FUNGuild CONSTAX Dataset creation####

#Example formating

#OTU ID	sample1	sample2	sample3	sample4	sample5	taxonomy
#OTU_100	0	1	0	0	0	93.6%|Laetisaria_fuciformis|EU118639|SH012042.06FU|reps_singleton|k__Fungi;p__Basidiomycota;c__Agaricomycetes;o__Corticiales;f__Corticiaceae;g__Laetisaria;s__Laetisaria_fuciformis
FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8 = data.frame(tax_table(FunRainLeaf2019_ZOTU_full.fung_cons))
head(FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8)
nrow(FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8)
#15504
#remove all of the unknown
FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8[FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8=="Unknown"|
                                                          is.na(FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8)]=""

FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8$taxonomy=paste("k__",FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8$Kingdom,";",
                                                                 "p__",FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8$Phylum,";",
                                                                 "c__",FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8$Class,";",
                                                                 "o__",FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8$Order,";",
                                                                 "f__",FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8$Family,";",
                                                                 "g__",FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8$Genus,";",
                                                                 "s__",FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8$Species,sep = "")
FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8$taxonomy=str_replace_all(FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8$taxonomy, " ", "_")
FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8$taxonomy=str_replace_all(FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8$taxonomy, ";p__;c__;o__;f__;g__;s__", "")
FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8$taxonomy=str_replace_all(FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8$taxonomy, ";c__;o__;f__;g__;s__", "")
FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8$taxonomy=str_replace_all(FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8$taxonomy, ";o__;f__;g__;s__", "")
FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8$taxonomy=str_replace_all(FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8$taxonomy, ";f__;g__;s__", "")
FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8$taxonomy=str_replace_all(FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8$taxonomy, ";g__;s__", "")
FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8$taxonomy=str_replace_all(FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8$taxonomy, ";s__\n", "")
FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8$taxonomy=str_replace_all(FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8$taxonomy, " ", "_")
head(FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8$taxonomy)
FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8$OTU_ID=row.names(FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8)
FunRainLeaf2019_ZOTU_CONSTAX_ZOTU=merge(otu_table(FunRainLeaf2019_ZOTU_full.fung_cons),FunRainLeaf2019.taxa_ZOTU_full_raw_CONSTAX_UNITE8[,c("OTU_ID","taxonomy")], 
                                    by="row.names",by.y="OTU_ID")
head(FunRainLeaf2019_ZOTU_CONSTAX_ZOTU)


colnames(FunRainLeaf2019_ZOTU_CONSTAX_ZOTU)[1]="OTUID"


colnames(FunRainLeaf2019_ZOTU_CONSTAX_ZOTU)

write.csv(FunRainLeaf2019_ZOTU_CONSTAX_ZOTU, here::here("R_file","FunRainLeaf2019_ZOTU_CONSTAXv4.2.2020_ZOTU_for_FunGuild.csv"), row.names = F) 
#I had to modify this in excel (i.e. turn "OTU.ID "to "OTU ID")

#RUN in shell 

#cd HardDrive/FungiRainLeaf2019/R_file/
#python Guilds_v1.1.py -otu FunRainLeaf2019_ZOTU_CONSTAXv4.2.2020_ZOTU_for_FunGuild.csv -db fungi -m -u

#Found 12894 matching taxonomy records in the database.
#Dereplicating and sorting the result...
#FunGuild tried to assign function to 15504 OTUs in 'FunRainLeaf2019_ZOTU_CONSTAXv4.2.2020_ZOTU_for_FunGuild.csv'.
#FUNGuild made assignments on 8307 OTUs.
#Result saved to 'FunRainLeaf2019_ZOTU_CONSTAXv4.2.2020_ZOTU_for_FunGuild.guilds.txt'

#Additional output:
#  FUNGuild made assignments on 8307 OTUs, these have been saved to FunRainLeaf2019_ZOTU_CONSTAXv4.2.2020_ZOTU_for_FunGuild.guilds_matched.txt.
#7197 OTUs were unassigned, these are saved to FunRainLeaf2019_ZOTU_CONSTAXv4.2.2020_ZOTU_for_FunGuild.guilds_unmatched.txt.

#Total calculating time: 185.96 seconds.







#####Preprocessing####

#####Leaf Contaminant removal####


#First Let's look at the leaf communities 
head(sample_data(FunRainLeaf2019_ZOTU_full.fung_cons))
unique(sample_data(FunRainLeaf2019_ZOTU_full.fung_cons)$extract_kit)
Mar_leaf.fung=subset_samples(FunRainLeaf2019_ZOTU_full.fung_cons,sample_type=="Leaf"|extract_kit=="PCR")#extract_kit=="Mock"
unique(sample_data(Mar_leaf.fung)$Sample_or_Control)
Mar_leaf.fung<-prune_taxa(taxa_sums(Mar_leaf.fung) > 0, Mar_leaf.fung)
nsamples(Mar_leaf.fung)
#75
ntaxa(Mar_leaf.fung)
#2456
sum(otu_table(Mar_leaf.fung))
#876222
max(sample_sums(Mar_leaf.fung))
#37843
min(sample_sums(Mar_leaf.fung))
#3
sort(sample_sums(Mar_leaf.fung))
#Leaf33        Leaf50        Leaf17        Leaf64 BlankRainPCR2        Leaf72        Leaf58        Leaf63        Leaf07        Leaf05        Leaf16        Leaf74        Leaf57 
#3            13            14            65           202           593           670          1074          1341          1375          1396          1469          1473 
#Leaf32        Leaf61        Leaf41        Leaf09        Leaf73        Leaf70        Leaf14        Leaf69        Leaf21        Leaf56        Leaf45        Leaf13        Leaf34 
#1553          1747          1760          2076          2295          2573          2849          2944          2985          3055          3196          3225          3270 
#Leaf30        Leaf19        Leaf25        Leaf66        Leaf22        Leaf71        Leaf27        Leaf10        Leaf03        Leaf36        Leaf06        Leaf40        Leaf38 
#3574          3858          4149          4586          4816          5223          5363          5765          6022          6360          6958          7231          7341 
#Leaf48        Leaf54        Leaf49        Leaf04        Leaf46        Leaf52        Leaf31        Leaf18        Leaf44        Leaf60        Leaf28        Leaf42 BlankRainPCR1 
#7414          7916          8528          8668          8718          9461         10132         10913         10944         11559         12877         13206         15129 
#Leaf35        Leaf11        Leaf65        Leaf08        Leaf26        Leaf53        Leaf51        Leaf47        Leaf43        Leaf68        Leaf39        Leaf29        Leaf37 
#18460         18967         19238         19857         21017         21297         22148         22911         23842         23927         24216         26147         28019 
#Leaf67        Leaf20        Leaf12        Leaf15        Leaf59        Leaf75        Leaf24        Leaf23        Leaf55        Leaf62 
#28350         29446         29917         32304         33187         35557         36253         36693         36699         37843 




#Let's remove any contaminants from the OTU table before we move forward
#I am going to use https://github.com/donaldtmcknight/microDecon
library(microDecon)



Mar_leaf.fung_ZOTU=data.frame(otu_table(Mar_leaf.fung))
#blanks need to be the first columns 
head(sample_data(Mar_leaf.fung_ZOTU))
L.Blanks=row.names(sample_data(subset_samples(Mar_leaf.fung, Sample_or_Control=="Control")))
length(L.Blanks)
#8

unique(sample_data(Mar_leaf.fung)$extract_kit)
#I also need to group by sample/extraction orders

Group4 = row.names(sample_data(subset_samples(Mar_leaf.fung, extract_kit=="Plant_1"&Sample_or_Control!="Control")))
length(Group4)
#45
Group5 = row.names(sample_data(subset_samples(Mar_leaf.fung, extract_kit=="Plant_2"&Sample_or_Control!="Control")))
length(Group5)
#22
#Group6 = row.names(sample_data(subset_samples(Mar_leaf.fung, extract_kit=="Mock")))
#length(Group6)
#3


Mar_leaf.fung_ZOTU$OTUs=row.names(Mar_leaf.fung_ZOTU)
colnames(Mar_leaf.fung_ZOTU)
Mar_leaf.fung_ZOTU_reorder=Mar_leaf.fung_ZOTU[,c("OTUs",L.Blanks,Group4,Group5)]#,Group6
colnames(Mar_leaf.fung_ZOTU_reorder)
nrow(Mar_leaf.fung_ZOTU_reorder)
#2456



Mar_leaf.fung_ZOTU_reorder_decon_obj=decon(Mar_leaf.fung_ZOTU_reorder,numb.blanks = 8,numb.ind=c(45,22), taxa=F)#,3
nrow(Mar_leaf.fung_ZOTU_reorder_decon_obj$decon.table)
#2241

L.decon_names=Mar_leaf.fung_ZOTU_reorder_decon_obj$OTUs.removed$OTUs
length(L.decon_names)
#260



Mar_leaf.fung_decon_table=Mar_leaf.fung_ZOTU_reorder_decon_obj$decon.table

head(Mar_leaf.fung_decon_table)

#reformat for use in phyloseq

row.names(Mar_leaf.fung_decon_table)=Mar_leaf.fung_decon_table$OTUs
Mar_leaf.fung_decon_table[,c("OTUs", "Mean.blank")]=NULL
head(Mar_leaf.fung_decon_table)

Mar_leaf.fung_decon=phyloseq(otu_table(Mar_leaf.fung_decon_table, taxa_are_rows = T), sample_data(Mar_leaf.fung), tax_table(Mar_leaf.fung))
nsamples(Mar_leaf.fung_decon)
#67
ntaxa(Mar_leaf.fung_decon)
#2241

sum(otu_table(Mar_leaf.fung_decon))
#685304

max(sample_sums(Mar_leaf.fung_decon))
#35971

min(sample_sums(Mar_leaf.fung_decon))
#1

sort(sample_sums(Mar_leaf.fung_decon))
#Leaf33 Leaf17 Leaf64 Leaf72 Leaf58 Leaf63 Leaf07 Leaf05 Leaf57 Leaf16 Leaf32 Leaf61 Leaf41 Leaf73 Leaf09 Leaf70 Leaf14 Leaf21 Leaf69 Leaf56 Leaf13 Leaf45 Leaf25 Leaf19 Leaf34 Leaf30 
#1      6     51    488    632    877   1150   1203   1209   1349   1434   1623   1630   2009   2010   2344   2643   2712   2859   2962   2993   3038   3073   3104   3118   3240 
#Leaf71 Leaf66 Leaf22 Leaf03 Leaf27 Leaf10 Leaf36 Leaf04 Leaf38 Leaf40 Leaf48 Leaf49 Leaf54 Leaf52 Leaf46 Leaf31 Leaf18 Leaf44 Leaf60 Leaf28 Leaf42 Leaf15 Leaf35 Leaf11 Leaf65 Leaf08 
#4426   4546   4557   4803   5013   5318   5601   6703   6834   6886   7129   7626   7738   8020   8068   9376   9638   9834  10929  11092  12841  15458  17077  17522  18904  19381 
#Leaf47 Leaf53 Leaf43 Leaf39 Leaf26 Leaf12 Leaf68 Leaf37 Leaf29 Leaf67 Leaf24 Leaf23 Leaf59 Leaf75 Leaf55 
#19824  20524  20771  20904  20927  22770  23089  25616  26041  27223  29671  32297  32418  34180  35971 



Mar_leaf.fung_decon_pr=prune_samples(sample_sums(Mar_leaf.fung_decon) > 1000, Mar_leaf.fung_decon)
Mar_leaf.fung_decon_pr=prune_taxa(taxa_sums(Mar_leaf.fung_decon_pr) > 0, Mar_leaf.fung_decon_pr)


nsamples(Mar_leaf.fung_decon_pr)
#61
ntaxa(Mar_leaf.fung_decon_pr)
#2071

sum(otu_table(Mar_leaf.fung_decon_pr))
#683249

max(sample_sums(Mar_leaf.fung_decon_pr))
#35971

min(sample_sums(Mar_leaf.fung_decon_pr))
#1150


#####Leaf Contaminant removal####




#####Rain Contaminant removal####


#Now let's look at rain
head(sample_data(FunRainLeaf2019_ZOTU_full.fung_cons))
unique(sample_data(FunRainLeaf2019_ZOTU_full.fung_cons)$extract_kit)
unique(sample_data(FunRainLeaf2019_ZOTU_full.fung_cons)$sample_type)
Mar_rain.fungi=subset_samples(FunRainLeaf2019_ZOTU_full.fung_cons,sample_type=="Rain"|extract_kit=="PCR")#|extract_kit=="Mock"
unique(sample_data(FunRainLeaf2019_ZOTU_full.fung_cons)$Sample_or_Control)
Mar_rain.fungi<-prune_taxa(taxa_sums(Mar_rain.fungi) > 0, Mar_rain.fungi)
nsamples(Mar_rain.fungi)
#124
ntaxa(Mar_rain.fungi)
#14921
sum(otu_table(Mar_rain.fungi))
#3907117
max(sample_sums(Mar_rain.fungi))
#59235
min(sample_sums(Mar_rain.fungi))
#32
sort(sample_sums(Mar_rain.fungi))
#RX026GLBRC26 BlankRainPCR2  RX037GLBRC37  RX065GLBRC65  RX063GLBRC63  RX054Blank54 BlankRainPCR1  RX064Marsh64  RX032GLBRC32  RX044GLBRC44   RX005Marsh5  RX076Marsh76  RX093Marsh93 
#32           202           781          1007          1665          9614         15129         17348         17711         18198         19112         19224         19253
head(sample_data(Mar_rain.fungi))


#####microDecon

#Let try to remove any contaminants from the OTU table before we move forward
#I am going to use https://github.com/donaldtmcknight/microDecon
library(microDecon)


Mar_rain.fungi_ZOTU=data.frame(otu_table(Mar_rain.fungi))
#blanks need to be the first columns 
head(sample_data(Mar_rain.fungi_ZOTU))
R.Blanks=row.names(sample_data(subset_samples(Mar_rain.fungi, Sample_or_Control=="Control")))
length(R.Blanks)
#11

unique(sample_data(Mar_leaf.fung)$extract_kit)
#I also need to group by sample/extraction orders

#I also need to group by sample/ extraction orders
Group1 = row.names(sample_data(subset_samples(Mar_rain.fungi, extract_kit=="Water_1"&Sample_or_Control!="Control")))
length(Group1)
#47
Group2 = row.names(sample_data(subset_samples(Mar_rain.fungi, extract_kit=="Water_2"&Sample_or_Control!="Control")))
length(Group2)
#47
Group3 = row.names(sample_data(subset_samples(Mar_rain.fungi, extract_kit=="Water_3"&Sample_or_Control!="Control")))
length(Group3)
#19
#Group6 = row.names(sample_data(subset_samples(Mar_rain.fungi, extract_kit=="Mock")))
#length(Group6)
#3


Mar_rain.fungi_ZOTU$OTUs=row.names(Mar_rain.fungi_ZOTU)
colnames(Mar_rain.fungi_ZOTU)
Mar_rain.fungi_ZOTU_reorder=Mar_rain.fungi_ZOTU[,c("OTUs",R.Blanks,Group1,Group2,Group3)]#,Group6
colnames(Mar_rain.fungi_ZOTU_reorder)
nrow(Mar_rain.fungi_ZOTU_reorder)
#14921



Mar_rain.fungi_ZOTU_reorder_decon_obj=decon(Mar_rain.fungi_ZOTU_reorder,numb.blanks = 11,numb.ind=c(47,47,19), taxa=F)#,3
nrow(Mar_rain.fungi_ZOTU_reorder_decon_obj$decon.table)
#10420

R.decon_names=Mar_rain.fungi_ZOTU_reorder_decon_obj$OTUs.removed$OTUs
length(R.decon_names)
#4815



Mar_rain.fung_decon_table=Mar_rain.fungi_ZOTU_reorder_decon_obj$decon.table

head(Mar_rain.fung_decon_table)

#reformat for use in phyloseq

row.names(Mar_rain.fung_decon_table)=Mar_rain.fung_decon_table$OTUs
Mar_rain.fung_decon_table[,c("OTUs", "Mean.blank")]=NULL
head(Mar_rain.fung_decon_table)

Mar_rain.fung_decon=phyloseq(otu_table(Mar_rain.fung_decon_table, taxa_are_rows = T), sample_data(Mar_rain.fungi), tax_table(Mar_rain.fungi))
nsamples(Mar_rain.fung_decon)
#113
ntaxa(Mar_rain.fung_decon)
#10420

sum(otu_table(Mar_rain.fung_decon))
#2840603

max(sample_sums(Mar_rain.fung_decon))
#49118

min(sample_sums(Mar_rain.fung_decon))
#21

sort(sample_sums(Mar_rain.fung_decon))

#RX026GLBRC26  RX037GLBRC37  RX065GLBRC65  RX063GLBRC63 RX117GLBRC117  RX030GLBRC30  RX061GLBRC61  RX048GLBRC48  RX066GLBRC66  RX041GLBRC41  RX039GLBRC39  RX043GLBRC43  RX057GLBRC57 
#21           473           768          1402          7911          9449         11684         13234         13463         13845         14605         15121         15283


Mar_rain.fung_decon_pr=prune_samples(sample_sums(Mar_rain.fung_decon) > 1000, Mar_rain.fung_decon)
Mar_rain.fung_decon_pr=prune_taxa(taxa_sums(Mar_rain.fung_decon_pr) > 0, Mar_rain.fung_decon_pr)


nsamples(Mar_rain.fung_decon_pr)
#110
ntaxa(Mar_rain.fung_decon_pr)
#10024

sum(otu_table(Mar_rain.fung_decon_pr))
#2839341

max(sample_sums(Mar_rain.fung_decon_pr))
#49118

min(sample_sums(Mar_rain.fung_decon_pr))
#1402



#####Rain Contaminant removal####


#####Combine Leaf and Rain post decontamination####


unique(sample_data(Mar_leaf.fung_decon_pr)$sample_type)
unique(sample_data(Mar_rain.fung_decon_pr)$sample_type)

unique(sample_data(Mar_leaf.fung_decon_pr)$sampleID_bact)
unique(sample_data(Mar_rain.fung_decon_pr)$sampleID_bact)

Mar_leaf_rain.fung_decon=merge_phyloseq(Mar_leaf.fung_decon_pr,Mar_rain.fung_decon_pr)
nsamples(Mar_leaf_rain.fung_decon)
#171

ntaxa(Mar_leaf_rain.fung_decon)
#10990

sum(otu_table(Mar_leaf_rain.fung_decon))
#3522590

max(sample_sums(Mar_leaf_rain.fung_decon))
#49118

min(sample_sums(Mar_leaf_rain.fung_decon))
#1150

sort(sample_sums(Mar_leaf_rain.fung_decon))
#Leaf07        Leaf05        Leaf57        Leaf16  RX063GLBRC63        Leaf32        Leaf61        Leaf41        Leaf73        Leaf09        Leaf70        Leaf14        Leaf21 
#1150          1203          1209          1349          1402          1434          1623          1630          2009          2010          2344          2643          2712 
#Leaf69        Leaf56        Leaf13        Leaf45        Leaf25        Leaf19        Leaf34        Leaf30        Leaf71        Leaf66        Leaf22        Leaf03        Leaf27 
#2859          2962          2993          3038          3073          3104          3118          3240          4426          4546          4557          4803          5013 
#Leaf10        Leaf36        Leaf04        Leaf38        Leaf40        Leaf48        Leaf49        Leaf54 RX117GLBRC117        Leaf52        Leaf46        Leaf31  RX030GLBRC30 
#5318          5601          6703          6834          6886          7129          7626          7738          7911          8020          8068          9376          9449 
#Leaf18        Leaf44        Leaf60        Leaf28  RX061GLBRC61        Leaf42  RX048GLBRC48  RX066GLBRC66  RX041GLBRC41  RX039GLBRC39  RX043GLBRC43  RX057GLBRC57        Leaf15 
#9638          9834         10929         11092         11684         12841         13234         13463         13845         14605         15121         15283         15458





#####Leaf and Rain Rarefied Community Analyses####

###
Mar_leaf_rain.fung_decon_rar=rarefy_even_depth(Mar_leaf_rain.fung_decon, sample.size =1000, rngseed=101,replace=F)
#`set.seed(101)` was used to initialize repeatable random subsampling.
#Please record this for your records so others can reproduce.
#Try `set.seed(101); .Random.seed` for the full vector
#...
#6511OTUs were removed because they are no longer
#present in any sample after random subsampling

save(Mar_leaf_rain.fung_decon_rar, file = here::here("R_file","Mar_leaf_rain_ZOTU_full_v4.2.2020.fung_decon_rar_pruned_phyloseq_obj.RData"))
load(here::here("R_file","Mar_leaf_rain_ZOTU_full_v4.2.2020.fung_decon_rar_pruned_phyloseq_obj.RData"))

#I need to the rep set here for comparison with previous published sequences


rep_set.FunRainLeaf2019_ZOTU_rar1000=rep_set.FunRainLeaf2019_ZOTU_full[names(rep_set.FunRainLeaf2019_ZOTU_full) %in% taxa_names(Mar_leaf_rain.fung_decon_rar)]
head(rep_set.FunRainLeaf2019_ZOTU_rar1000)
length(rep_set.FunRainLeaf2019_ZOTU_rar1000)
#4479


write.fasta(sequences =rep_set.FunRainLeaf2019_ZOTU_rar1000, names = names(rep_set.FunRainLeaf2019_ZOTU_rar1000), 
            file.out =here::here("R_file","rep_set.FunRainLeaf2019_v4.2.2020_ZOTU_rar1000.fna"))

#####Bray Rarefied Pairwise community distance####
#I am only interested in the Marshall Proj Rain
head(sample_data(Mar_leaf_rain.fung_decon_rar))
Mar_leaf_rain.fung_decon_rar_Mar=subset_samples(Mar_leaf_rain.fung_decon_rar, sub_proj=="Marshall")
nsamples(Mar_leaf_rain.fung_decon_rar_Mar)
#117

# pairwise measure 
#now I need to make the betapart object
Mar_leaf_rain.fung_decon_rar_Mar_map=sample_data(Mar_leaf_rain.fung_decon_rar_Mar)
Mar_leaf_rain.fung_decon_rar_Mar_tbl=data.frame(t(otu_table(Mar_leaf_rain.fung_decon_rar_Mar)))
summary(Mar_leaf_rain.fung_decon_rar_Mar_tbl[,1:3])
Mar_leaf_rain.fung_decon_rar_Mar.core <- betapart.core.abund(Mar_leaf_rain.fung_decon_rar_Mar_tbl)

# multiple site measures
(Mar_leaf_rain.fung_decon_rar_Mar.multi <- beta.multi.abund(Mar_leaf_rain.fung_decon_rar_Mar.core))
#$beta.BRAY.BAL
#[1] 0.9819595

#$beta.BRAY.GRA
#[1] 0

#$beta.BRAY
#[1] 0.9819595

Mar_leaf_rain.fung_decon_rar_Mar.pair.s <- beta.pair.abund(Mar_leaf_rain.fung_decon_rar_Mar.core)
Mar_leaf_rain.fung_decon_rar_Mar.bray_M <- matrixConvert(Mar_leaf_rain.fung_decon_rar_Mar.pair.s$beta.bray, 
                                                         colname = c("sample1", "sample2", "bray"))#total distance 

Mar_leaf_rain.fung_decon_rar_Mar.bray_M$s1_s2=with(Mar_leaf_rain.fung_decon_rar_Mar.bray_M, interaction(sample1,sample2))
nrow(Mar_leaf_rain.fung_decon_rar_Mar.bray_M)
#6786
head(Mar_leaf_rain.fung_decon_rar_Mar.bray_M)
Mar_leaf_rain.fung_decon_rar_Mar.grad_M <- matrixConvert(Mar_leaf_rain.fung_decon_rar_Mar.pair.s$beta.bray.gra, 
                                                         colname = c("sample1", "sample2", "gradient"))#total distance 

Mar_leaf_rain.fung_decon_rar_Mar.grad_M$s1_s2=with(Mar_leaf_rain.fung_decon_rar_Mar.grad_M, interaction(sample1,sample2))
nrow(Mar_leaf_rain.fung_decon_rar_Mar.grad_M)
#6786
head(Mar_leaf_rain.fung_decon_rar_Mar.grad_M)
Mar_leaf_rain.fung_decon_rar_Mar.grad_M[,c("sample1", "sample2")]=NULL

Mar_leaf_rain.fung_decon_rar_Mar.bal_M <- matrixConvert(Mar_leaf_rain.fung_decon_rar_Mar.pair.s$beta.bray.bal, 
                                                        colname = c("sample1", "sample2", "balancing"))#total distance 

Mar_leaf_rain.fung_decon_rar_Mar.bal_M$s1_s2=with(Mar_leaf_rain.fung_decon_rar_Mar.bal_M, interaction(sample1,sample2))
nrow(Mar_leaf_rain.fung_decon_rar_Mar.bal_M)
#6786
head(Mar_leaf_rain.fung_decon_rar_Mar.bal_M)
Mar_leaf_rain.fung_decon_rar_Mar.bal_M[,c("sample1", "sample2")]=NULL

Mar_leaf_rain.fung_decon_rar_Mar.betapair0=merge(Mar_leaf_rain.fung_decon_rar_Mar.bray_M,Mar_leaf_rain.fung_decon_rar_Mar.grad_M, by="s1_s2")
Mar_leaf_rain.fung_decon_rar_Mar.betapair=merge(Mar_leaf_rain.fung_decon_rar_Mar.betapair0,Mar_leaf_rain.fung_decon_rar_Mar.bal_M, by="s1_s2")

Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt0=merge(Mar_leaf_rain.fung_decon_rar_Mar.betapair,Mar_leaf_rain.fung_decon_rar_Mar_map[,c("sample_type","collect_date","sub_proj",
                                                                                                                                       "plant_type","gh_block","rain_trt","UTM_lat","UTM_long")], 
                                                     by.x = "sample1",by.y = "row.names")
nrow(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt0)
#6786
head(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt0)
colnames(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt0)
colnames(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt0)[7:14]=c("s1_sample_type","s1_collect_date","s1_sub_proj","s1_plant_type",
                                                                 "s1_gh_block","s1_rain_trt","s1_UTM_lat","s1_UTM_long")
head(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt0)


Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt=merge(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt0,Mar_leaf_rain.fung_decon_rar_Mar_map[,c("sample_type","collect_date","sub_proj",
                                                                                                                                           "plant_type","gh_block","rain_trt","UTM_lat","UTM_long")], 
                                                    by.x = "sample2",by.y = "row.names")


nrow(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt)
#6786
ncol(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt)
head(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt)
colnames(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt)[15:22]=c("s2_sample_type","s2_collect_date","s2_sub_proj","s2_plant_type",
                                                                 "s2_gh_block","s2_rain_trt","s2_UTM_lat","s2_UTM_long")

summary(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt)

Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt$s1_s2_plant_type=with(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt,interaction(s1_plant_type,s2_plant_type))

Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt$s1_s2_collect_date=with(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt,interaction(s1_collect_date,s2_collect_date))

Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt$s1_collect_date_plant_type=with(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt,interaction(s1_collect_date,s1_plant_type))

Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt$s2_collect_date_plant_type=with(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt,interaction(s2_collect_date,s2_plant_type))

head(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt)


#I want to add in some raw numbers of species share and not shared

summary(Mar_leaf_rain.fung_decon_rar_Mar.core)


Mar_leaf_rain.fung_decon_rar_Mar_B.shared_M <- matrixConvert(Mar_leaf_rain.fung_decon_rar_Mar.core$pair.shared.abund, 
                                                             colname = c("sample1", "sample2", "shared_abun_tax"))#number of shared taxa
head(Mar_leaf_rain.fung_decon_rar_Mar_B.shared_M)
Mar_leaf_rain.fung_decon_rar_Mar_B.not_shared_M <- matrixConvert(Mar_leaf_rain.fung_decon_rar_Mar.core$pair.sum.not.shared.abund, 
                                                                 colname = c("sample1", "sample2", "not_shared_abun_tax"))#number of not shared taxa
head(Mar_leaf_rain.fung_decon_rar_Mar_B.not_shared_M)


Mar_leaf_rain.fung_decon_rar_Mar_B.shared_M$s1_s2=with(Mar_leaf_rain.fung_decon_rar_Mar_B.shared_M, interaction(sample1,sample2))
Mar_leaf_rain.fung_decon_rar_Mar_B.shared_M[,c("sample1", "sample2")]=NULL

Mar_leaf_rain.fung_decon_rar_Mar_B.not_shared_M$s1_s2=with(Mar_leaf_rain.fung_decon_rar_Mar_B.not_shared_M, interaction(sample1,sample2))
Mar_leaf_rain.fung_decon_rar_Mar_B.not_shared_M[,c("sample1", "sample2")]=NULL

Mar_leaf_rain.fung_decon_rar_Mar_B.shared_or_not_M=merge(Mar_leaf_rain.fung_decon_rar_Mar_B.shared_M,Mar_leaf_rain.fung_decon_rar_Mar_B.not_shared_M, by="s1_s2")
Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt=merge(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt,Mar_leaf_rain.fung_decon_rar_Mar_B.shared_or_not_M, by="s1_s2")
head(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt)
Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt$similarity_B=1-Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt$bray
Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt$s1_s2_rain_trt=with(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt,interaction(s1_rain_trt,s2_rain_trt))

Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt$s1_s2_gh_block=with(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt,
                                                                  interaction(s1_gh_block,s2_gh_block))
head(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt)
write.csv(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt,here::here("R_file","Pairwise_turnover_grad_distance_fung_v4.2.2020_decon_pruned_ZOTU_full_rar_Marshall.csv"))

bray_paiwise=Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt %>% group_by(s1_s2_plant_type,s1_s2_collect_date) %>% summarise_at("bray", 
                                                                                                                              list(~mean(.),se=~sd(.)/sqrt(n()),~n()))


#####Jaccard Rarefied Pairwise community distance####
# pairwise measure 
#now I need to make the betapart object
Mar_leaf_rain.fung_decon_rar_Mar_map=sample_data(Mar_leaf_rain.fung_decon_rar_Mar)

Mar_leaf_rain.fung_decon_rar_Mar_PA=transform_sample_counts(Mar_leaf_rain.fung_decon_rar_Mar,pa)
Mar_leaf_rain.fung_decon_rar_Mar_PA_tbl=data.frame(t(otu_table(Mar_leaf_rain.fung_decon_rar_Mar_PA)))
summary(Mar_leaf_rain.fung_decon_rar_Mar_PA_tbl[,1:3])
Mar_leaf_rain.fung_decon_rar_Mar_PA.core <- betapart.core(Mar_leaf_rain.fung_decon_rar_Mar_PA_tbl)

# multiple site measures
(Mar_leaf_rain.fung_decon_rar_Mar_PA.multi <- beta.multi(Mar_leaf_rain.fung_decon_rar_Mar_PA.core,index.family="jaccard"))
#$beta.JTU
#[1] 0.9800176

#$beta.JNE
#[1] 0.00963593

#$beta.JAC
#[1] 0.9896535

Mar_leaf_rain.fung_decon_rar_Mar_PA.pair.s <- beta.pair(Mar_leaf_rain.fung_decon_rar_Mar_PA.core,index.family="jaccard")
Mar_leaf_rain.fung_decon_rar_Mar_PA.jac_M <- matrixConvert(Mar_leaf_rain.fung_decon_rar_Mar_PA.pair.s$beta.jac, 
                                                           colname = c("sample1", "sample2", "jaccard"))#total distance 

Mar_leaf_rain.fung_decon_rar_Mar_PA.jac_M$s1_s2=with(Mar_leaf_rain.fung_decon_rar_Mar_PA.jac_M, interaction(sample1,sample2))
nrow(Mar_leaf_rain.fung_decon_rar_Mar_PA.jac_M)
#6786
head(Mar_leaf_rain.fung_decon_rar_Mar_PA.jac_M)
Mar_leaf_rain.fung_decon_rar_Mar_PA.net_M <- matrixConvert(Mar_leaf_rain.fung_decon_rar_Mar_PA.pair.s$beta.jne, 
                                                           colname = c("sample1", "sample2", "nestedness"))#total distance 

Mar_leaf_rain.fung_decon_rar_Mar_PA.net_M$s1_s2=with(Mar_leaf_rain.fung_decon_rar_Mar_PA.net_M, interaction(sample1,sample2))
nrow(Mar_leaf_rain.fung_decon_rar_Mar_PA.net_M)
#6786
head(Mar_leaf_rain.fung_decon_rar_Mar_PA.net_M)
Mar_leaf_rain.fung_decon_rar_Mar_PA.net_M[,c("sample1","sample2")]=NULL


Mar_leaf_rain.fung_decon_rar_Mar_PA.tur_M <- matrixConvert(Mar_leaf_rain.fung_decon_rar_Mar_PA.pair.s$beta.jtu, 
                                                           colname = c("sample1", "sample2", "turnover"))#total distance 

Mar_leaf_rain.fung_decon_rar_Mar_PA.tur_M$s1_s2=with(Mar_leaf_rain.fung_decon_rar_Mar_PA.tur_M, interaction(sample1,sample2))
nrow(Mar_leaf_rain.fung_decon_rar_Mar_PA.tur_M)
#6786
head(Mar_leaf_rain.fung_decon_rar_Mar_PA.tur_M)
Mar_leaf_rain.fung_decon_rar_Mar_PA.tur_M[,c("sample1","sample2")]=NULL
Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair0=merge(Mar_leaf_rain.fung_decon_rar_Mar_PA.jac_M,Mar_leaf_rain.fung_decon_rar_Mar_PA.net_M, by="s1_s2")
Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair=merge(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair0,Mar_leaf_rain.fung_decon_rar_Mar_PA.tur_M, by="s1_s2")

Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt0=merge(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair,Mar_leaf_rain.fung_decon_rar_Mar_map[,c("sample_type","collect_date","sub_proj",
                                                                                                                                             "plant_type","gh_block","rain_trt","UTM_lat","UTM_long")], 
                                                        by.x = "sample1",by.y = "row.names")
nrow(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt0)
#6786
head(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt0)
colnames(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt0)[7:14]=c("s1_sample_type","s1_collect_date","s1_sub_proj","s1_plant_type",
                                                                    "s1_gh_block","s1_rain_trt","s1_UTM_lat","s1_UTM_long")
head(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt0)


Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt=merge(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt0,Mar_leaf_rain.fung_decon_rar_Mar_map[,c("sample_type","collect_date","sub_proj",
                                                                                                                                                 "plant_type","gh_block","rain_trt","UTM_lat","UTM_long")], 
                                                       by.x = "sample2",by.y = "row.names")


nrow(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt)
#6786
ncol(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt)
head(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt)
colnames(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt)[15:22]=c("s2_sample_type","s2_collect_date","s2_sub_proj","s2_plant_type",
                                                                    "s2_gh_block","s2_rain_trt","s2_UTM_lat","s2_UTM_long")

summary(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt)

Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt$s1_s2_plant_type=with(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt,interaction(s1_plant_type,s2_plant_type))

Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt$s1_s2_collect_date=with(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt,interaction(s1_collect_date,s2_collect_date))

Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt$s1_collect_date_plant_type=with(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt,interaction(s1_collect_date,s1_plant_type))

Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt$s2_collect_date_plant_type=with(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt,interaction(s2_collect_date,s2_plant_type))

head(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt)


#I want to add in some raw numbers of species share and not shared

summary(Mar_leaf_rain.fung_decon_rar_Mar_PA.core)


Mar_leaf_rain.fung_decon_rar_Mar_PA_B.shared_M <- matrixConvert(Mar_leaf_rain.fung_decon_rar_Mar_PA.core$shared, 
                                                                colname = c("sample1", "sample2", "shared"))#number of shared taxa
head(Mar_leaf_rain.fung_decon_rar_Mar_PA_B.shared_M)
Mar_leaf_rain.fung_decon_rar_Mar_PA_B.not_shared_M <- matrixConvert(Mar_leaf_rain.fung_decon_rar_Mar_PA.core$not.shared, 
                                                                    colname = c("sample1", "sample2", "not.shared"))#number of not shared taxa
head(Mar_leaf_rain.fung_decon_rar_Mar_PA_B.not_shared_M)


Mar_leaf_rain.fung_decon_rar_Mar_PA_B.shared_M$s1_s2=with(Mar_leaf_rain.fung_decon_rar_Mar_PA_B.shared_M, interaction(sample1,sample2))
Mar_leaf_rain.fung_decon_rar_Mar_PA_B.shared_M[,c("sample1", "sample2")]=NULL

Mar_leaf_rain.fung_decon_rar_Mar_PA_B.not_shared_M$s1_s2=with(Mar_leaf_rain.fung_decon_rar_Mar_PA_B.not_shared_M, interaction(sample1,sample2))
Mar_leaf_rain.fung_decon_rar_Mar_PA_B.not_shared_M[,c("sample1", "sample2")]=NULL

Mar_leaf_rain.fung_decon_rar_Mar_PA_B.shared_or_not_M=merge(Mar_leaf_rain.fung_decon_rar_Mar_PA_B.shared_M,Mar_leaf_rain.fung_decon_rar_Mar_PA_B.not_shared_M, by="s1_s2")
Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt=merge(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt,Mar_leaf_rain.fung_decon_rar_Mar_PA_B.shared_or_not_M, by="s1_s2")
head(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt)
Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt$similarity=1-Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt$jaccard
Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt$s1_s2_rain_trt=with(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt,interaction(s1_rain_trt,s2_rain_trt))

Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt$s1_s2_gh_block=with(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt,
                                                                     interaction(s1_gh_block,s2_gh_block))
head(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt)
write.csv(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt,here::here("R_file","Pairwise_turnover_Jaccard_fung_v4.2.2020_decon_pruned_ZOTU_full_rar_Marshall.csv"))





#####Begin analyses#####



#####Petri Bray analyses No Nano####

load(here::here("R_file","Mar_leaf_rain_ZOTU_full_v4.2.2020.fung_decon_rar_pruned_phyloseq_obj.RData"))
#I am only interested in the Marshall Proj Rain
head(sample_data(Mar_leaf_rain.fung_decon_rar))
Mar_leaf_rain.fung_decon_rar_Mar=subset_samples(Mar_leaf_rain.fung_decon_rar, sub_proj=="Marshall")
nsamples(Mar_leaf_rain.fung_decon_rar_Mar)
#117

#Let's look at Rain and the petri experiment
head(sample_data(Mar_leaf_rain.fung_decon_rar))
Mar_leaf_rain.fung_decon_rar_Mar_petri=subset_samples(Mar_leaf_rain.fung_decon_rar_Mar, plant_type!="Adult"&plant_type!="Seedling")
nsamples(Mar_leaf_rain.fung_decon_rar_Mar_petri)
#92
Mar_leaf.fung_decon_pr_petr_S.rar=subset_samples(Mar_leaf_rain.fung_decon_rar_Mar_petri,plant_type!="Rain"|collect_date=="8/21/2018"|
                                                     collect_date=="7/21/2018")
nsamples(Mar_leaf.fung_decon_pr_petr_S.rar)
#43

Mar_leaf.fung_decon_pr_petr_S_nan.rar=subset_samples(Mar_leaf.fung_decon_pr_petr_S.rar,rain_trt!="Nano")
nsamples(Mar_leaf.fung_decon_pr_petr_S_nan.rar)
#35


Mar_leaf.fung_decon_pr_petr_S_nan.rar=prune_taxa(taxa_sums(Mar_leaf.fung_decon_pr_petr_S_nan.rar) > 0, Mar_leaf.fung_decon_pr_petr_S_nan.rar)
ntaxa(Mar_leaf.fung_decon_pr_petr_S_nan.rar)
#944
sum(taxa_sums(Mar_leaf.fung_decon_pr_petr_S_nan.rar))
#35000
min(sample_sums(Mar_leaf.fung_decon_pr_petr_S_nan.rar))
#1000
max(sample_sums(Mar_leaf.fung_decon_pr_petr_S_nan.rar))
#1000



Mar_leaf.fung_decon_pr_petr_S_nan.rar.ord=ordinate(Mar_leaf.fung_decon_pr_petr_S_nan.rar,method = "NMDS",distance = "bray")
#*** Solution reached
#0.1793605       


#Creation of files for PERMANOVA analyses in Primer v6

Mar_leaf.fung_decon_pr_petr_S_nan.rar_dis=phyloseq::distance(Mar_leaf.fung_decon_pr_petr_S_nan.rar,method = "bray")
Mar_leaf.fung_decon_pr_petr_S_nan.rar_map=sample_data(Mar_leaf.fung_decon_pr_petr_S_nan.rar)

write.csv(as.matrix(Mar_leaf.fung_decon_pr_petr_S_nan.rar_dis),here::here("R_file","Bray_Mar_leaf.fung_decon_pr_petr_S_nan.rar_dis.csv"))
write.csv(as.matrix(Mar_leaf.fung_decon_pr_petr_S_nan.rar_map),here::here("R_file","Treatments_Mar_leaf.fung_decon_pr_petr_S_nan.rar_map.csv"))


####Petri Bray Betadispersion#####

Mar_leaf.fung_decon_pr_petr_S_nan.rar_betamod=betadisper(Mar_leaf.fung_decon_pr_petr_S_nan.rar_dis, 
                                         with(Mar_leaf.fung_decon_pr_petr_S_nan.rar_map,interaction(plant_type,collect_date,rain_trt)))

Mar_leaf.fung_decon_pr_petr_S_nan.rar_betamod_dist=as.data.frame(Mar_leaf.fung_decon_pr_petr_S_nan.rar_betamod$distances)
colnames(Mar_leaf.fung_decon_pr_petr_S_nan.rar_betamod_dist)="betadisp"

Mar_leaf.fung_decon_pr_petr_S_nan.rar_betamod_dist=merge(Mar_leaf.fung_decon_pr_petr_S_nan.rar_betamod_dist,Mar_leaf.fung_decon_pr_petr_S_nan.rar_map,by="row.names")


#####Graphing the Bray Betadispersion of Petri No nano####

brewer.pal(n = 12, name = "Paired")
head(Mar_leaf.fung_decon_pr_petr_S_nan.rar_betamod_dist)


Mar_leaf.fung_decon_pr_petr_S_nan.rar_betamod_dist$round=ifelse(Mar_leaf.fung_decon_pr_petr_S_nan.rar_betamod_dist$rain_trt=="No_rain","Seed",ifelse(Mar_leaf.fung_decon_pr_petr_S_nan.rar_betamod_dist$collect_date=="9/18/2018"|Mar_leaf.fung_decon_pr_petr_S_nan.rar_betamod_dist$collect_date=="8/21/2018","Round2","Round1"))

Mar_leaf.fung_decon_pr_petr_S_nan.rar_betamod_dist_sum=Mar_leaf.fung_decon_pr_petr_S_nan.rar_betamod_dist%>%group_by(round,collect_date,rain_trt)%>%
  summarise_at(vars("betadisp"),list(~mean(.),~n(),se=~sd(.)/sqrt(n())))
unique(with(Mar_leaf.fung_decon_pr_petr_S_nan.rar_betamod_dist,interaction(collect_date,rain_trt)))
date_trt=c("7/10/2018.No_rain","8/14/2018.Sterile_rain","8/14/2018.Live_Rain","7/21/2018.Rain",
           "9/18/2018.Sterile_rain","9/18/2018.Live_Rain","8/21/2018.Rain")



(petri_betadisp_nan_p=ggplot(Mar_leaf.fung_decon_pr_petr_S_nan.rar_betamod_dist)+
    geom_point(aes(x=factor(rain_trt, levels = c("No_rain","Sterile_rain","Live_Rain","Rain"),
                            labels = c("Seed","Sterile","Live","Rain")),
                   y=betadisp,shape=rain_trt,fill=factor(rain_trt, levels = c("No_rain","Sterile_rain","Live_Rain","Rain")),
                   color=factor(rain_trt, levels = c("No_rain","Sterile_rain","Live_Rain","Rain"))),size=3,alpha=0.5)+
    geom_errorbar(data=Mar_leaf.fung_decon_pr_petr_S_nan.rar_betamod_dist_sum,aes(x=factor(rain_trt, levels = c("No_rain","Sterile_rain","Live_Rain","Rain"),
                                                                       labels = c("Seed","Sterile","Live","Rain")),
                                                              y=mean,ymax=mean+se,ymin=mean-se),
                  width=0.5)+facet_grid(~factor(round, levels= c("Seed","Round1","Round2"),
                                                labels = c("","Round1","Round2")),scale = "free_x",space = "free_x")+
    geom_point(data=Mar_leaf.fung_decon_pr_petr_S_nan.rar_betamod_dist_sum,aes(x=factor(rain_trt, levels = c("No_rain","Sterile_rain","Live_Rain","Rain"),
                                                                    labels = c("Seed","Sterile","Live","Rain")),y=mean,shape=rain_trt,
                                                           fill=factor(rain_trt, levels = c("No_rain","Sterile_rain","Live_Rain","Rain")),
                                                           color=factor(rain_trt, levels = c("No_rain","Sterile_rain","Live_Rain","Rain"))),
               size=10, stroke=1.5)+
    scale_shape_manual(values = c(24,13,25,21))+ylab("Beta-dispersion (Bray-Curtis)")+
    scale_color_manual(values = c("#E31A1C","black", "black","black"),name=NULL)+
    scale_fill_manual(values = c("#E31A1C","#6A3D9A","#1F78B4","black"),name=NULL)+
    scale_x_discrete(name=NULL)+theme_cowplot()+
    theme(axis.title = element_text(size = 32),axis.text = element_text(size = 28),legend.title = element_blank(),legend.position = "none",
          strip.background = element_rect(fill = NA),strip.text = element_text(size = 32)))

#1200x800
#geom_line(data=Mar_leaf_rar_petri_div_map_nan_sum,aes(x=factor(interaction(collect_date,rain_trt),levels = date_trt),group=collect_date,
#y=Observed_mean),size=2)+

Mar_leaf.fung_decon_pr_petr_S_nan.rar_betamod_dist$trt_date=with(Mar_leaf.fung_decon_pr_petr_S_nan.rar_betamod_dist,interaction(rain_trt,collect_date))

#Bray Betadispersion Analyses

Petri_betadisp_mod_nan=lmer((betadisp)^(-1)~trt_date+(1|gh_block),Mar_leaf.fung_decon_pr_petr_S_nan.rar_betamod_dist)
plot(Petri_betadisp_mod_nan)
hist(resid(Petri_betadisp_mod_nan))
qqPlot(resid(Petri_betadisp_mod_nan))
shapiro.test(resid(Petri_betadisp_mod_nan))
#W = 0.96579, p-value = 0.3384

anova(Petri_betadisp_mod_nan)
#trt_date 8.7538   1.459     6    28   2.098 0.08538 .


emmeans(Petri_betadisp_mod_nan, pairwise~trt_date,adjust="fdr")


ggplot(Mar_leaf.fung_decon_pr_petr_S_nan.rar_betamod_dist,aes(x=plant_type,y=betadisp))+geom_point()+geom_boxplot()




#####LFE Analyses the betadisp of Petri No nano####


#Petri 

Mar_leaf.fung_decon_pr_petr_S_nan.rar_betamod_dist_p=subset(Mar_leaf.fung_decon_pr_petr_S_nan.rar_betamod_dist,plant_type=="Petri")

#Let's make a round factor
unique(Mar_leaf.fung_decon_pr_petr_S_nan.rar_betamod_dist_p$collect_date)
Mar_leaf.fung_decon_pr_petr_S_nan.rar_betamod_dist_p$round=as.factor(ifelse(Mar_leaf.fung_decon_pr_petr_S_nan.rar_betamod_dist_p$collect_date=="9/18/2018"|Mar_leaf.fung_decon_pr_petr_S_nan.rar_betamod_dist_p$collect_date=="8/21/2018",
                                                        "Round2","Round1"))
summary(Mar_leaf.fung_decon_pr_petr_S_nan.rar_betamod_dist_p)
#Observed Richness

Petri_LFE_betadisp_mod=lmer((betadisp)^-1~rain_trt*round+(1|gh_block),data=Mar_leaf.fung_decon_pr_petr_S_nan.rar_betamod_dist_p)
plot(Petri_LFE_betadisp_mod)
hist(resid(Petri_LFE_betadisp_mod))
qqPlot(resid(Petri_LFE_betadisp_mod))
shapiro.test(resid(Petri_LFE_betadisp_mod))
#W = 0.93014, p-value = 0.09817

####Table S5 Beta-dispersion (Bray)####
anova(Petri_LFE_betadisp_mod)
#rain_trt       0.063509 0.063509     1 17.308  0.1043 0.7506
#round          0.024095 0.024095     1 19.451  0.0396 0.8444
#rain_trt:round 0.078410 0.078410     1 17.308  0.1287 0.7241


emmeans(Petri_LFE_betadisp_mod, pairwise~rain_trt|round)


ggplot(Mar_leaf.fung_decon_pr_petr_S_nan.rar_betamod_dist_p,aes(x=interaction(rain_trt,start_date),y=betadisp))+geom_point()+geom_boxplot()


#####Petri NMDS graph####

#Hull for the Bray NMDS
Mar_leaf.fung_decon_pr_petr_S_nan.rar_map=sample_data(Mar_leaf.fung_decon_pr_petr_S_nan.rar)
summary(Mar_leaf.fung_decon_pr_petr_S_nan.rar_map)
Mar_leaf.fung_decon_pr_petr_S_nan.rar_coordinates=merge(Mar_leaf.fung_decon_pr_petr_S_nan.rar.ord$points,Mar_leaf.fung_decon_pr_petr_S_nan.rar_map, by="row.names")
Mar_leaf.fung_decon_pr_petr_S_nan.rar_coordinates$plant_data_grp=ifelse(Mar_leaf.fung_decon_pr_petr_S_nan.rar_coordinates$plant_type=="Rain"&
                                                                      Mar_leaf.fung_decon_pr_petr_S_nan.rar_coordinates$collect_date=="7/21/2018",
                                                                           "Round1_Rain",ifelse(Mar_leaf.fung_decon_pr_petr_S_nan.rar_coordinates$plant_type=="Rain"&
                                                                                                  Mar_leaf.fung_decon_pr_petr_S_nan.rar_coordinates$collect_date=="8/21/2018",
                                                                                                "Round2_Rain", 
                                                                                          ifelse(Mar_leaf.fung_decon_pr_petr_S_nan.rar_coordinates$rain_trt=="Live_Rain",
                                                                                                 ifelse(Mar_leaf.fung_decon_pr_petr_S_nan.rar_coordinates$collect_date=="9/18/2018","Round2_Live","Round1_Live"),
                                                                                                 ifelse(Mar_leaf.fung_decon_pr_petr_S_nan.rar_coordinates$rain_trt=="Sterile_rain",
                                                                                                        ifelse(Mar_leaf.fung_decon_pr_petr_S_nan.rar_coordinates$collect_date=="9/18/2018","Round2_Sterile","Round1_Sterile"),"Seed"))))


unique(Mar_leaf.fung_decon_pr_petr_S_nan.rar_coordinates$plant_data_grp)

grp1000_Round1_Live <- Mar_leaf.fung_decon_pr_petr_S_nan.rar_coordinates[Mar_leaf.fung_decon_pr_petr_S_nan.rar_coordinates$plant_data_grp == "Round1_Live", 
                                                                            ][chull(Mar_leaf.fung_decon_pr_petr_S_nan.rar_coordinates[
                                                                              Mar_leaf.fung_decon_pr_petr_S_nan.rar_coordinates$
                                                                                plant_data_grp =="Round1_Live", c("MDS1", "MDS2")]), ] 
nrow(grp1000_Round1_Live)
#5
grp1000_Round2_Live <- Mar_leaf.fung_decon_pr_petr_S_nan.rar_coordinates[Mar_leaf.fung_decon_pr_petr_S_nan.rar_coordinates$plant_data_grp == "Round2_Live", 
                                                                             ][chull(Mar_leaf.fung_decon_pr_petr_S_nan.rar_coordinates[
                                                                               Mar_leaf.fung_decon_pr_petr_S_nan.rar_coordinates$
                                                                                 plant_data_grp =="Round2_Live", c("MDS1", "MDS2")]), ]  

nrow(grp1000_Round2_Live)
#4


grp1000_Round1_Sterile <- Mar_leaf.fung_decon_pr_petr_S_nan.rar_coordinates[Mar_leaf.fung_decon_pr_petr_S_nan.rar_coordinates$plant_data_grp == "Round1_Sterile", 
                                                                          ][chull(Mar_leaf.fung_decon_pr_petr_S_nan.rar_coordinates[
                                                                            Mar_leaf.fung_decon_pr_petr_S_nan.rar_coordinates$
                                                                              plant_data_grp =="Round1_Sterile", c("MDS1", "MDS2")]), ]  

nrow(grp1000_Round1_Sterile)
#5

grp1000_Round2_Sterile <- Mar_leaf.fung_decon_pr_petr_S_nan.rar_coordinates[Mar_leaf.fung_decon_pr_petr_S_nan.rar_coordinates$plant_data_grp == "Round2_Sterile", 
                                                                             ][chull(Mar_leaf.fung_decon_pr_petr_S_nan.rar_coordinates[
                                                                               Mar_leaf.fung_decon_pr_petr_S_nan.rar_coordinates$
                                                                                 plant_data_grp =="Round2_Sterile", c("MDS1", "MDS2")]), ]  

nrow(grp1000_Round2_Sterile)
#5


grp1000_P_Seed <- Mar_leaf.fung_decon_pr_petr_S_nan.rar_coordinates[Mar_leaf.fung_decon_pr_petr_S_nan.rar_coordinates$plant_data_grp == "Seed", 
                                                                          ][chull(Mar_leaf.fung_decon_pr_petr_S_nan.rar_coordinates[
                                                                            Mar_leaf.fung_decon_pr_petr_S_nan.rar_coordinates$
                                                                              plant_data_grp =="Seed", c("MDS1", "MDS2")]), ]  

nrow(grp1000_P_Seed)
#4


grp1000_P_Round1_Rain <- Mar_leaf.fung_decon_pr_petr_S_nan.rar_coordinates[Mar_leaf.fung_decon_pr_petr_S_nan.rar_coordinates$plant_data_grp == "Round1_Rain", 
                                                                     ][chull(Mar_leaf.fung_decon_pr_petr_S_nan.rar_coordinates[
                                                                       Mar_leaf.fung_decon_pr_petr_S_nan.rar_coordinates$
                                                                         plant_data_grp =="Round1_Rain", c("MDS1", "MDS2")]), ]  

nrow(grp1000_P_Round1_Rain)
#3

grp1000_P_Round2_Rain <- Mar_leaf.fung_decon_pr_petr_S_nan.rar_coordinates[Mar_leaf.fung_decon_pr_petr_S_nan.rar_coordinates$plant_data_grp == "Round2_Rain", 
                                                                       ][chull(Mar_leaf.fung_decon_pr_petr_S_nan.rar_coordinates[
                                                                         Mar_leaf.fung_decon_pr_petr_S_nan.rar_coordinates$
                                                                           plant_data_grp =="Round2_Rain", c("MDS1", "MDS2")]), ]  

nrow(grp1000_P_Round2_Rain)
#4

hull1000_Petri=rbind(grp1000_Round1_Live,grp1000_Round2_Live,grp1000_Round1_Sterile,grp1000_Round2_Sterile,
                     grp1000_P_Round1_Rain, grp1000_P_Round2_Rain,grp1000_P_Seed)
nrow(hull1000_Petri)
#30

brewer.pal(n = 12, name = "Paired")

#In order to use facet I need to duplicate seed
hull1000_Petri_seed=subset(hull1000_Petri,plant_type=="Seed")
Mar_leaf.fung_decon_pr_petr_S_nan.rar_coordinates_seed=subset(Mar_leaf.fung_decon_pr_petr_S_nan.rar_coordinates,plant_type=="Seed")

Mar_leaf.fung_decon_pr_petr_S_nan.rar_coordinates$round=ifelse(Mar_leaf.fung_decon_pr_petr_S_nan.rar_coordinates$collect_date=="9/18/2018"|
                                                                 Mar_leaf.fung_decon_pr_petr_S_nan.rar_coordinates$collect_date=="8/21/2018",
                            "Round2","Round1")
hull1000_Petri$round=ifelse(hull1000_Petri$collect_date=="9/18/2018"|hull1000_Petri$collect_date=="8/21/2018",
                            "Round2","Round1")

Mar_leaf.fung_decon_pr_petr_S_nan.rar_coordinates_seed$round=rep("Round2")

hull1000_Petri_seed$round=rep("Round2")


Mar_leaf.fung_decon_pr_petr_S_nan.rar_coordinates2=rbind(Mar_leaf.fung_decon_pr_petr_S_nan.rar_coordinates_seed,
                                                              Mar_leaf.fung_decon_pr_petr_S_nan.rar_coordinates)

hull1000_Petri2=rbind(hull1000_Petri,hull1000_Petri_seed)




unique(Mar_leaf.fung_decon_pr_petr_S_nan.rar_coordinates$plant_data_grp)

P_plant_data_grp_order=c("Seed","Round1_Sterile","Round2_Sterile","Round1_Live","Round2_Live","Round1_Rain","Round2_Rain")
length(P_plant_data_grp_order)
unique(Mar_leaf.fung_decon_pr_petr_S_nan.rar_coordinates$rain_trt)

(petri_nmds_bray=ggplot(Mar_leaf.fung_decon_pr_petr_S_nan.rar_coordinates2,aes(x= MDS1, y=MDS2))+
  geom_polygon(data=hull1000_Petri2,aes(x=MDS1,y=MDS2,fill=factor(plant_data_grp,levels = P_plant_data_grp_order),
                                         group=factor(plant_data_grp,levels = P_plant_data_grp_order)), alpha=0.5)+
  geom_point(size=4,aes(fill=factor(plant_data_grp,levels = P_plant_data_grp_order),
                        shape=rain_trt, color=factor(plant_data_grp,levels = P_plant_data_grp_order)),stroke=2)+
  scale_shape_manual(values = c(24,13,21,25),name=NULL)+
  scale_color_manual(values = c("#E31A1C","black", "black","black","black","black","black","black","black"),name=NULL)+
  scale_fill_manual(values = c("#E31A1C","#6A3D9A","#6A3D9A","#1F78B4","#1F78B4","black","black"),name=NULL)+
  facet_nested(.~round)+theme_classic()+theme(axis.title.x = element_blank(),axis.title.y = element_text(size = 28),
                                              axis.text = element_text(size = 30),panel.border =  element_rect(fill = NA,linetype="solid"),
                                              legend.position = "none",strip.text = element_text(size = 36),
                                              strip.placement = "outside",strip.background = element_rect(linetype="blank"),
                                              panel.spacing = unit(0, "lines")))

##


#####Petri Stacked Taxon Bargraph####



sample_data(Mar_leaf.fung_decon_pr_petr_S_nan.rar)$date_trt=droplevels(with(sample_data(Mar_leaf.fung_decon_pr_petr_S_nan.rar), interaction(collect_date, rain_trt)))


unique(sample_data(Mar_leaf.fung_decon_pr_petr_S_nan.rar)$date_trt)

ntaxa(Mar_leaf.fung_decon_pr_petr_S_nan.rar)
#944

#merge OTUs by the soil and precipitation treatment
Mar_leaf.fung_decon_pr_petr_S_nan.rar_fact=merge_samples(Mar_leaf.fung_decon_pr_petr_S_nan.rar, "date_trt")
sample_names(Mar_leaf.fung_decon_pr_petr_S_nan.rar_fact)     

#combine the reads at Phylum level
get_taxa_unique(Mar_leaf.fung_decon_pr_petr_S_nan.rar_fact, taxonomic.rank="Phylum")
#"Basidiomycota" "Unknown"       "Ascomycota"    "Mucoromycota"  "Glomeromycota"
(Mar_leaf.fung_decon_pr_petr_S_nan.rar_fact.phylum<-tax_glom(Mar_leaf.fung_decon_pr_petr_S_nan.rar_fact, taxrank="Phylum"))
#5


#subset so there is only the top ten most abundant phyla

Mar_leaf.fung_decon_pr_petr_S_nan.rar_Names=(data.frame(tax_table(Mar_leaf.fung_decon_pr_petr_S_nan.rar_fact.phylum))$Phylum)

#Transform the read counts to prop of total reads

Mar_leaf.fung_decon_pr_petr_S_nan.rar_fact.phylum.prop=transform_sample_counts(Mar_leaf.fung_decon_pr_petr_S_nan.rar_fact.phylum, function(x)x/sum(x))

petri_positions=c("7/10/2018.No_rain","8/14/2018.Sterile_rain","8/14/2018.Live_Rain","7/21/2018.Rain","9/18/2018.Sterile_rain",
                  "9/18/2018.Live_Rain","8/21/2018.Rain")


Mar_leaf.fung_decon_pr_petr_S_nan.rar_fact.phylum.prop_otu=as.data.frame(t(otu_table(Mar_leaf.fung_decon_pr_petr_S_nan.rar_fact.phylum.prop)))
Mar_leaf.fung_decon_pr_petr_S_nan.rar_fact.phylum.prop_otu[,"Phylum"]=Mar_leaf.fung_decon_pr_petr_S_nan.rar_Names
Mar_leaf.fung_decon_pr_petr_S_nan.rar_fact.phylum.prop_otu_M=melt(Mar_leaf.fung_decon_pr_petr_S_nan.rar_fact.phylum.prop_otu,id="Phylum")




Mar_leaf.fung_decon_pr_petr_S_nan.rar_fact.phylum.prop_otu_M2=Mar_leaf.fung_decon_pr_petr_S_nan.rar_fact.phylum.prop_otu_M%>%separate(variable,c("collectDate","rain_trt"),
                                                                                                                                       remove=F,sep = "[.]")

Mar_leaf.fung_decon_pr_petr_S_nan.rar_fact.phylum.prop_otu_M2$round=ifelse(Mar_leaf.fung_decon_pr_petr_S_nan.rar_fact.phylum.prop_otu_M2$collectDate=="7/10/2018","Seed",
                                                                           ifelse(Mar_leaf.fung_decon_pr_petr_S_nan.rar_fact.phylum.prop_otu_M2$collectDate=="8/14/2018"|
                                                                                  Mar_leaf.fung_decon_pr_petr_S_nan.rar_fact.phylum.prop_otu_M2$collectDate=="7/21/2018",
                                                                                  "Round1","Round2"))

#####PLOT: Figure S7####
(p_petri_color_F=ggplot(Mar_leaf.fung_decon_pr_petr_S_nan.rar_fact.phylum.prop_otu_M2,aes(x=factor(rain_trt,levels = c("No_rain","Sterile_rain","Live_Rain","Rain"),
                                                                                                   labels=c("Seeds","Sterile","Live","Rain")),
                                                                                          y=value,fill=Phylum))+
    geom_bar(aes( fill=Phylum), stat="identity", position="stack",color="black")+theme_cowplot()+
    theme(axis.text.y=element_text(size=26),axis.text.x=element_text(size=26),legend.position = "none",
          axis.title=element_text(size=32),panel.grid.major=element_blank(),legend.text = element_text(size=16), legend.title = element_text(size=20),
          strip.background = element_rect(fill = NA),strip.text = element_text(size=28))+xlab(NULL)+ylab("Proportion")+
    facet_grid(cols = vars(factor(round,levels = c("Seed","Round1","Round2"),labels=c("","Round1","Round2"))),space = "free",scales = "free")+
    guides(fill=guide_legend(title="Phyla"))+scale_fill_brewer(palette="Paired"))


(p_petri_color_F_leg=ggplot(Mar_leaf.fung_decon_pr_petr_S_nan.rar_fact.phylum.prop_otu_M2,aes(x=factor(rain_trt,levels = c("No_rain","Sterile_rain","Live_Rain","Rain"),
                                                                                                   labels=c("Seeds","Sterile","Live","Rain")),
                                                                                          y=value,fill=Phylum))+
    geom_bar(aes( fill=Phylum), stat="identity", position="stack",color="black")+theme_cowplot()+
    theme(axis.text.y=element_text(size=26),axis.text.x=element_text(size=26),
          axis.title=element_text(size=32),panel.grid.major=element_blank(),legend.text = element_text(size=16), legend.title = element_text(size=20),
          strip.background = element_rect(fill = NA),strip.text = element_text(size=28))+xlab(NULL)+ylab("Proportion")+
    facet_grid(cols = vars(factor(round,levels = c("Seed","Round1","Round2"),labels=c("","Round1","Round2"))),space = "free",scales = "free")+
    guides(fill=guide_legend(title="Phyla"))+scale_fill_brewer(palette="Paired"))

#
#class Level
#combine the reads at Class level
get_taxa_unique(Mar_leaf.fung_decon_pr_petr_S_nan.rar_fact, taxonomic.rank="Class")
#[1] "Tremellomycetes"      "Unknown"              "Dothideomycetes"      "Taphrinomycetes"      "Sordariomycetes"      "Exobasidiomycetes"    "Microbotryomycetes"  
#[8] "Cystobasidiomycetes"  "Agaricomycetes"       "Ustilaginomycetes"    "Saccharomycetes"      "Lecanoromycetes"      "Eurotiomycetes"       "Leotiomycetes"       
#[15] "Wallemiomycetes"      "Agaricostilbomycetes" "Spiculogloeomycetes"  "Orbiliomycetes"       "Classiculomycetes"    "Mucoromycetes"        "Glomeromycetes"      
#[22] "Laboulbeniomycetes"
(Mar_leaf.fung_decon_pr_petr_S_nan.rar.class<-tax_glom(Mar_leaf.fung_decon_pr_petr_S_nan.rar_fact, taxrank="Class"))
#24 taxa

#How many 

Mar_leaf.fung_decon_pr_petr_S_nan.rar.class_Names_raw=data.frame("Phylum"=data.frame(tax_table(Mar_leaf.fung_decon_pr_petr_S_nan.rar.class))$Phylum,
                                                                     "Class"=data.frame(tax_table(Mar_leaf.fung_decon_pr_petr_S_nan.rar.class))$Class,
                                                                     "OTU"=taxa_names(Mar_leaf.fung_decon_pr_petr_S_nan.rar.class))

Mar_leaf.fung_decon_pr_petr_S_nan.rar.class_Names_raw$P_C=with(Mar_leaf.fung_decon_pr_petr_S_nan.rar.class_Names_raw,interaction(Phylum,Class))
#Let's take the top 10 taxa 

P_TopCLASS_Field = names(sort(taxa_sums(Mar_leaf.fung_decon_pr_petr_S_nan.rar.class), TRUE)[1:10])
TOP10_P_Class_Names=data.frame("Phylum"=data.frame(tax_table(prune_taxa(P_TopCLASS_Field, Mar_leaf.fung_decon_pr_petr_S_nan.rar.class)))$Phylum,
                             "Class"=data.frame(tax_table(prune_taxa(P_TopCLASS_Field, Mar_leaf.fung_decon_pr_petr_S_nan.rar.class)))$Class,
                             "OTU"=P_TopCLASS_Field)



#Need to make an other species section

Mar_leaf.fung_decon_pr_petr_S_nan.rar.class_Names=merge(Mar_leaf.fung_decon_pr_petr_S_nan.rar.class_Names_raw,TOP10_P_Class_Names, by="OTU", all.x = T)

Mar_leaf.fung_decon_pr_petr_S_nan.rar.class_Names$Real_P_C=ifelse(is.na(Mar_leaf.fung_decon_pr_petr_S_nan.rar.class_Names$Class.y),
                                                                      "Other_spp",
                                                                      ifelse(Mar_leaf.fung_decon_pr_petr_S_nan.rar.class_Names$Class.x=="Unknown",
                                                                             paste(Mar_leaf.fung_decon_pr_petr_S_nan.rar.class_Names$Phylum.x,
                                                                                   "Other_spp",sep="."),                              
                                                                             as.character(Mar_leaf.fung_decon_pr_petr_S_nan.rar.class_Names$P_C)
                                                                      ))

#Transform the read counts to prop of total reads

Mar_leaf.fung_decon_pr_petr_S_nan.rar.class.prop=transform_sample_counts(Mar_leaf.fung_decon_pr_petr_S_nan.rar.class, function(x)x/sum(x))




Mar_leaf.fung_decon_pr_petr_S_nan.rar.class.prop_otu=as.data.frame(t(otu_table(Mar_leaf.fung_decon_pr_petr_S_nan.rar.class.prop)))
Mar_leaf.fung_decon_pr_petr_S_nan.rar.class.prop_otu2=merge(Mar_leaf.fung_decon_pr_petr_S_nan.rar.class_Names[,c("Real_P_C","OTU")],
                                                        Mar_leaf.fung_decon_pr_petr_S_nan.rar.class.prop_otu,by.y = "row.names",by.x = "OTU")
summary(Mar_leaf.fung_decon_pr_petr_S_nan.rar.class.prop_otu2)
Mar_leaf.fung_decon_pr_petr_S_nan.rar.class.prop_otu2$OTU=NULL
Mar_leaf.fung_decon_pr_petr_S_nan.rar.class.prop_otu_sum=Mar_leaf.fung_decon_pr_petr_S_nan.rar.class.prop_otu2%>%group_by(Real_P_C)%>%
  summarise_all(~sum(.))

Mar_leaf.fung_decon_pr_petr_S_nan.rar.class.prop_otu_sum_M=melt(Mar_leaf.fung_decon_pr_petr_S_nan.rar.class.prop_otu_sum,id="Real_P_C")

petri_positions=c("7/10/2018.No_rain","8/14/2018.Sterile_rain","8/14/2018.Live_Rain","7/21/2018.Rain","9/18/2018.Sterile_rain",
                  "9/18/2018.Live_Rain","8/21/2018.Rain")



Mar_leaf.fung_decon_pr_petr_S_nan.rar.class.prop_otu_sum_M2=Mar_leaf.fung_decon_pr_petr_S_nan.rar.class.prop_otu_sum_M%>%separate(variable,c("collectDate","rain_trt"),
                                                                                                                                      remove=F,sep = "[.]")

Mar_leaf.fung_decon_pr_petr_S_nan.rar.class.prop_otu_sum_M2$round=ifelse(Mar_leaf.fung_decon_pr_petr_S_nan.rar.class.prop_otu_sum_M2$collectDate=="7/10/2018","Seed",
                                                                           ifelse(Mar_leaf.fung_decon_pr_petr_S_nan.rar.class.prop_otu_sum_M2$collectDate=="8/14/2018"|
                                                                                    Mar_leaf.fung_decon_pr_petr_S_nan.rar.class.prop_otu_sum_M2$collectDate=="7/21/2018",
                                                                                  "Round1","Round2"))
(p_petri_color_class_F=ggplot(Mar_leaf.fung_decon_pr_petr_S_nan.rar.class.prop_otu_sum_M2,aes(x=factor(rain_trt,levels = c("No_rain","Sterile_rain","Live_Rain","Rain"),
                                                                                                   labels=c("Seeds","Sterile","Live","Rain")),
                                                                                          y=value,fill=Real_P_C))+
    geom_bar(aes( fill=Real_P_C), stat="identity", position="stack",color="black")+theme_cowplot()+
    theme(axis.text.y=element_text(size=26),axis.text.x=element_text(size=26),legend.position = "none",
          axis.title=element_text(size=32),panel.grid.major=element_blank(),legend.text = element_text(size=16), legend.title = element_text(size=20),
          strip.background = element_rect(fill = NA),strip.text = element_text(size=28))+xlab(NULL)+ylab("Proportion")+
    facet_grid(cols = vars(factor(round,levels = c("Seed","Round1","Round2"),labels=c("","Round1","Round2"))),space = "free",scales = "free")+
    guides(fill=guide_legend(title="Phyla"))+scale_fill_brewer(palette="Paired"))

#Stack_Class_rain_petri_nan
#1000x700


(p_petri_color_class_F_leg=ggplot(Mar_leaf.fung_decon_pr_petr_S_nan.rar.class.prop_otu_sum_M2,aes(x=factor(rain_trt,levels = c("No_rain","Sterile_rain","Live_Rain","Rain"),
                                                                                                       labels=c("Seeds","Sterile","Live","Rain")),
                                                                                              y=value,fill=Real_P_C))+
    geom_bar(aes( fill=Real_P_C), stat="identity", position="stack",color="black")+theme_cowplot()+
    theme(axis.text.y=element_text(size=26),axis.text.x=element_text(size=26),
          axis.title=element_text(size=32),panel.grid.major=element_blank(),legend.text = element_text(size=16), legend.title = element_text(size=20),
          strip.background = element_rect(fill = NA),strip.text = element_text(size=28))+xlab(NULL)+ylab("Proportion")+
    facet_grid(cols = vars(factor(round,levels = c("Seed","Round1","Round2"),labels=c("","Round1","Round2"))),space = "free",scales = "free")+
    guides(fill=guide_legend(title="Phyla"))+scale_fill_brewer(palette="Paired"))

#####Jaccard Petri Community Analyses######
Mar_leaf.fung_decon_pr_petr_S_nan.rar_J.ord=ordinate(Mar_leaf.fung_decon_pr_petr_S_nan.rar,method = "NMDS",distance =  "jaccard",binary = TRUE)
#Run 20 stress 0.2229492 
#*** Solution reached
#0.1753713      

#Creation of files for PERMANOVA analyses in Primer v6
Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_dis=phyloseq::distance(Mar_leaf.fung_decon_pr_petr_S_nan.rar,method =  "jaccard",binary = TRUE)
Mar_leaf.fung_decon_pr_petr_S_nan.rar_map=sample_data(Mar_leaf.fung_decon_pr_petr_S_nan.rar)


write.csv(as.matrix(Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_dis),here::here("R_file","Jaccard_Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_dis.csv"))

####Petri Jaccard Betadispersion#####

Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_betamod=betadisper(Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_dis, 
                                                         with(Mar_leaf.fung_decon_pr_petr_S_nan.rar_map,interaction(plant_type,collect_date,rain_trt)))

Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_betamod_dist=as.data.frame(Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_betamod$distances)
colnames(Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_betamod_dist)="betadisp"

Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_betamod_dist=merge(Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_betamod_dist,Mar_leaf.fung_decon_pr_petr_S_nan.rar_map,by="row.names")


#####Graphing the Betadispersion of Petri No nano####

brewer.pal(n = 12, name = "Paired")
head(Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_betamod_dist)


Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_betamod_dist$round=ifelse(Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_betamod_dist$rain_trt=="No_rain","Seed",ifelse(Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_betamod_dist$collect_date=="9/18/2018"|Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_betamod_dist$collect_date=="8/21/2018","Round2","Round1"))

Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_betamod_dist_sum=Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_betamod_dist%>%group_by(round,collect_date,rain_trt)%>%
  summarise_at(vars("betadisp"),list(~mean(.),~n(),se=~sd(.)/sqrt(n())))
unique(with(Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_betamod_dist,interaction(collect_date,rain_trt)))
date_trt=c("7/10/2018.No_rain","8/14/2018.Sterile_rain","8/14/2018.Live_Rain","7/21/2018.Rain",
           "9/18/2018.Sterile_rain","9/18/2018.Live_Rain","8/21/2018.Rain")



(petri_betadisp_nan_p=ggplot(Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_betamod_dist)+
    geom_point(aes(x=factor(rain_trt, levels = c("No_rain","Sterile_rain","Live_Rain","Rain"),
                            labels = c("Seed","Sterile","Live","Rain")),
                   y=betadisp,shape=rain_trt,fill=factor(rain_trt, levels = c("No_rain","Sterile_rain","Live_Rain","Rain")),
                   color=factor(rain_trt, levels = c("No_rain","Sterile_rain","Live_Rain","Rain"))),size=3,alpha=0.5)+
    geom_errorbar(data=Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_betamod_dist_sum,aes(x=factor(rain_trt, levels = c("No_rain","Sterile_rain","Live_Rain","Rain"),
                                                                                           labels = c("Seed","Sterile","Live","Rain")),
                                                                                  y=mean,ymax=mean+se,ymin=mean-se),
                  width=0.5)+facet_grid(~factor(round, levels= c("Seed","Round1","Round2"),
                                                labels = c("","Round1","Round2")),scale = "free_x",space = "free_x")+
    geom_point(data=Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_betamod_dist_sum,aes(x=factor(rain_trt, levels = c("No_rain","Sterile_rain","Live_Rain","Rain"),
                                                                                        labels = c("Seed","Sterile","Live","Rain")),y=mean,shape=rain_trt,
                                                                               fill=factor(rain_trt, levels = c("No_rain","Sterile_rain","Live_Rain","Rain")),
                                                                               color=factor(rain_trt, levels = c("No_rain","Sterile_rain","Live_Rain","Rain"))),
               size=10, stroke=1.5)+
    scale_shape_manual(values = c(24,13,25,21))+ylab("Beta-dispersion (Jaccard)")+
    scale_color_manual(values = c("#E31A1C","black", "black","black"),name=NULL)+
    scale_fill_manual(values = c("#E31A1C","#6A3D9A","#1F78B4","black"),name=NULL)+
    scale_x_discrete(name=NULL)+theme_cowplot()+
    theme(axis.title = element_text(size = 32),axis.text = element_text(size = 28),legend.title = element_blank(),legend.position = "none",
          strip.background = element_rect(fill = NA),strip.text = element_text(size = 32)))

#1200x800




Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_betamod_dist$trt_date=with(Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_betamod_dist,interaction(rain_trt,collect_date))



#Jaccarc bet disp

Petri_betadisp_J_mod_nan=lmer((betadisp)~trt_date+(1|gh_block),Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_betamod_dist)
plot(Petri_betadisp_J_mod_nan)
hist(resid(Petri_betadisp_J_mod_nan))
qqPlot(resid(Petri_betadisp_J_mod_nan))
shapiro.test(resid(Petri_betadisp_J_mod_nan))
#W = 0.97358, p-value = 0.5486

anova(Petri_betadisp_J_mod_nan)
#trt_date 0.092204 0.015367     6 14.208  4.3456 0.01076 *

#####POSTHOC TEST: Figure S6b####
emmeans(Petri_betadisp_J_mod_nan, pairwise~trt_date,adjust="fdr")





#####LFE Analyses the betadisp of Petri No nano####


#Petri 

Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_betamod_dist_p=subset(Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_betamod_dist,plant_type=="Petri")

#Let's make a round factor
unique(Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_betamod_dist_p$collect_date)
Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_betamod_dist_p$round=as.factor(ifelse(Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_betamod_dist_p$collect_date=="9/18/2018"|Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_betamod_dist_p$collect_date=="8/21/2018",
                                                                            "Round2","Round1"))
summary(Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_betamod_dist_p)
#Observed Richness

Petri_LFE_betadisp_J_mod=lmer((betadisp)~rain_trt*round+(1|gh_block),data=Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_betamod_dist_p)
plot(Petri_LFE_betadisp_J_mod)
hist(resid(Petri_LFE_betadisp_J_mod))
qqPlot(resid(Petri_LFE_betadisp_J_mod))
shapiro.test(resid(Petri_LFE_betadisp_J_mod))
#W = 0.96512, p-value = 0.5496

####Table S5 Beta-dispersion (Jaccard)####
anova(Petri_LFE_betadisp_J_mod)
#rain_trt       0.0160336 0.0160336     1 17.114  4.3017 0.05347 .
#round          0.0003384 0.0003384     1 19.711  0.0908 0.76632  
#rain_trt:round 0.0041289 0.0041289     1 17.114  1.1077 0.30721 


emmeans(Petri_LFE_betadisp_J_mod, pairwise~rain_trt|round)




#####PLOT: Figure S6####
#####Plot Petri combined betadisp graph####

(B_petri_betadisp_nan_p2=ggplot(Mar_leaf.fung_decon_pr_petr_S_nan.rar_betamod_dist)+
   geom_point(aes(x=factor(rain_trt, levels = c("No_rain","Sterile_rain","Live_Rain","Rain"),
                           labels = c("Seed","Sterile","Live","Rain")),
                  y=betadisp,shape=rain_trt,fill=factor(rain_trt, levels = c("No_rain","Sterile_rain","Live_Rain","Rain")),
                  color=factor(rain_trt, levels = c("No_rain","Sterile_rain","Live_Rain","Rain"))),size=3,alpha=0.5)+
   geom_errorbar(data=Mar_leaf.fung_decon_pr_petr_S_nan.rar_betamod_dist_sum,aes(x=factor(rain_trt, levels = c("No_rain","Sterile_rain","Live_Rain","Rain"),
                                                                                          labels = c("Seed","Sterile","Live","Rain")),
                                                                                 y=mean,ymax=mean+se,ymin=mean-se),
                 width=0.5)+facet_grid(~factor(round, levels= c("Seed","Round1","Round2"),
                                               labels = c("","Round1","Round2")),scale = "free_x",space = "free_x")+
   geom_point(data=Mar_leaf.fung_decon_pr_petr_S_nan.rar_betamod_dist_sum,aes(x=factor(rain_trt, levels = c("No_rain","Sterile_rain","Live_Rain","Rain"),
                                                                                       labels = c("Seed","Sterile","Live","Rain")),y=mean,shape=rain_trt,
                                                                              fill=factor(rain_trt, levels = c("No_rain","Sterile_rain","Live_Rain","Rain")),
                                                                              color=factor(rain_trt, levels = c("No_rain","Sterile_rain","Live_Rain","Rain"))),
              size=10, stroke=1.5)+
   scale_shape_manual(values = c(24,13,25,21))+ylab("Beta-dispersion\n(Bray-Curtis)")+
   scale_color_manual(values = c("#E31A1C","black", "black","black"),name=NULL)+
   scale_fill_manual(values = c("#E31A1C","#6A3D9A","#1F78B4","black"),name=NULL)+
   scale_x_discrete(name=NULL)+theme_cowplot()+
   theme(axis.title = element_text(size = 32),axis.text = element_text(size = 28),legend.title = element_blank(),legend.position = "none",
         strip.background = element_rect(fill = NA),strip.text = element_text(size = 32)))

(J_petri_betadisp_nan_p2=ggplot(Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_betamod_dist)+
    geom_point(aes(x=factor(rain_trt, levels = c("No_rain","Sterile_rain","Live_Rain","Rain"),
                            labels = c("Seed","Sterile","Live","Rain")),
                   y=betadisp,shape=rain_trt,fill=factor(rain_trt, levels = c("No_rain","Sterile_rain","Live_Rain","Rain")),
                   color=factor(rain_trt, levels = c("No_rain","Sterile_rain","Live_Rain","Rain"))),size=3,alpha=0.5)+
    geom_errorbar(data=Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_betamod_dist_sum,aes(x=factor(rain_trt, levels = c("No_rain","Sterile_rain","Live_Rain","Rain"),
                                                                                             labels = c("Seed","Sterile","Live","Rain")),
                                                                                    y=mean,ymax=mean+se,ymin=mean-se),
                  width=0.5)+facet_grid(~factor(round, levels= c("Seed","Round1","Round2"),
                                                labels = c("","Round1","Round2")),scale = "free_x",space = "free_x")+
    geom_point(data=Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_betamod_dist_sum,aes(x=factor(rain_trt, levels = c("No_rain","Sterile_rain","Live_Rain","Rain"),
                                                                                          labels = c("Seed","Sterile","Live","Rain")),y=mean,shape=rain_trt,
                                                                                 fill=factor(rain_trt, levels = c("No_rain","Sterile_rain","Live_Rain","Rain")),
                                                                                 color=factor(rain_trt, levels = c("No_rain","Sterile_rain","Live_Rain","Rain"))),
               size=10, stroke=1.5)+
    scale_shape_manual(values = c(24,13,25,21))+ylab("Beta-dispersion\n(Jaccard)")+
    scale_color_manual(values = c("#E31A1C","black", "black","black"),name=NULL)+
    scale_fill_manual(values = c("#E31A1C","#6A3D9A","#1F78B4","black"),name=NULL)+
    scale_x_discrete(name=NULL)+theme_cowplot()+
    theme(axis.title = element_text(size = 32),axis.text = element_text(size = 28),legend.title = element_blank(),legend.position = "none",
          strip.background = element_rect(fill = NA),strip.text = element_text(size = 32)))

plot_grid(B_petri_betadisp_nan_p2,J_petri_betadisp_nan_p2,nrow = 2,align = "v")

#1000*1000
#Petri_betadsp_Bray_Jacc_raw


#####Graphing Jaccard Petri NMDS######
#Hull for the Jaccard NMDS
Mar_leaf.fung_decon_pr_petr_S_nan.rar_map=sample_data(Mar_leaf.fung_decon_pr_petr_S_nan.rar)
summary(Mar_leaf.fung_decon_pr_petr_S_nan.rar_map)
Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_coordinates=merge(Mar_leaf.fung_decon_pr_petr_S_nan.rar_J.ord$points,Mar_leaf.fung_decon_pr_petr_S_nan.rar_map, by="row.names")

Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_coordinates$plant_data_grp=ifelse(Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_coordinates$plant_type=="Rain"&
                                                                      Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_coordinates$collect_date=="7/21/2018",
                                                                    "Round1_Rain",ifelse(Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_coordinates$plant_type=="Rain"&
                                                                                           Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_coordinates$collect_date=="8/21/2018",
                                                                                         "Round2_Rain", 
                                                                                                               ifelse(Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_coordinates$rain_trt=="Live_Rain",
                                                                                                                      ifelse(Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_coordinates$collect_date=="9/18/2018","Round2_Live","Round1_Live"),
                                                                                                                      ifelse(Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_coordinates$rain_trt=="Sterile_rain",
                                                                                                                             ifelse(Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_coordinates$collect_date=="9/18/2018","Round2_Sterile","Round1_Sterile"),"Seed"))))


unique(Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_coordinates$plant_data_grp)

J_grp1000_Round1_Live <- Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_coordinates[Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_coordinates$plant_data_grp == "Round1_Live", 
                                                                     ][chull(Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_coordinates[
                                                                       Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_coordinates$
                                                                         plant_data_grp =="Round1_Live", c("MDS1", "MDS2")]), ] 
nrow(J_grp1000_Round1_Live)
#5
J_grp1000_Round2_Live <- Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_coordinates[Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_coordinates$plant_data_grp == "Round2_Live", 
                                                                     ][chull(Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_coordinates[
                                                                       Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_coordinates$
                                                                         plant_data_grp =="Round2_Live", c("MDS1", "MDS2")]), ]  

nrow(J_grp1000_Round2_Live)
#4


J_grp1000_Round1_Sterile <- Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_coordinates[Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_coordinates$plant_data_grp == "Round1_Sterile", 
                                                                        ][chull(Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_coordinates[
                                                                          Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_coordinates$
                                                                            plant_data_grp =="Round1_Sterile", c("MDS1", "MDS2")]), ]  

nrow(J_grp1000_Round1_Sterile)
#5

J_grp1000_Round2_Sterile <- Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_coordinates[Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_coordinates$plant_data_grp == "Round2_Sterile", 
                                                                        ][chull(Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_coordinates[
                                                                          Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_coordinates$
                                                                            plant_data_grp =="Round2_Sterile", c("MDS1", "MDS2")]), ]  

nrow(J_grp1000_Round2_Sterile)
#4


J_grp1000_P_Seed <- Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_coordinates[Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_coordinates$plant_data_grp == "Seed", 
                                                                ][chull(Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_coordinates[
                                                                  Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_coordinates$
                                                                    plant_data_grp =="Seed", c("MDS1", "MDS2")]), ]  

nrow(J_grp1000_P_Seed)
#3


J_grp1000_P_Round1_Rain <- Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_coordinates[Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_coordinates$plant_data_grp == "Round1_Rain", 
                                                                       ][chull(Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_coordinates[
                                                                         Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_coordinates$
                                                                           plant_data_grp =="Round1_Rain", c("MDS1", "MDS2")]), ]  

nrow(J_grp1000_P_Round1_Rain)
#3

J_grp1000_P_Round2_Rain <- Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_coordinates[Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_coordinates$plant_data_grp == "Round2_Rain", 
                                                                       ][chull(Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_coordinates[
                                                                         Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_coordinates$
                                                                           plant_data_grp =="Round2_Rain", c("MDS1", "MDS2")]), ]  

nrow(J_grp1000_P_Round2_Rain)
#4

J_hull1000_Petri=rbind(J_grp1000_Round1_Live,J_grp1000_Round2_Live,J_grp1000_Round1_Sterile,J_grp1000_Round2_Sterile,
                     J_grp1000_P_Round1_Rain, J_grp1000_P_Round2_Rain,J_grp1000_P_Seed)
nrow(J_hull1000_Petri)
#28


#In order to use facet I need to duplicate seed
J_hull1000_Petri_seed=subset(J_hull1000_Petri,plant_type=="Seed")
Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_coordinates_seed=subset(Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_coordinates,plant_type=="Seed")

Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_coordinates$round=ifelse(Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_coordinates$collect_date=="9/18/2018"|
                                                                   Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_coordinates$collect_date=="8/21/2018",
                                                               "Round2","Round1")
J_hull1000_Petri$round=ifelse(J_hull1000_Petri$collect_date=="9/18/2018"|J_hull1000_Petri$collect_date=="8/21/2018",
                            "Round2","Round1")

Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_coordinates_seed$round=rep("Round2")

J_hull1000_Petri_seed$round=rep("Round2")


Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_coordinates2=rbind(Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_coordinates_seed,
                                                         Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_coordinates)

J_hull1000_Petri2=rbind(J_hull1000_Petri,J_hull1000_Petri_seed)






unique(Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_coordinates$plant_data_grp)

P_plant_data_grp_order=c("Seed","Round1_Sterile","Round2_Sterile","Round1_Live","Round2_Live","Round1_Rain","Round2_Rain")
length(P_plant_data_grp_order)
unique(Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_coordinates$rain_trt)

(petri_nmds_jacc=ggplot(Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_coordinates2,aes(x= MDS1, y=MDS2))+
  geom_polygon(data=J_hull1000_Petri2,aes(x=MDS1,y=MDS2,fill=factor(plant_data_grp,levels = P_plant_data_grp_order),
                                       group=factor(plant_data_grp,levels = P_plant_data_grp_order)), alpha=0.5)+
  geom_point(size=4,aes(fill=factor(plant_data_grp,levels = P_plant_data_grp_order),
                        shape=rain_trt, color=factor(plant_data_grp,levels = P_plant_data_grp_order)),stroke=2)+
  scale_shape_manual(values = c(24,13,21,25),name=NULL)+
  scale_color_manual(values = c("#E31A1C","black", "black","black","black","black","black","black","black"),name=NULL)+
  scale_fill_manual(values = c("#E31A1C","#6A3D9A","#6A3D9A","#1F78B4","#1F78B4","black","black"),name=NULL)+
  facet_nested(.~round)+theme_classic()+theme(axis.title.x = element_blank(),axis.title.y = element_text(size = 28),
                                              axis.text = element_text(size = 30),panel.border =  element_rect(fill = NA,linetype="solid"),
                                              legend.position = "none",strip.text = element_text(size = 36),
                                              strip.placement = "outside",strip.background = element_rect(linetype="blank"),
                                              panel.spacing = unit(0, "lines")))

##900x800



#####PLOT: Figure S1####
#####Plot Petri combined NMDS####

(petri_nmds_bray=ggplot(Mar_leaf.fung_decon_pr_petr_S_nan.rar_coordinates2,aes(x= MDS1, y=MDS2))+
    geom_polygon(data=hull1000_Petri2,aes(x=MDS1,y=MDS2,fill=factor(plant_data_grp,levels = P_plant_data_grp_order),
                                          group=factor(plant_data_grp,levels = P_plant_data_grp_order)), alpha=0.5)+
    geom_point(size=4,aes(fill=factor(plant_data_grp,levels = P_plant_data_grp_order),
                          shape=rain_trt, color=factor(plant_data_grp,levels = P_plant_data_grp_order)),stroke=2)+
    scale_shape_manual(values = c(24,13,21,25),name=NULL)+
    scale_x_continuous(limits = c(-1,1.1),breaks = c(-0.6,0,0.6))+
    scale_y_continuous(limits = c(-1.1,1.1),name = "NMDS2",breaks = c(-0.6,0,0.6))+
    scale_color_manual(values = c("#E31A1C","black", "black","black","black","black","black","black","black"),name=NULL)+
    scale_fill_manual(values = c("#E31A1C","#6A3D9A","#6A3D9A","#1F78B4","#1F78B4","black","black"),name=NULL)+
    facet_nested(.~round)+theme_classic()+theme(axis.title.x = element_blank(),axis.title.y = element_text(size = 28),
                                                axis.text.x = element_blank(),axis.text.y = element_text(size = 30),
                                                panel.border =  element_rect(fill = NA,linetype="solid"),
                                                legend.position = "none",strip.text = element_text(size = 36),
                                                strip.placement = "outside",strip.background = element_rect(linetype="blank"),
                                                panel.spacing = unit(0, "lines")))
(petri_nmds_jacc=ggplot(Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_coordinates2,aes(x= MDS1, y=MDS2))+
    geom_polygon(data=J_hull1000_Petri2,aes(x=MDS1,y=MDS2,fill=factor(plant_data_grp,levels = P_plant_data_grp_order),
                                            group=factor(plant_data_grp,levels = P_plant_data_grp_order)), alpha=0.5)+
    geom_point(size=4,aes(fill=factor(plant_data_grp,levels = P_plant_data_grp_order),
                          shape=rain_trt, color=factor(plant_data_grp,levels = P_plant_data_grp_order)),stroke=2)+
    scale_shape_manual(values = c(24,13,21,25),name=NULL)+
    scale_color_manual(values = c("#E31A1C","black", "black","black","black","black","black","black","black"),name=NULL)+
    scale_fill_manual(values = c("#E31A1C","#6A3D9A","#6A3D9A","#1F78B4","#1F78B4","black","black"),name=NULL)+
    scale_x_continuous(limits = c(-1,1.1),name = "NMDS1", breaks = c(-0.6,0,0.6))+
    scale_y_continuous(limits = c(-1.1,1.1),name = "NMDS2",breaks = c(-0.6,0,0.6))+
    facet_nested(.~round)+theme_classic()+theme(axis.title = element_text(size = 28),
                                                axis.text = element_text(size = 30),panel.border =  element_rect(fill = NA,linetype="solid"),
                                                legend.position = "none",strip.text = element_blank(),
                                                strip.placement = "none",strip.background = element_rect(linetype="blank"),
                                                panel.spacing = unit(0, "lines")))


plot_grid(petri_nmds_bray,petri_nmds_jacc, nrow = 2, rel_heights = c(.9,1))
#Petri_NMDS_comb_raw

#####Petri Diversity####



ntaxa(Mar_leaf.fung_decon_pr_petr_S_nan.rar)
#944
sum(taxa_sums(Mar_leaf.fung_decon_pr_petr_S_nan.rar))
#35000
min(sample_sums(Mar_leaf.fung_decon_pr_petr_S_nan.rar))
#1000
max(sample_sums(Mar_leaf.fung_decon_pr_petr_S_nan.rar))
#1000

Mar_leaf.fung_decon_pr_petr_S_nan.rar_div=estimate_richness(Mar_leaf.fung_decon_pr_petr_S_nan.rar)
Mar_leaf_rar_petri_div_map_nan=merge(Mar_leaf.fung_decon_pr_petr_S_nan.rar_div,sample_data(Mar_leaf.fung_decon_pr_petr_S_nan.rar),by="row.names")


#####PLOT: Figure S4####
#####Graphing the diversity of Petri No nano####


Mar_leaf_rar_petri_div_map_nan$round=ifelse(Mar_leaf_rar_petri_div_map_nan$rain_trt=="No_rain","Seed",ifelse(Mar_leaf_rar_petri_div_map_nan$collect_date=="9/18/2018"|Mar_leaf_rar_petri_div_map_nan$collect_date=="8/21/2018","Round2","Round1"))

Mar_leaf_rar_petri_div_map_nan_sum=Mar_leaf_rar_petri_div_map_nan%>%group_by(round,collect_date,rain_trt)%>%
  summarise_at(vars("Observed","InvSimpson"),list(~mean(.),~n(),se=~sd(.)/sqrt(n())))
unique(with(Mar_leaf_rar_petri_div_map_nan,interaction(collect_date,rain_trt)))
date_trt=c("7/10/2018.No_rain","8/14/2018.Sterile_rain","8/14/2018.Live_Rain","7/21/2018.Rain",
           "9/18/2018.Sterile_rain","9/18/2018.Live_Rain","8/21/2018.Rain")


#Richness
(petri_rich_nan=ggplot(Mar_leaf_rar_petri_div_map_nan)+geom_point(aes(x=factor(rain_trt, levels = c("No_rain","Sterile_rain","Live_Rain","Rain"),
                                                               labels = c("Seed","Sterile","Live","Rain")),
                                                  y=Observed,shape=rain_trt,
                                                  fill=factor(rain_trt, levels = c("No_rain","Sterile_rain","Live_Rain","Rain")),
                                                  color=factor(rain_trt, levels = c("No_rain","Sterile_rain","Live_Rain","Rain"))),size=3,alpha=0.5)+
  geom_errorbar(data=Mar_leaf_rar_petri_div_map_nan_sum,aes(x=factor(rain_trt, levels = c("No_rain","Sterile_rain","Live_Rain","Rain"),
                                                                     labels = c("Seed","Sterile","Live","Rain")),
                                                        y=Observed_mean,ymax=Observed_mean+Observed_se,ymin=Observed_mean-Observed_se),
                width=0.5)+facet_grid(~factor(round, levels= c("Seed","Round1","Round2"),
                                              labels = c("","Round1","Round2")),scale = "free_x",space = "free_x")+
  geom_point(data=Mar_leaf_rar_petri_div_map_nan_sum,aes(x=factor(rain_trt, levels = c("No_rain","Sterile_rain","Live_Rain","Rain"),
                                                                  labels = c("Seed","Sterile","Live","Rain")),y=Observed_mean,shape=rain_trt,
                                                                              fill=factor(rain_trt, levels = c("No_rain","Sterile_rain","Live_Rain","Rain")),
                                                                              color=factor(rain_trt, levels = c("No_rain","Sterile_rain","Live_Rain","Rain"))),
             size=10, stroke=1.5)+
  scale_shape_manual(values = c(24,13,25,21))+ylab("Richness (ZOTUs)")+
  scale_color_manual(values = c("#E31A1C","black", "black","black"),name=NULL)+
  scale_fill_manual(values = c("#E31A1C","#6A3D9A","#1F78B4","black"),name=NULL)+
  scale_x_discrete(name=NULL)+theme_cowplot()+
  theme(axis.title = element_text(size = 32),axis.text = element_text(size = 28),legend.title = element_blank(),legend.position = "none",
        strip.background = element_rect(fill = NA),strip.text = element_text(size = 32)))

#1200x800
#geom_line(data=Mar_leaf_rar_petri_div_map_nan_sum,aes(x=factor(interaction(collect_date,rain_trt),levels = date_trt),group=collect_date,
#y=Observed_mean),size=2)+


#Inverse Simpson



(petri_invSimp_nan=ggplot(Mar_leaf_rar_petri_div_map_nan)+geom_point(aes(x=factor(rain_trt, levels = c("No_rain","Sterile_rain","Live_Rain","Rain"),
                                                                               labels = c("Seed","Sterile","Live","Rain")),
                                                                      y=InvSimpson,shape=rain_trt,
                                                                      fill=factor(rain_trt, levels = c("No_rain","Sterile_rain","Live_Rain","Rain")),
                                                                      color=factor(rain_trt, levels = c("No_rain","Sterile_rain","Live_Rain","Rain"))),size=3,alpha=0.5)+
    geom_errorbar(data=Mar_leaf_rar_petri_div_map_nan_sum,aes(x=factor(rain_trt, levels = c("No_rain","Sterile_rain","Live_Rain","Rain"),
                                                                       labels = c("Seed","Sterile","Live","Rain")),
                                                              y=InvSimpson_mean,ymax=InvSimpson_mean+InvSimpson_se,ymin=InvSimpson_mean-InvSimpson_se),
                  width=0.5)+facet_grid(~factor(round, levels= c("Seed","Round1","Round2"),
                                                labels = c("","Round1","Round2")),scale = "free_x",space = "free_x")+
    geom_point(data=Mar_leaf_rar_petri_div_map_nan_sum,aes(x=factor(rain_trt, levels = c("No_rain","Sterile_rain","Live_Rain","Rain"),
                                                                    labels = c("Seed","Sterile","Live","Rain")),y=InvSimpson_mean,shape=rain_trt,
                                                           fill=factor(rain_trt, levels = c("No_rain","Sterile_rain","Live_Rain","Rain")),
                                                           color=factor(rain_trt, levels = c("No_rain","Sterile_rain","Live_Rain","Rain"))),
               size=10, stroke=1.5)+
    scale_shape_manual(values = c(24,13,25,21))+ylab("inverse Simpson")+
    scale_color_manual(values = c("#E31A1C","black", "black","black"),name=NULL)+
    scale_fill_manual(values = c("#E31A1C","#6A3D9A","#1F78B4","black"),name=NULL)+
    scale_x_discrete(name=NULL)+theme_cowplot()+
    theme(axis.title = element_text(size = 32),axis.text = element_text(size = 28),legend.title = element_blank(),legend.position = "none",
          strip.background = element_rect(fill = NA),strip.text = element_text(size = 32)))


#1200x800

plot_grid(petri_rich_nan,petri_invSimp_nan,nrow = 2, align="v")
#1000*1000

Mar_leaf_rar_petri_div_map_nan$trt_date=with(Mar_leaf_rar_petri_div_map_nan,interaction(collect_date,rain_trt))
#Observed Richness

Petri_seed_rich_mod_nan=lmer(sqrt(Observed)~trt_date+(1|gh_block),Mar_leaf_rar_petri_div_map_nan)
plot(Petri_seed_rich_mod_nan)
hist(resid(Petri_seed_rich_mod_nan))
qqPlot(resid(Petri_seed_rich_mod_nan))
shapiro.test(resid(Petri_seed_rich_mod_nan))
#W = 0.95952, p-value = 0.2213

anova(Petri_seed_rich_mod_nan)
#trt_date 340.22  56.703     6    28  27.491 1.648e-10 ***
#no change
#####POSTHOC TEST: Figure S4a####
emmeans(Petri_seed_rich_mod_nan, pairwise~trt_date,adjust="fdr")





#invSimpson

Petri_seed_invS_mod_nan=lmer((InvSimpson)~trt_date+(1|gh_block),Mar_leaf_rar_petri_div_map_nan)
plot(Petri_seed_invS_mod_nan)
hist(resid(Petri_seed_invS_mod_nan))
qqPlot(resid(Petri_seed_invS_mod_nan))
shapiro.test(resid(Petri_seed_invS_mod_nan))
#W = 0.9816, p-value = 0.8101

anova(Petri_seed_invS_mod_nan)
#trt_date 3155.9  525.98     6 18.209  6.5779 0.0007946 ***

#####POSTHOC TEST: Figure S4b####
emmeans(Petri_seed_invS_mod_nan, pairwise~trt_date)


#####LFE Analyses the diversity of Petri No nano####


#Petri 

Mar_leaf_rar_petri_div_map_P_nan=subset(Mar_leaf_rar_petri_div_map_nan,plant_type=="Petri"&rain_trt!="Nano")

#Let's make a round factor
unique(Mar_leaf_rar_petri_div_map_P_nan$collect_date)
Mar_leaf_rar_petri_div_map_P_nan$round=as.factor(ifelse(Mar_leaf_rar_petri_div_map_P_nan$collect_date=="9/18/2018"|Mar_leaf_rar_petri_div_map_P_nan$collect_date=="8/21/2018",
                                          "Round2","Round1"))
summary(Mar_leaf_rar_petri_div_map_P_nan)
#Observed Richness

Petri_rich_mod_nan=lmer((Observed)~rain_trt*round+(1|gh_block),data=Mar_leaf_rar_petri_div_map_P_nan)
#boundary (singular) fit: see ?isSingular
hist(resid(Petri_rich_mod_nan))
qqPlot(resid(Petri_rich_mod_nan))
shapiro.test(resid(Petri_rich_mod_nan))
#W = 0.9362, p-value = 0.1342

####Table S5 Richness####
anova(Petri_rich_mod_nan)
#rain_trt       925.04  925.04     1    20  2.1376 0.1593
#round          145.04  145.04     1    20  0.3352 0.5691
#rain_trt:round  26.04   26.04     1    20  0.0602 0.8087
#no change

emmeans(Petri_rich_mod_nan, pairwise~rain_trt|round)


ggplot(Mar_leaf_rar_petri_div_map_P_nan,aes(x=interaction(rain_trt,start_date),y=Observed))+geom_point()+geom_boxplot()


#invSimp

Petri_invSimp_mod_nan=lmer((InvSimpson)~rain_trt*round+(1|gh_block),data=Mar_leaf_rar_petri_div_map_P_nan)
##boundary (singular) fit: see ?isSingular
hist(resid(Petri_invSimp_mod_nan))
qqPlot(resid(Petri_invSimp_mod_nan))
shapiro.test(resid(Petri_invSimp_mod_nan))
#W = 0.9821, p-value = 0.9312

####Table S5 Inverse Simpson####
anova(Petri_invSimp_mod_nan)
#rain_trt       91.446  91.446     1    20  2.0565 0.1670
#round          20.833  20.833     1    20  0.4685 0.5015
#rain_trt:round 93.170  93.170     1    20  2.0953 0.1632

#no change
emmeans(Petri_invSimp_mod_nan, pairwise~rain_trt|round)



#Shannon

Petri_Shan_mod_nan=lmer((Shannon)~rain_trt*round+(1|gh_block),data=Mar_leaf_rar_petri_div_map_P_nan)
#boundary (singular) fit: see ?isSingular
hist(resid(Petri_Shan_mod_nan))
qqPlot(resid(Petri_Shan_mod_nan))
shapiro.test(resid(Petri_Shan_mod_nan))
#W = 0.93728, p-value = 0.1419

anova(Petri_Shan_mod_nan)
#rain_trt       0.39536 0.39536     1    20  0.8129 0.3780
#round          0.03611 0.03611     1    20  0.0742 0.7880
#rain_trt:round 0.93046 0.93046     1    20  1.9130 0.1819

#no change

emmeans(Petri_Shan_mod_nan, pairwise~rain_trt|round)


ggplot(Mar_leaf_rar_petri_div_map_P_nan,aes(x=interaction(rain_trt,start_date),y=Shannon))+geom_point()+geom_boxplot()



#####Petri Plant Analyses####
#Germination and Leaf fungal colonization

Petri_plant=read.csv(here::here("USEARCHv11","Marshall_petri_exp_plant_data.csv"),header = T)
summary(Petri_plant)


#final germination

Petri_plant_end=subset(Petri_plant,measure_date=="9/18/2018"|measure_date=="8/14/2018")
nrow(Petri_plant_end)
#32


#no nano

Petri_plant_end_nan=subset(Petri_plant_end,rain_trt!="Nano")

Petri_plant_end_sum_nan=subset(Petri_plant_end_sum,rain_trt!="Nano")



#Germination
hist(Petri_plant_end_nan$germ_seed)
qqPlot(Petri_plant_end_nan$germ_seed)

Petri_germ_end_mod_nan=lmer((germ_seed)~rain_trt*collect_date+(1|field_plot),data = Petri_plant_end_nan)
qqPlot(resid(Petri_germ_end_mod_nan))
hist(resid(Petri_germ_end_mod_nan))

shapiro.test(resid(Petri_germ_end_mod_nan))
#0.6724

####Table S6 Germination rate####
anova(Petri_germ_end_mod_nan, type = 3)
#Nada sig



(petri_germ=ggplot(Petri_plant_end_nan)+geom_point(aes(x=factor(rain_trt,levels = c("Sterile_rain","Live_Rain")),shape=rain_trt,
                                                       y=germ_seed, fill=factor(rain_trt,levels = c("Sterile_rain","Live_Rain"))),size=3, alpha=0.1)+
    geom_errorbar(data=Petri_plant_end_sum_nan,aes(x=factor(rain_trt,levels = c("Sterile_rain","Live_Rain")),y=germ_seed_mean,
                                                   ymax=germ_seed_mean+germ_seed_se,ymin=germ_seed_mean-germ_seed_se),width=0.5)+
    geom_point(data=Petri_plant_end_sum_nan,aes(x=factor(rain_trt,levels = c("Sterile_rain","Live_Rain")),shape=rain_trt,
                                                y=germ_seed_mean, fill=factor(rain_trt,levels = c("Sterile_rain","Live_Rain"))),size=10)+
    scale_shape_manual(values = c(24,25))+facet_grid(cols = vars(factor(collect_date,labels = c("Round1","Round2"))))+
    scale_fill_manual(values = c("#6A3D9A","#1F78B4","#6A3D9A","#1F78B4"),name=NULL)+
    scale_x_discrete(labels=c("Sterile","Live"),name=NULL)+scale_y_continuous(name="Number of\ngerminants")+
    theme_cowplot()+theme(axis.text = element_text(size = 28),axis.title = element_text(size = 32),legend.position = "none",
                          strip.background = element_rect(fill = NA),strip.text = element_text(size = 32)))

#petri_germination_dot_error_nan
#800x600



#Fungal colonization
hist(Petri_plant_end_nan$fung_seed)
qqPlot(Petri_plant_end_nan$fung_seed)
hist((Petri_plant_end_nan$fung_seed+1)^(1/3))
qqPlot((Petri_plant_end_nan$fung_seed+1)^(1/3))

Petri_fung_col_end_mod_nan=lmer(sqrt(fung_seed)~rain_trt*collect_date+(1|field_plot),data = Petri_plant_end_nan)
#boundary (singular) fit: see ?isSingular
plot(Petri_fung_col_end_mod_nan)
qqPlot(resid(Petri_fung_col_end_mod_nan))
hist(resid(Petri_fung_col_end_mod_nan))
shapiro.test(resid(Petri_fung_col_end_mod_nan))
#0.0591

####Table S6 Fungal colonization####
anova(Petri_fung_col_end_mod_nan, type = 3)
#collect_date          2.29044 2.29044     1    20  4.3761 0.04941 *

emmeans(Petri_fung_col_end_mod_nan,pairwise~collect_date|rain_trt)

emmeans(Petri_fung_col_end_mod_nan,pairwise~rain_trt|collect_date)

emmeans(Petri_fung_col_end_mod_nan,pairwise~rain_trt*collect_date, adjust="fdr")

#no nano

Petri_plant_end_nan=subset(Petri_plant_end,rain_trt!="Nano")

Petri_plant_end_sum_nan=subset(Petri_plant_end_sum,rain_trt!="Nano")



(petri_fung_col=ggplot(Petri_plant_end_nan)+geom_point(aes(x=factor(rain_trt,levels = c("Sterile_rain","Live_Rain")),shape=rain_trt,
                                                       y=fung_seed, fill=factor(rain_trt,levels = c("Sterile_rain","Live_Rain"))),size=3,alpha=0.1)+
    geom_errorbar(data=Petri_plant_end_sum_nan,aes(x=factor(rain_trt,levels = c("Sterile_rain","Live_Rain")),y=fung_seed_mean,
                                                   ymax=fung_seed_mean+fung_seed_se,ymin=fung_seed_mean-fung_seed_se),width=0.5)+
    geom_point(data=Petri_plant_end_sum_nan,aes(x=factor(rain_trt,levels = c("Sterile_rain","Live_Rain")),shape=rain_trt,
                                                y=fung_seed_mean, fill=factor(rain_trt,levels = c("Sterile_rain","Live_Rain"))),size=10)+
    scale_shape_manual(values = c(24,25))+facet_grid(cols = vars(factor(collect_date,labels = c("Round1","Round2"))))+
    scale_fill_manual(values = c("#6A3D9A","#1F78B4","#6A3D9A","#1F78B4"),name=NULL)+
    scale_x_discrete(labels=c("Sterile","Live"),name=NULL)+scale_y_continuous(name="Number of\nseeds colonized")+
    theme_cowplot()+theme(axis.text = element_text(size = 28),axis.title = element_text(size = 32),legend.position = "none",
                          strip.background = element_rect(fill = NA),strip.text = element_text(size = 32)))

#petri_fungi_colon_dot_error_nan
#800x600
#####PLOT: Figure S8####
#####Seed Germ and Colon Combined Graphs#####

plot_grid(petri_germ,petri_fung_col,nrow = 2, align="v")
#900*900




#####Petri Published leaf endophyte sequence analyses####



#Run using Usearch 

#cd HardDrive/FungiRainLeaf2019/R_file/
#The program needs for the codons to be uppercase
#awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' rep_set.FunRainLeaf2019_v4.2.2020_ZOTU_rar1000.fna.fung_decon_rar_pruned_phyloseq_obj.RData > CAP_rep_set.FunRainLeaf2019_v4.2.2020_ZOTU_rar1000.fna.fung_decon_rar_pruned_phyloseq_obj.RData
#~/HardDrive/Sciencey_Program/Old_usearch_v/usearch10.0.240_win32.exe -usearch_global CAP_rep_set.FunRainLeaf2019_v4.2.2020_ZOTU_rar1000.fna.fung_decon_rar_pruned_phyloseq_obj.RData -db switchgrass_compiled_ITS_database_v2.fasta -id 0.97 -strand both -maxaccepts 0 -maxhits 10 -matched MATCHED_pub_seq_rep_set.FunRainLeaf2019_ZOTU_rar1000.fa -notmatched NOT_matched_pub_seq_rep_set.FunRainLeaf2019_ZOTU_rar1000.fa -userout TBL_pub_seq_rep_set.FunRainLeaf2019_ZOTU_rar1000.txt -userfields query+target+id+mid+bits+evalue+ql+ts+qlor+qhir+tlor+thir
#00:02 14Mb    100.0% Searching, 5.1% matched

publ_matched_rar_OTUs=read.delim(here::here("R_file","TBL_pub_seq_rep_set.FunRainLeaf2019_ZOTU_v4.2.2020_rar1000.txt"), header = F)
head(publ_matched_rar_OTUs)

colnames(publ_matched_rar_OTUs)=c("query","target","id","mid","bits","evalue","ql","ts","qlor","qhir","tlor","thir")

nrow(publ_matched_rar_OTUs)
#468
#Subset to have only the top hiting taxa 
publ_matched_rar_OTUs_top_hits=publ_matched_rar_OTUs %>% group_by(query) %>% slice(which(id == max(id)))

nrow(publ_matched_rar_OTUs_top_hits)
#258



Mar_leaf_rain.fung_decon_rar_Mar=subset_samples(Mar_leaf_rain.fung_decon_rar, sub_proj=="Marshall")
Mar_leaf_rain.fung_decon_rar_Mar_petri=subset_samples(Mar_leaf_rain.fung_decon_rar_Mar, plant_type!="Adult"&plant_type!="Seedling")

Mar_leaf.fung_decon_pr_petr_S.rar=subset_samples(Mar_leaf_rain.fung_decon_rar_Mar_petri,plant_type!="Rain"|collect_date=="8/21/2018"|
                                                   collect_date=="7/21/2018")

nsamples(Mar_leaf.fung_decon_pr_petr_S.rar)
#43


Mar_leaf.fung_decon_pr_petr_S.rar=prune_taxa(taxa_sums(Mar_leaf.fung_decon_pr_petr_S.rar) > 0, Mar_leaf.fung_decon_pr_petr_S.rar)
ntaxa(Mar_leaf.fung_decon_pr_petr_S.rar)
#1027
sum(taxa_sums(Mar_leaf.fung_decon_pr_petr_S.rar))
#43000
min(sample_sums(Mar_leaf.fung_decon_pr_petr_S.rar))
#1000
max(sample_sums(Mar_leaf.fung_decon_pr_petr_S.rar))
#1000


Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo=merge(otu_table(Mar_leaf.fung_decon_pr_petr_S.rar),
                                                      publ_matched_rar_OTUs_top_hits,
                                                      by.x="row.names", by.y="query",all.x=T)
head(Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo)
colnames(Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo)
Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo$OTU_sum=rowSums(Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo[,2:44])
Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_sort=Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo[order(-Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo$OTU_sum),]

unique(Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo$target)
#30

#Rplace NAs with not matched 

Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo[is.na(Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo)] <- "Not_matched"

#Let's look at main taxa

colnames(Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo)
summary(Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo)
sapply(Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo, class)

#Group by interaction and indicator type

#Load in the names table
publ_matchednames=read.csv(here::here("R_file","Matched_pub_seq_names_TBL.csv"), header = T)
colnames(publ_matchednames)
Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_names=merge(Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo,publ_matchednames,by = "target")
nrow(Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_names)
#1043
colnames(Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_names)

#First, let' look at the groups

Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_names_derup=
  Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_names[!duplicated(Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_names[,c("Row.names","Group")]),]

nrow(Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_names_derup)
#1034
colnames(Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_names_derup)

Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_names_derup[,c("OTU_sum","Group")]%>%group_by(Group)%>%summarise_all(~sum(.))

Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_names_derup_sum=
  Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_names_derup[,c(3:45,60)]%>%group_by(Group)%>%summarise_all(~sum(.))




Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_names_derup_sum_M=melt(Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_names_derup_sum)

#Merge with the metadata
Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_names_derup_sum_M_trt=merge(Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_names_derup_sum_M,
                                                                       sample_data(Mar_leaf.fung_decon_pr_petr_S.rar),by.x = "variable",
                                                                            by.y = "row.names")


Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_names_derup_sum_M_trt2=merge(Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_names_derup_sum_M_trt,
                                                                             sample_sums(Mar_leaf.fung_decon_pr_petr_S.rar),by.x = "variable",
                                                                             by.y = "row.names")

Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_names_derup_sum_M_trt2$prop_abun=Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_names_derup_sum_M_trt2$value/
  Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_names_derup_sum_M_trt2$y


#Second, let' look at the Interactions

Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_names_derup_int=
  Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_names[!duplicated(Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_names[,c("Row.names","Interaction")]),]
unique(Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_names_derup_int$Interaction)
nrow(Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_names_derup_int)
#1034
colnames(Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_names_derup_int)

Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_names_derup_int_sum=
  Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_names_derup_int[,c(3:45,61)]%>%group_by(Interaction)%>%summarise_all(~sum(.))


unique(Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_names_derup_int_sum$Interaction)

Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_names_derup_int_sum_M=melt(Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_names_derup_int_sum)

#Merge with the metadata
Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_names_derup_int_sum_trt=merge(Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_names_derup_int_sum_M,
                                                                              sample_data(Mar_leaf.fung_decon_pr_petr_S.rar),by.x = "variable",
                                                                              by.y = "row.names")


Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_names_derup_int_sum_trt2=merge(Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_names_derup_int_sum_trt,
                                                                               sample_sums(Mar_leaf.fung_decon_pr_petr_S.rar),by.x = "variable",
                                                                               by.y = "row.names")

Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_names_derup_int_sum_trt2$prop_abun=Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_names_derup_int_sum_trt2$value/
  Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_names_derup_int_sum_trt2$y


#Let's combine the interaction groups with the indicators
colnames(Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_names_derup_int_sum_trt2)[colnames(Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_names_derup_int_sum_trt2)=="Interaction"]="grp_inter"
colnames(Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_names_derup_sum_M_trt2)[colnames(Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_names_derup_sum_M_trt2)=="Group"]="grp_inter"

Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_comb=rbind(subset(Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_names_derup_int_sum_trt2,
                                                                  grp_inter!="Unknown"),
                                                           subset(Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_names_derup_sum_M_trt2,
                                                                  grp_inter=="Cave-in-Rock Indicator"|grp_inter=="Prairie Indicator"))

Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_comb_not_match=subset(Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_comb,grp_inter=="Not_matched")
Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_comb_not_match$grp_inter=rep("Endophyte")
Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_comb_not_match$value=1000-Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_comb_not_match$value
Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_comb_not_match$prop_abun=1-Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_comb_not_match$prop_abun
Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_comb=rbind(Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_comb_not_match,
                                                           Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_comb)
unique(Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_comb$grp_inter)
colnames(Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_comb)
unique(with(Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_comb,interaction(collect_date,plant_type)))
Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_comb$round=ifelse(Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_comb$collect_date=="7/21/2018"|
                                                               Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_comb$collect_date=="8/14/2018",
                                                                           "Round1", ifelse(Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_comb$collect_date=="9/18/2018"|
                                                                                              Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_comb$collect_date=="8/21/2018",
                                                                                            "Round2","Seed"))

Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_comb_E=subset(Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_comb,grp_inter!="Not_matched")
Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_comb_sum=Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_comb_E%>%group_by(round,rain_trt,grp_inter)%>%
  summarise_at("prop_abun",list(~mean(.),se=~sd(.)/sqrt(n())))

#####PLOT: Figure S5####


grp_inter_order=c("Endophyte","Pathogen","Context Mutualist","Mutualist","Prairie Indicator","Cave-in-Rock Indicator")
ggplot(subset(Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_comb_E,rain_trt!="Nano"))+geom_point(aes(x=factor(interaction(round,rain_trt), levels = 
                                                                 c("Seed.No_rain","Round1.Sterile_rain","Round1.Live_Rain","Round1.Rain",
                                                                   "Round2.Sterile_rain","Round2.Live_Rain","Round2.Rain")),
                                                      y=prop_abun,shape=rain_trt,
                                                      fill=factor(interaction(round,rain_trt),levels = 
                                                                    c("Seed.No_rain","Round1.Sterile_rain","Round1.Live_Rain","Round1.Rain",
                                                                      "Round2.Sterile_rain","Round2.Live_Rain","Round2.Rain")),
                                                      color=factor(interaction(round,rain_trt),levels = 
                                                                     c("Seed.No_rain","Round1.Sterile_rain","Round1.Live_Rain","Round1.Rain",
                                                                       "Round2.Sterile_rain","Round2.Live_Rain","Round2.Rain"))),size=3,alpha=0.5)+
  geom_errorbar(data=subset(Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_comb_sum,rain_trt!="Nano"),aes(x=factor(interaction(round,rain_trt),levels = 
                                                                                                                 c("Seed.No_rain","Round1.Sterile_rain","Round1.Live_Rain", 
                                                                                                                   "Round1.Rain","Round2.Sterile_rain","Round2.Live_Rain",
                                                                                                                   "Round2.Rain")),
                                                                                                      y=mean,ymax=mean+se,ymin=mean-se),width=0.5)+
  geom_point(data=subset(Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_comb_sum,rain_trt!="Nano"),aes(x=factor(interaction(round,rain_trt),
                                                                                   levels = c("Seed.No_rain","Round1.Sterile_rain","Round1.Live_Rain","Round1.Rain",
                                                                                              "Round2.Sterile_rain","Round2.Live_Rain","Round2.Rain")),
                                                                                   y=mean,shape=rain_trt,
                                                         fill=factor(interaction(round,rain_trt)),
                                                         color=factor(interaction(round,rain_trt))),size=10, stroke=1.5)+
  facet_wrap(vars(factor(grp_inter,levels = grp_inter_order)),scales="free_y")+
  scale_shape_manual(values = c(24,13,21,25))+ylab("Proportion of read\nper sample")+
  scale_color_manual(values = c("#E31A1C","black", "black","black","black","black","black","black"),name=NULL)+
  scale_fill_manual(values = c("#E31A1C","#6A3D9A","#1F78B4","black","#6A3D9A","#1F78B4","black"),name=NULL)+
  scale_x_discrete(labels=c("Seed","Sterile","Live","Rain","Sterile","Live","Rain"),name=NULL)+theme_cowplot()+
  theme(axis.title = element_text(size = 32),axis.text.y = element_text(size = 26),axis.text.x = element_text(size = 22),
        legend.title = element_blank(),strip.text = element_text(size = 28),legend.position = "none",strip.background = element_rect(fill = NA,colour = "black"))

#2100*900
#Petri_pub_seq_grp_prop_dot_nan





#####Stats Petri Published leaf endophyte#####

Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_comb_E_nan = subset(Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_comb_E,rain_trt!="Nano")
head(Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_comb_E_nan)

Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_comb_E_nan$rain_trt_round_nan=with(Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_comb_E_nan,
                                                                          interaction(round,rain_trt))

unique(Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_comb_E_nan$grp_inter)

#"Endophyte"

hist(subset(Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_comb_E_nan,
            grp_inter=="Endophyte")$prop_abun)

End_pub_seq_mod=lmer(sqrt(prop_abun)~rain_trt_round_nan+(1|gh_block),subset(Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_comb_E_nan,
                                                                                   grp_inter=="Endophyte"))

hist(resid(End_pub_seq_mod))
qqPlot(resid(End_pub_seq_mod))
shapiro.test(resid(End_pub_seq_mod))
#W = 0.96464, p-value = 0.3134

anova(End_pub_seq_mod)
#rain_trt_round_nan 0.24507 0.040844     6    28  0.7583 0.6084
#no change

#emmeans(End_pub_seq_mod, pairwise~trt_date,adjust="fdr")


#Context Mutualist
hist(subset(Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_comb_E_nan,
            grp_inter=="Context Mutualist")$prop_abun)

Con_Mut_pub_seq_mod=lmer(log(prop_abun+0.0001)~rain_trt_round_nan+(1|gh_block),subset(Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_comb_E_nan,
                                                                                       grp_inter=="Context Mutualist"))

hist(resid(Con_Mut_pub_seq_mod))
qqPlot(resid(Con_Mut_pub_seq_mod))
shapiro.test(resid(Con_Mut_pub_seq_mod))
plot(Con_Mut_pub_seq_mod)
#W = 0.75566, p-value = 3.079e-06

anova(Con_Mut_pub_seq_mod)
#interaction(round, rain_trt) 66.878  11.146     6    28  3.5498 0.009686 **

emmeans(Con_Mut_pub_seq_mod, pairwise~rain_trt_round_nan,adjust="fdr")


#There is an extreme outlier 

#Context Mutualist

hist(subset(Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_comb_E_nan,
            grp_inter=="Context Mutualist"&prop_abun<0.2)$prop_abun)

Con_Mut_pub_seq_mod_out=lmer(log(prop_abun+0.001)~rain_trt_round_nan+(1|gh_block),subset(Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_comb_E_nan,
                                                                                    grp_inter=="Context Mutualist"&prop_abun<0.2))

hist(resid(Con_Mut_pub_seq_mod_out))
qqPlot(resid(Con_Mut_pub_seq_mod_out))
shapiro.test(resid(Con_Mut_pub_seq_mod_out))
plot(Con_Mut_pub_seq_mod)
#W = 0.81129, p-value = 4.27e-05

anova(Con_Mut_pub_seq_mod_out)
#rain_trt_round_nan  9.102   1.517     6    27  4.8191 0.001842 **

emmeans(Con_Mut_pub_seq_mod_out, pairwise~rain_trt_round_nan,adjust="fdr")
#Round1.Live_Rain - Round1.Rain             -1.4256 0.408 26.5 -3.498  0.0127 
#Round1.Live_Rain - Round2.Rain             -1.3702 0.367 25.4 -3.733  0.0127 
#Round2.Live_Rain - Round1.Rain             -1.3100 0.397 24.3 -3.302  0.0127 
#Round2.Live_Rain - Round2.Rain             -1.2546 0.368 25.4 -3.410  0.0127
#Seed.No_rain - Round1.Rain                 -1.4256 0.433 18.6 -3.289  0.0138 
#Seed.No_rain - Round2.Rain                 -1.3702 0.397 17.0 -3.454  0.0127 
#Round1.Rain - Round1.Sterile_rain           1.0594 0.408 26.5  2.599  0.0396
#Round2.Rain - Round1.Sterile_rain           1.0040 0.367 25.4  2.735  0.0336 

ggplot(subset(Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_comb_E_nan,
              grp_inter=="Context Mutualist"&prop_abun<0.2),aes(x=rain_trt_round_nan,y=prop_abun))+geom_boxplot()


#Zero inflated, let's try negative binomial with the total reads
Con_Mut_pub_seq_mod_out_nb=glmer.nb((value)~rain_trt_round_nan+(1|gh_block),subset(Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_comb_E_nan,
                                                                                     grp_inter=="Context Mutualist"&prop_abun<0.2))

hist(resid(Con_Mut_pub_seq_mod_out_nb))
qqPlot(resid(Con_Mut_pub_seq_mod_out_nb))
shapiro.test(resid(Con_Mut_pub_seq_mod_out_nb))
plot(Con_Mut_pub_seq_mod_out_nb)
#W = 0.81129, p-value = 4.27e-05

Anova(Con_Mut_pub_seq_mod_out_nb,type = 3)

#Occurrences to sparse and too variable to runs stats



#Mutualist
hist(subset(Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_comb_E_nan,
            grp_inter=="Mutualist")$prop_abun)

Mut_pub_seq_mod=lmer((prop_abun+0.001)^-1~rain_trt_round_nan+(1|gh_block),subset(Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_comb_E_nan,
                                                                    grp_inter=="Mutualist"))

hist(resid(Mut_pub_seq_mod))
qqPlot(resid(Mut_pub_seq_mod))
shapiro.test(resid(Mut_pub_seq_mod))
plot(Mut_pub_seq_mod)
#W = 0.88443, p-value = 0.001549

anova(Mut_pub_seq_mod)


#Zero inflated, let's try negative binomial with the total reads
Mut_pub_seq_mod_nb=glmer.nb((value)~rain_trt_round_nan+(1|gh_block),subset(Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_comb_E_nan,
                                                                                    grp_inter=="Mutualist"))

hist(resid(Mut_pub_seq_mod_nb))
qqPlot(resid(Mut_pub_seq_mod_nb))
shapiro.test(resid(Mut_pub_seq_mod_nb))
plot(Mut_pub_seq_mod_nb)
#W = 0.91433, p-value = 0.009822

Anova(Mut_pub_seq_mod_nb,type = 3)
#

emmeans(Mut_pub_seq_mod_nb, pairwise~rain_trt_round_nan,adjust="fdr")
#nada
#Occurrences to sparse and too variable to runs stats


#Pathogen
hist(subset(Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_comb_E_nan,
            grp_inter=="Pathogen")$prop_abun)
Path_pub_seq_mod=lmer(logit(prop_abun)~rain_trt_round_nan+(1|gh_block),subset(Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_comb_E_nan,
                                                                                grp_inter=="Pathogen"))

hist(resid(Path_pub_seq_mod))
qqPlot(resid(Path_pub_seq_mod))
shapiro.test(resid(Path_pub_seq_mod))
plot(Path_pub_seq_mod)
#W = 0.93539, p-value = 0.04058

anova(Path_pub_seq_mod)
#rain_trt_round_nan 41.348  6.8913     6    28   1.707 0.1563
#no change

#####POSTHOC TEST: Figure S5b####
emmeans(Path_pub_seq_mod, pairwise~rain_trt_round_nan,adjust="fdr")


#Zero inflated, let's try negative binomial with the total reads

Path_pub_seq_mod_nb=glmer.nb((value)~rain_trt_round_nan+(1|gh_block),subset(Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_comb_E_nan,
                                                                          grp_inter=="Pathogen"))

hist(resid(Path_pub_seq_mod_nb))
qqPlot(resid(Path_pub_seq_mod_nb))
shapiro.test(resid(Path_pub_seq_mod_nb))
plot(Path_pub_seq_mod_nb)
#W = 0.95475, p-value = 0.1585


Anova(Path_pub_seq_mod_nb,type = 3)
#rain_trt_round_nan  27.215  6   0.000132 ***
emmeans(Path_pub_seq_mod_nb, pairwise~rain_trt_round_nan,adjust="fdr")

#Prairie Indicator


hist(subset(Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_comb_E_nan,
            grp_inter=="Prairie Indicator")$prop_abun)

Prair_pub_seq_mod=lmer(log(prop_abun+0.001)~rain_trt_round_nan+(1|gh_block),subset(Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_comb_E_nan,
                                                                            grp_inter=="Prairie Indicator"))

hist(resid(Prair_pub_seq_mod))
qqPlot(resid(Prair_pub_seq_mod))
shapiro.test(resid(Prair_pub_seq_mod))
plot(Prair_pub_seq_mod)
#W = 0.87764, p-value = 0.08591

anova(Prair_pub_seq_mod)
#rain_trt_round_nan 10.139  1.6898     6    28  0.7161 0.6398
#No change
#emmeans(Path_pub_seq_mod, pairwise~rain_trt_round_nan,adjust="fdr")

#Zero inflated, let's try negative binomial with the total reads

Prair_pub_seq_mod_nb=glmer.nb((value)~rain_trt_round_nan+(1|gh_block),subset(Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_comb_E_nan,
                                                                               grp_inter=="Prairie Indicator"))

hist(resid(Prair_pub_seq_mod_nb))
qqPlot(resid(Prair_pub_seq_mod_nb))
shapiro.test(resid(Prair_pub_seq_mod_nb))
plot(Prair_pub_seq_mod_nb)
#W = 0.87764, p-value = 0.08591

Anova(Prair_pub_seq_mod_nb,type = 3)
#rain_trt_round_nan 5.5458  6    0.47593
 
#Cave-in-Rock Indicator

hist(subset(Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_comb_E_nan,
            grp_inter=="Cave-in-Rock Indicator"&prop_abun<0.003)$prop_abun)

Cave_pub_seq_mod=lmer(log(prop_abun+0.001)~rain_trt_round_nan+(1|gh_block),subset(Mar_leaf.fung_decon_pr_petr_S.rar_pub_Endo_comb_E_nan,
                                                                            grp_inter=="Cave-in-Rock Indicator"))

hist(resid(Cave_pub_seq_mod))
qqPlot(resid(Cave_pub_seq_mod))
shapiro.test(resid(Cave_pub_seq_mod))
plot(Cave_pub_seq_mod)
#W = 0.75167, p-value = 2.624e-06

anova(Cave_pub_seq_mod)
#rain_trt_round_nan 0.80793 0.13465     6    28  1.0111 0.4381

emmeans(Cave_pub_seq_mod, pairwise~rain_trt_round_nan,adjust="fdr")


#####Bray Petri Pariwise####
Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt=read.csv(here::here("R_file","Pairwise_turnover_grad_distance_fung_v4.2.2020_decon_pruned_ZOTU_full_rar_Marshall.csv"))
head(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt)
nrow(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt)
#6786


#Subset to only petri dishes between rounds within trt

Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_dishes=subset(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt,
                                                                     s1_s2_plant_type=="Petri.Petri")
nrow(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_dishes)
#496
unique(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_dishes$s1_s2_collect_date)

Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_dishes_rounds=subset(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_dishes,
                                                                            s1_s2_collect_date=="9/18/2018.8/14/2018"|
                                                                              s1_s2_collect_date=="8/14/2018.9/18/2018")
nrow(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_dishes_rounds)
#256

Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_dishes_rounds_w_trt=subset(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_dishes_rounds,
                                                                                  s1_s2_rain_trt=="Live_Rain.Live_Rain"|s1_s2_rain_trt=="Nano.Nano"|
                                                                                    s1_s2_rain_trt=="Sterile_rain.Sterile_rain")

unique(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_dishes_rounds_w_trt$s1_s2_rain_trt)

Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_dishes_rounds_w_trt$treatment=rep("Between_Rounds")


#Subset Seed to trt and round

unique(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt$s1_s2_plant_type)

Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_seed=subset(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt, 
                                                                   s1_s2_plant_type=="Petri.Seed"|s1_s2_plant_type=="Seed.Petri")
nrow(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_seed)
#128

unique(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_seed$s1_s2_rain_trt)


Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_seed$treatment=rep("Seed")


#Subset Rain to trt and round

unique(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt$s1_s2_plant_type)

Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_rain0=subset(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt, 
                                                                   s1_s2_plant_type=="Petri.Rain"|s1_s2_plant_type=="Rain.Petri")


nrow(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_rain0)
#1792

unique(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_rain0$s1_s2_collect_date)

#Petri 
#start_date collect_dates 
#8/21/2018	9/18/2018
#7/21/2018	8/14/2018

Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_rain=subset(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_rain0,s1_s2_collect_date=="8/21/2018.9/18/2018"|
                                                                  s1_s2_collect_date=="9/18/2018.8/21/2018"|s1_s2_collect_date=="7/21/2018.8/14/2018"|
                                                                  s1_s2_collect_date=="8/14/2018.7/21/2018")


unique(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_rain$s1_s2_collect_date)
nrow(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_rain)
#112

Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_rain$treatment=rep("Rain")


#combine the matrix 
Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_pairwise=rbind(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_dishes_rounds_w_trt,
                                                                      Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_seed,
                                                                      Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_rain)


Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_pairwise$leaf_sample=ifelse(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_pairwise$
                                                                                     s1_plant_type=="Petri",as.character(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_pairwise$
                                                                                                                           sample1),as.character(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_pairwise$
                                                                                                                                                   sample2))
unique(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_pairwise$leaf_sample)
Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_pairwise$unique_rain_trt=ifelse(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_pairwise$
                                                                                         s1_plant_type=="Petri",as.character(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_pairwise$
                                                                                                                               s1_rain_trt),as.character(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_pairwise$
                                                                                                                                                           s2_rain_trt))

unique(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_pairwise$unique_rain_trt)

Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_pairwise$round=ifelse(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_pairwise$
                                                                               s1_collect_date=="9/18/2018"|Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_pairwise$
                                                                               s2_collect_date=="9/18/2018","Round2","Round1")



head(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_pairwise)
nrow(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_pairwise)
#328

Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_pairwise$gh_block_comb=ifelse(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_pairwise$
                                                                            s1_plant_type=="Petri",as.character(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_pairwise$
                                                                                                                  s1_gh_block),
                                                                          as.character(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_pairwise$
                                                                                         s2_gh_block))

head(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_pairwise)
nrow(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_pairwise)
#328
unique(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_pairwise$gh_block_comb)

#Drop the between the rounds analyses
Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_pair_sub=subset(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_pairwise,treatment!="Between_Rounds")
nrow(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_pair_sub)
#240


#Bray

Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_pair_sub_sum=Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_pair_sub%>%
  group_by(gh_block_comb,leaf_sample,round,unique_rain_trt,treatment)%>%
  summarise_at("bray",list(~mean(.),se=~sd(.)/sqrt(n()),~n()))


Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_pair_sub_sum1=data.frame(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_pair_sub_sum
                                                                                [,c("round","unique_rain_trt","treatment","mean")])

colnames(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_pair_sub_sum1)=c("round","unique_rain_trt","treatment","bray")
Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_pair_sub_sum_sum=Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_pair_sub_sum1%>%
  group_by(round,unique_rain_trt,treatment)%>%
  summarise_all(list(~mean(.),se=~sd(.)/sqrt(n()),~n(),~sd(.)))

#No nanopure samples
Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_pair_sub_sum_nan=subset(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_pair_sub_sum, unique_rain_trt!="Nano")
Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_pair_sub_sum_sum_nan=subset(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_pair_sub_sum_sum, unique_rain_trt!="Nano")

Rain_order_nan=c("Sterile_rain","Live_Rain")
rain_trt_round_nan=c("Sterile_rain.Round1","Live_Rain.Round1","Sterile_rain.Round2","Live_Rain.Round2")




(bray1000_petri_nan_p=ggplot(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_pair_sub_sum_nan,aes(x=factor(treatment,levels = c("Seed","Rain")),y=mean))+
    geom_point(aes(shape=factor(unique_rain_trt,levels = Rain_order_nan),fill=factor(unique_rain_trt,levels = Rain_order_nan),
                   group=factor(interaction(unique_rain_trt,round),levels = rain_trt_round_nan)), 
               color="black",size=4,position = position_dodge(0.65),alpha=.15)+
    geom_errorbar(data=Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_pair_sub_sum_sum_nan, 
                  aes(x=factor(treatment,levels = c("Seed","Rain")),group=factor(interaction(unique_rain_trt,round),levels = rain_trt_round_nan),
                      y=mean,ymin=mean-se,ymax=mean+se),width=.5,
                  color="black",position = position_dodge(0.65))+
    geom_point(data=Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_pair_sub_sum_sum_nan,
               aes(shape=factor(unique_rain_trt,levels = Rain_order_nan),group=factor(interaction(unique_rain_trt,round),levels = rain_trt_round_nan),
                   fill=factor(unique_rain_trt,levels = Rain_order_nan),x=factor(treatment,levels = c("Seed","Rain"))), size=10,
               color="black",position = position_dodge(0.65))+facet_nested(.~round)+
    scale_shape_manual(values = c(25,24),name=NULL)+scale_x_discrete(labels=c("Seed","Rain"))+
    scale_fill_manual(values = c("#6A3D9A","#1F78B4"),name=NULL)+
    ylab("Bray-Curtis distance")+theme_classic()+
    theme(axis.title.x = element_blank(),axis.title.y = element_text(size = 28),axis.text = element_text(size = 30),panel.border =  element_rect(fill = NA,linetype="solid"),
          legend.position = "none",strip.text = element_text(size = 36),strip.placement = "outside",strip.background = element_rect(linetype="blank"),
          panel.spacing = unit(0, "lines")))


#distance to seeds

subset(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_pair_sub_sum_sum_nan,treatment=="Seed"&unique_rain_trt=="Live_Rain")
subset(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_pair_sub_sum_sum_nan,treatment=="Seed"&unique_rain_trt=="Sterile_rain")
#round1
(0.886-0.770)/0.770
#round2
(0.846-0.782)/0.782

(((0.886-0.770)/0.770)+((0.846-0.782)/0.782))/2

#distance to rain

subset(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_pair_sub_sum_sum_nan,treatment=="Rain"&unique_rain_trt=="Live_Rain")
subset(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_pair_sub_sum_sum_nan,treatment=="Rain"&unique_rain_trt=="Sterile_rain")
#round1
(0.890-0.866)/0.866
#round2
(0.709-0.845)/0.845

(((0.890-0.866)/0.866)+((0.709-0.845)/0.845))/2

Petri_pair_dist_1000_nan_mod=lmer(sqrt(mean)~round*unique_rain_trt*treatment+(1|leaf_sample)+(1|gh_block_comb),data=Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_pair_sub_sum_nan)
plot(Petri_pair_dist_1000_nan_mod)
hist(resid(Petri_pair_dist_1000_nan_mod))
qqPlot(resid(Petri_pair_dist_1000_nan_mod))
shapiro.test(resid(Petri_pair_dist_1000_nan_mod))
#W = 0.95926, p-value = 0.09421

#####Table S3 Bray-Curtis####
anova(Petri_pair_dist_1000_nan_mod)
#round                           0.0029206 0.0029206     1 19.366  4.3619  0.050183 .
#round:unique_rain_trt           0.0025231 0.0025231     1 17.109  3.7682  0.068891 .  
#round:treatment                 0.0073274 0.0073274     1 20.000 10.9433  0.003511 ** 
#unique_rain_trt:treatment       0.0207200 0.0207200     1 20.000 30.9448 1.912e-05 ***
#round:unique_rain_trt:treatment 0.0028535 0.0028535     1 20.000  4.2617  0.052194 .  
#CHANGED

emmeans(Petri_pair_dist_1000_nan_mod,pairwise~unique_rain_trt*treatment*round,adjust="fdr")

emmeans(Petri_pair_dist_1000_nan_mod,pairwise~unique_rain_trt|treatment)
emmeans(Petri_pair_dist_1000_nan_mod,pairwise~treatment|round)
emmeans(Petri_pair_dist_1000_nan_mod,pairwise~unique_rain_trt|treatment|round)

#####POSTHOC TEST: Figure 2a####
emmeans(Petri_pair_dist_1000_nan_mod,pairwise~unique_rain_trt*treatment|round,adjust="fdr")

#####Overlap in Taxa between Petri Taxa####


Mar_leaf_rain.fung_decon_rar_Mar=subset_samples(Mar_leaf_rain.fung_decon_rar, sub_proj=="Marshall")
Mar_leaf_rain.fung_decon_rar_Mar_petri=subset_samples(Mar_leaf_rain.fung_decon_rar_Mar, plant_type!="Adult"&plant_type!="Seedling")

Mar_leaf.fung_decon_pr_petr_S.rar=subset_samples(Mar_leaf_rain.fung_decon_rar_Mar_petri,plant_type!="Rain"|collect_date=="8/21/2018"|
                                                   collect_date=="7/21/2018")

nsamples(Mar_leaf.fung_decon_pr_petr_S.rar)
#43

unique(sample_data(Mar_leaf.fung_decon_pr_petr_S.rar)$sub_proj)

Mar_leaf.fung_decon_pr_petr_S.rar=prune_taxa(taxa_sums(Mar_leaf.fung_decon_pr_petr_S.rar) > 0, Mar_leaf.fung_decon_pr_petr_S.rar)
ntaxa(Mar_leaf.fung_decon_pr_petr_S.rar)
#1027
sum(taxa_sums(Mar_leaf.fung_decon_pr_petr_S.rar))
#43000
min(sample_sums(Mar_leaf.fung_decon_pr_petr_S.rar))
#1000
max(sample_sums(Mar_leaf.fung_decon_pr_petr_S.rar))
#1000

ntaxa(Mar_leaf.fung_decon_pr_petr_S.rar)
#1027
nsamples(Mar_leaf.fung_decon_pr_petr_S.rar)
#43

#what OTUs are in Seeds 
Mar_leaf.fung_decon_pr_petr_S.rar_SEED=subset_samples(Mar_leaf.fung_decon_pr_petr_S.rar, plant_type=="Seed")
Mar_leaf.fung_decon_pr_petr_S.rar_SEED=prune_taxa(taxa_sums(Mar_leaf.fung_decon_pr_petr_S.rar_SEED)>0,Mar_leaf.fung_decon_pr_petr_S.rar_SEED)

nsamples(Mar_leaf.fung_decon_pr_petr_S.rar_SEED)
#4
ntaxa(Mar_leaf.fung_decon_pr_petr_S.rar_SEED)
#152

Petri_SEED=taxa_names(Mar_leaf.fung_decon_pr_petr_S.rar_SEED)

unique(sample_data(Mar_leaf.fung_decon_pr_petr_S.rar)$collect_date)

#Round 1

#what OTUs are in Rain  
Mar_leaf.fung_decon_pr_petr_S.rar_RAIN_R1=subset_samples(Mar_leaf.fung_decon_pr_petr_S.rar, plant_type=="Rain"&collect_date=="7/21/2018")
Mar_leaf.fung_decon_pr_petr_S.rar_RAIN_R1=prune_taxa(taxa_sums(Mar_leaf.fung_decon_pr_petr_S.rar_RAIN_R1)>0,Mar_leaf.fung_decon_pr_petr_S.rar_RAIN_R1)

nsamples(Mar_leaf.fung_decon_pr_petr_S.rar_RAIN_R1)
#3
ntaxa(Mar_leaf.fung_decon_pr_petr_S.rar_RAIN_R1)
#417

Petri_RAIN_R1=taxa_names(Mar_leaf.fung_decon_pr_petr_S.rar_RAIN_R1)

unique(sample_data(Mar_leaf.fung_decon_pr_petr_S.rar)$collect_date)
unique(sample_data(Mar_leaf.fung_decon_pr_petr_S.rar)$rain_trt)
unique(with(sample_data(Mar_leaf.fung_decon_pr_petr_S.rar),interaction(collect_date,rain_trt)))
#what OTUs are in Nano
Mar_leaf.fung_decon_pr_petr_S.rar_Nano_R1=subset_samples(Mar_leaf.fung_decon_pr_petr_S.rar, rain_trt=="Nano"&collect_date=="8/14/2018")
Mar_leaf.fung_decon_pr_petr_S.rar_Nano_R1=prune_taxa(taxa_sums(Mar_leaf.fung_decon_pr_petr_S.rar_Nano_R1)>0,Mar_leaf.fung_decon_pr_petr_S.rar_Nano_R1)

nsamples(Mar_leaf.fung_decon_pr_petr_S.rar_Nano_R1)
#4
ntaxa(Mar_leaf.fung_decon_pr_petr_S.rar_Nano_R1)
#142

Mar_leaf.fung_decon_pr_petr_S.rar_NANO_R1_overlap=data.frame(estimate_richness(Mar_leaf.fung_decon_pr_petr_S.rar_Nano_R1,measures ="Observed"),
                                                             estimate_richness(prune_taxa(Petri_SEED,Mar_leaf.fung_decon_pr_petr_S.rar_Nano_R1),
                                                                               measures ="Observed"),
                                                             estimate_richness(prune_taxa(Petri_RAIN_R1,Mar_leaf.fung_decon_pr_petr_S.rar_Nano_R1),
                                                                               measures ="Observed"),
                                                             "overlap_seed_sum"=sample_sums(prune_taxa(Petri_SEED,Mar_leaf.fung_decon_pr_petr_S.rar_Nano_R1)),
                                                             "overlap_rain_sum"=sample_sums(prune_taxa(Petri_RAIN_R1,Mar_leaf.fung_decon_pr_petr_S.rar_Nano_R1)),
                                                             sample_data(Mar_leaf.fung_decon_pr_petr_S.rar_Nano_R1))
colnames(Mar_leaf.fung_decon_pr_petr_S.rar_NANO_R1_overlap)[1:3]=c("tot_rich","overlap_seed_rich","overlap_rain_rich")

#what OTUs are in Sterile
Mar_leaf.fung_decon_pr_petr_S.rar_Sterile_R1=subset_samples(Mar_leaf.fung_decon_pr_petr_S.rar, rain_trt=="Sterile_rain"&collect_date=="8/14/2018")
Mar_leaf.fung_decon_pr_petr_S.rar_Sterile_R1=prune_taxa(taxa_sums(Mar_leaf.fung_decon_pr_petr_S.rar_Sterile_R1)>0,Mar_leaf.fung_decon_pr_petr_S.rar_Sterile_R1)

nsamples(Mar_leaf.fung_decon_pr_petr_S.rar_Sterile_R1)
#6
ntaxa(Mar_leaf.fung_decon_pr_petr_S.rar_Sterile_R1)
#163


Mar_leaf.fung_decon_pr_petr_S.rar_STERILE_R1_overlap=data.frame(estimate_richness(Mar_leaf.fung_decon_pr_petr_S.rar_Sterile_R1,measures ="Observed"),
                                                                       estimate_richness(prune_taxa(Petri_SEED,Mar_leaf.fung_decon_pr_petr_S.rar_Sterile_R1),
                                                                                         measures ="Observed"),
                                                                       estimate_richness(prune_taxa(Petri_RAIN_R1,Mar_leaf.fung_decon_pr_petr_S.rar_Sterile_R1),
                                                                                         measures ="Observed"),
                                                                       "overlap_seed_sum"=sample_sums(prune_taxa(Petri_SEED,Mar_leaf.fung_decon_pr_petr_S.rar_Sterile_R1)),
                                                                       "overlap_rain_sum"=sample_sums(prune_taxa(Petri_RAIN_R1,Mar_leaf.fung_decon_pr_petr_S.rar_Sterile_R1)),
                                                                       sample_data(Mar_leaf.fung_decon_pr_petr_S.rar_Sterile_R1))
colnames(Mar_leaf.fung_decon_pr_petr_S.rar_STERILE_R1_overlap)[1:3]=c("tot_rich","overlap_seed_rich","overlap_rain_rich")


#what OTUs are in Live
Mar_leaf.fung_decon_pr_petr_S.rar_Live_R1=subset_samples(Mar_leaf.fung_decon_pr_petr_S.rar, rain_trt=="Live_Rain"&collect_date=="8/14/2018")
Mar_leaf.fung_decon_pr_petr_S.rar_Live_R1=prune_taxa(taxa_sums(Mar_leaf.fung_decon_pr_petr_S.rar_Live_R1)>0,Mar_leaf.fung_decon_pr_petr_S.rar_Live_R1)

nsamples(Mar_leaf.fung_decon_pr_petr_S.rar_Live_R1)
#6
ntaxa(Mar_leaf.fung_decon_pr_petr_S.rar_Live_R1)
#154


Mar_leaf.fung_decon_pr_petr_S.rar_LIVE_R1_overlap=data.frame(estimate_richness(Mar_leaf.fung_decon_pr_petr_S.rar_Live_R1,measures ="Observed"),
                                                                estimate_richness(prune_taxa(Petri_SEED,Mar_leaf.fung_decon_pr_petr_S.rar_Live_R1),
                                                                                  measures ="Observed"),
                                                                estimate_richness(prune_taxa(Petri_RAIN_R1,Mar_leaf.fung_decon_pr_petr_S.rar_Live_R1),
                                                                                  measures ="Observed"),
                                                                "overlap_seed_sum"=sample_sums(prune_taxa(Petri_SEED,Mar_leaf.fung_decon_pr_petr_S.rar_Live_R1)),
                                                                "overlap_rain_sum"=sample_sums(prune_taxa(Petri_RAIN_R1,Mar_leaf.fung_decon_pr_petr_S.rar_Live_R1)),
                                                                sample_data(Mar_leaf.fung_decon_pr_petr_S.rar_Live_R1))
colnames(Mar_leaf.fung_decon_pr_petr_S.rar_LIVE_R1_overlap)[1:3]=c("tot_rich","overlap_seed_rich","overlap_rain_rich")



#Round 2

#what OTUs are in Rain  
Mar_leaf.fung_decon_pr_petr_S.rar_RAIN_R2=subset_samples(Mar_leaf.fung_decon_pr_petr_S.rar, plant_type=="Rain"&collect_date=="8/21/2018")
Mar_leaf.fung_decon_pr_petr_S.rar_RAIN_R2=prune_taxa(taxa_sums(Mar_leaf.fung_decon_pr_petr_S.rar_RAIN_R2)>0,Mar_leaf.fung_decon_pr_petr_S.rar_RAIN_R2)

nsamples(Mar_leaf.fung_decon_pr_petr_S.rar_RAIN_R2)
#4
ntaxa(Mar_leaf.fung_decon_pr_petr_S.rar_RAIN_R2)
#439

Petri_RAIN_R2=taxa_names(Mar_leaf.fung_decon_pr_petr_S.rar_RAIN_R2)

#what OTUs are in Nano
Mar_leaf.fung_decon_pr_petr_S.rar_Nano_R2=subset_samples(Mar_leaf.fung_decon_pr_petr_S.rar, rain_trt=="Nano"&collect_date=="9/18/2018")
Mar_leaf.fung_decon_pr_petr_S.rar_Nano_R2=prune_taxa(taxa_sums(Mar_leaf.fung_decon_pr_petr_S.rar_Nano_R2)>0,Mar_leaf.fung_decon_pr_petr_S.rar_Nano_R2)

nsamples(Mar_leaf.fung_decon_pr_petr_S.rar_Nano_R2)
#4
ntaxa(Mar_leaf.fung_decon_pr_petr_S.rar_Nano_R2)
#158

Mar_leaf.fung_decon_pr_petr_S.rar_NANO_R2_overlap=data.frame(estimate_richness(Mar_leaf.fung_decon_pr_petr_S.rar_Nano_R2,measures ="Observed"),
                                                             estimate_richness(prune_taxa(Petri_SEED,Mar_leaf.fung_decon_pr_petr_S.rar_Nano_R2),
                                                                               measures ="Observed"),
                                                             estimate_richness(prune_taxa(Petri_RAIN_R2,Mar_leaf.fung_decon_pr_petr_S.rar_Nano_R2),
                                                                               measures ="Observed"),
                                                             "overlap_seed_sum"=sample_sums(prune_taxa(Petri_SEED,Mar_leaf.fung_decon_pr_petr_S.rar_Nano_R2)),
                                                             "overlap_rain_sum"=sample_sums(prune_taxa(Petri_RAIN_R2,Mar_leaf.fung_decon_pr_petr_S.rar_Nano_R2)),
                                                             sample_data(Mar_leaf.fung_decon_pr_petr_S.rar_Nano_R2))
colnames(Mar_leaf.fung_decon_pr_petr_S.rar_NANO_R2_overlap)[1:3]=c("tot_rich","overlap_seed_rich","overlap_rain_rich")

#what OTUs are in Sterile
Mar_leaf.fung_decon_pr_petr_S.rar_Sterile_R2=subset_samples(Mar_leaf.fung_decon_pr_petr_S.rar, rain_trt=="Sterile_rain"&collect_date=="9/18/2018")
Mar_leaf.fung_decon_pr_petr_S.rar_Sterile_R2=prune_taxa(taxa_sums(Mar_leaf.fung_decon_pr_petr_S.rar_Sterile_R2)>0,Mar_leaf.fung_decon_pr_petr_S.rar_Sterile_R2)

nsamples(Mar_leaf.fung_decon_pr_petr_S.rar_Sterile_R2)
#6
ntaxa(Mar_leaf.fung_decon_pr_petr_S.rar_Sterile_R2)
#164


Mar_leaf.fung_decon_pr_petr_S.rar_STERILE_R2_overlap=data.frame(estimate_richness(Mar_leaf.fung_decon_pr_petr_S.rar_Sterile_R2,measures ="Observed"),
                                                                estimate_richness(prune_taxa(Petri_SEED,Mar_leaf.fung_decon_pr_petr_S.rar_Sterile_R2),
                                                                                  measures ="Observed"),
                                                                estimate_richness(prune_taxa(Petri_RAIN_R2,Mar_leaf.fung_decon_pr_petr_S.rar_Sterile_R2),
                                                                                  measures ="Observed"),
                                                                "overlap_seed_sum"=sample_sums(prune_taxa(Petri_SEED,Mar_leaf.fung_decon_pr_petr_S.rar_Sterile_R2)),
                                                                "overlap_rain_sum"=sample_sums(prune_taxa(Petri_RAIN_R2,Mar_leaf.fung_decon_pr_petr_S.rar_Sterile_R2)),
                                                                sample_data(Mar_leaf.fung_decon_pr_petr_S.rar_Sterile_R2))
colnames(Mar_leaf.fung_decon_pr_petr_S.rar_STERILE_R2_overlap)[1:3]=c("tot_rich","overlap_seed_rich","overlap_rain_rich")


#what OTUs are in Live
Mar_leaf.fung_decon_pr_petr_S.rar_Live_R2=subset_samples(Mar_leaf.fung_decon_pr_petr_S.rar, rain_trt=="Live_Rain"&collect_date=="9/18/2018")
Mar_leaf.fung_decon_pr_petr_S.rar_Live_R2=prune_taxa(taxa_sums(Mar_leaf.fung_decon_pr_petr_S.rar_Live_R2)>0,Mar_leaf.fung_decon_pr_petr_S.rar_Live_R2)

nsamples(Mar_leaf.fung_decon_pr_petr_S.rar_Live_R2)
#6
ntaxa(Mar_leaf.fung_decon_pr_petr_S.rar_Live_R2)
#136


Mar_leaf.fung_decon_pr_petr_S.rar_LIVE_R2_overlap=data.frame(estimate_richness(Mar_leaf.fung_decon_pr_petr_S.rar_Live_R2,measures ="Observed"),
                                                             estimate_richness(prune_taxa(Petri_SEED,Mar_leaf.fung_decon_pr_petr_S.rar_Live_R2),
                                                                               measures ="Observed"),
                                                             estimate_richness(prune_taxa(Petri_RAIN_R2,Mar_leaf.fung_decon_pr_petr_S.rar_Live_R2),
                                                                               measures ="Observed"),
                                                             "overlap_seed_sum"=sample_sums(prune_taxa(Petri_SEED,Mar_leaf.fung_decon_pr_petr_S.rar_Live_R2)),
                                                             "overlap_rain_sum"=sample_sums(prune_taxa(Petri_RAIN_R2,Mar_leaf.fung_decon_pr_petr_S.rar_Live_R2)),
                                                             sample_data(Mar_leaf.fung_decon_pr_petr_S.rar_Live_R2))
colnames(Mar_leaf.fung_decon_pr_petr_S.rar_LIVE_R2_overlap)[1:3]=c("tot_rich","overlap_seed_rich","overlap_rain_rich")
summary(Mar_leaf.fung_decon_pr_petr_S.rar_LIVE_R2_overlap)

Mar_leaf.fung_decon_pr_petr_S.rar_overlap=rbind(Mar_leaf.fung_decon_pr_petr_S.rar_NANO_R1_overlap,
                                                Mar_leaf.fung_decon_pr_petr_S.rar_STERILE_R1_overlap,
                                                Mar_leaf.fung_decon_pr_petr_S.rar_LIVE_R1_overlap,
                                                Mar_leaf.fung_decon_pr_petr_S.rar_NANO_R2_overlap,
                                                Mar_leaf.fung_decon_pr_petr_S.rar_STERILE_R2_overlap,
                                                Mar_leaf.fung_decon_pr_petr_S.rar_LIVE_R2_overlap)



Mar_leaf.fung_decon_pr_petr_S.rar_overlap$overlap_seed_rich_prop=Mar_leaf.fung_decon_pr_petr_S.rar_overlap$overlap_seed_rich/
  Mar_leaf.fung_decon_pr_petr_S.rar_overlap$tot_rich
Mar_leaf.fung_decon_pr_petr_S.rar_overlap$overlap_rain_rich_prop=Mar_leaf.fung_decon_pr_petr_S.rar_overlap$overlap_rain_rich/
  Mar_leaf.fung_decon_pr_petr_S.rar_overlap$tot_rich

Mar_leaf.fung_decon_pr_petr_S.rar_overlap$overlap_seed_sum_prop=Mar_leaf.fung_decon_pr_petr_S.rar_overlap$overlap_seed_sum/1000
Mar_leaf.fung_decon_pr_petr_S.rar_overlap$overlap_rain_sum_prop=Mar_leaf.fung_decon_pr_petr_S.rar_overlap$overlap_rain_sum/1000


#Richness overlap

Mar_leaf.fung_decon_pr_petr_S.rar_overlap_Rich_M=melt(Mar_leaf.fung_decon_pr_petr_S.rar_overlap[c("rain_trt","collect_date","sampleID_bact",
                                                                                                  "overlap_rain_rich_prop","overlap_seed_rich_prop",
                                                                                                  "gh_block")])


#Richness Overlap


Mar_leaf.fung_decon_pr_petr_S.rar_overlap_Rich_sum=Mar_leaf.fung_decon_pr_petr_S.rar_overlap_Rich_M%>%
  group_by(rain_trt,collect_date,variable)%>%
  summarise_at(vars("value"),list(~mean(.),se=~sd(.)/sqrt(n()),~n(),~sd(.)))
Rain_order=c("Nano","Sterile_rain","Live_Rain")
unique(with(Mar_leaf.fung_decon_pr_petr_S.rar_overlap_Rich_sum,interaction(collect_date,variable)))
round_trt_over_order=c("8/14/2018.overlap_seed_rich_prop","9/18/2018.overlap_seed_rich_prop","8/14/2018.overlap_rain_rich_prop","9/18/2018.overlap_rain_rich_prop")
with(Mar_leaf.fung_decon_pr_petr_S.rar_overlap_Rich_sum,interaction(rain_trt,collect_date))
rain_trt_over_round=c("Nano.8/14/2018","Nano.9/18/2018","Sterile_rain.8/14/2018","Sterile_rain.9/18/2018","Live_Rain.8/14/2018","Live_Rain.9/18/2018")
brewer.pal(n = 12, name = "Paired")



#No Nanopure samples
Mar_leaf.fung_decon_pr_petr_S.rar_overlap_Rich_M_nan=subset(Mar_leaf.fung_decon_pr_petr_S.rar_overlap_Rich_M,rain_trt!="Nano")
Mar_leaf.fung_decon_pr_petr_S.rar_overlap_Rich_sum_nan=subset(Mar_leaf.fung_decon_pr_petr_S.rar_overlap_Rich_sum,rain_trt!="Nano")
Rain_order_nan=c("Sterile_rain","Live_Rain")






(rich_over_petri_nan_p=ggplot(Mar_leaf.fung_decon_pr_petr_S.rar_overlap_Rich_M_nan,
                              aes(x=factor(variable,levels = c("overlap_seed_rich_prop","overlap_rain_rich_prop")),y=value))+
    geom_point(aes(shape=factor(rain_trt,levels = Rain_order_nan),fill=factor(rain_trt,levels = Rain_order_nan),
                   group=factor(rain_trt,levels = Rain_order_nan)), 
               color="black",size=4,position = position_dodge(0.65),alpha=.15)+
    geom_errorbar(data=Mar_leaf.fung_decon_pr_petr_S.rar_overlap_Rich_sum_nan, 
                  aes(x=factor(variable,levels = c("overlap_seed_rich_prop","overlap_rain_rich_prop")),group=factor(rain_trt,levels = Rain_order_nan),
                      y=mean,ymin=mean-se,ymax=mean+se),width=.5,
                  color="black",position = position_dodge(0.65))+
    geom_point(data=Mar_leaf.fung_decon_pr_petr_S.rar_overlap_Rich_sum_nan,
               aes(y=mean,shape=factor(rain_trt,levels = Rain_order_nan),group=factor(rain_trt,levels = Rain_order_nan),
                   fill=factor(rain_trt,levels = Rain_order_nan),x=factor(variable,levels = c("overlap_seed_rich_prop","overlap_rain_rich_prop"))), size=10,
               color="black",position = position_dodge(0.65))+facet_nested(.~(factor(collect_date,labels = c("Round1","Round2"))))+
    scale_shape_manual(values = c(25,24),name=NULL)+scale_x_discrete(labels=c("Seed","Rain"))+
    scale_fill_manual(values = c("#6A3D9A","#1F78B4"),name=NULL)+
    ylab("Proportion ASV overlap")+theme_classic()+
    theme(axis.title.x = element_blank(),axis.title.y = element_text(size = 28),axis.text = element_text(size = 30),panel.border =  element_rect(fill = NA,linetype="solid"),
          legend.position = "none",strip.text = element_text(size = 36),strip.placement = "outside",strip.background = element_rect(linetype="blank"),
          panel.spacing = unit(0, "lines")))


Petri_rich_over_1000_nan_mod=lmer((value)~collect_date*rain_trt*variable+(1|sampleID_bact)+(1|gh_block),data=data.frame(Mar_leaf.fung_decon_pr_petr_S.rar_overlap_Rich_M_nan))
plot(Petri_rich_over_1000_nan_mod)
hist(resid(Petri_rich_over_1000_nan_mod))
qqPlot(resid(Petri_rich_over_1000_nan_mod))
shapiro.test(resid(Petri_rich_over_1000_nan_mod))
#p-value = 0.6108

#####Table S3 ASV Overlap####
anova(Petri_rich_over_1000_nan_mod)
#nada sig 



#Abundance overlap

Mar_leaf.fung_decon_pr_petr_S.rar_overlap_Abun_M=melt(Mar_leaf.fung_decon_pr_petr_S.rar_overlap[c("rain_trt","collect_date","sampleID_bact",
                                                                                                  "overlap_seed_sum_prop","overlap_rain_sum_prop",
                                                                                                  "gh_block")])


#Reads


Mar_leaf.fung_decon_pr_petr_S.rar_overlap_Abun_sum=Mar_leaf.fung_decon_pr_petr_S.rar_overlap_Abun_M%>%
  group_by(rain_trt,collect_date,variable)%>%
  summarise_at(vars("value"),list(~mean(.),se=~sd(.)/sqrt(n()),~n(),~sd(.)))
Rain_order=c("Nano","Sterile_rain","Live_Rain")

with(Mar_leaf.fung_decon_pr_petr_S.rar_overlap_Abun_sum,interaction(rain_trt,collect_date))
rain_trt_over_round=c("Nano.8/14/2018","Nano.9/18/2018","Sterile_rain.8/14/2018","Sterile_rain.9/18/2018","Live_Rain.8/14/2018","Live_Rain.9/18/2018")
brewer.pal(n = 12, name = "Paired")



#No nanopure samples 
Rain_order_nan=c("Sterile_rain","Live_Rain")
Mar_leaf.fung_decon_pr_petr_S.rar_overlap_Abun_M_nan=subset(Mar_leaf.fung_decon_pr_petr_S.rar_overlap_Abun_M,rain_trt!="Nano")
Mar_leaf.fung_decon_pr_petr_S.rar_overlap_Abun_sum_nan=subset(Mar_leaf.fung_decon_pr_petr_S.rar_overlap_Abun_sum,rain_trt!="Nano")



#1200*600


(abun_over_petri_nan_p=ggplot(Mar_leaf.fung_decon_pr_petr_S.rar_overlap_Abun_M_nan,
                              aes(x=factor(variable,levels = c("overlap_seed_sum_prop","overlap_rain_sum_prop")),y=value))+
    geom_point(aes(shape=factor(rain_trt,levels = Rain_order_nan),fill=factor(rain_trt,levels = Rain_order_nan),
                   group=factor(rain_trt,levels = Rain_order_nan)), 
               color="black",size=4,position = position_dodge(0.65),alpha=.15)+
    geom_errorbar(data=Mar_leaf.fung_decon_pr_petr_S.rar_overlap_Abun_sum_nan, 
                  aes(x=factor(variable,levels = c("overlap_seed_sum_prop","overlap_rain_sum_prop")),group=factor(rain_trt,levels = Rain_order_nan),
                      y=mean,ymin=mean-se,ymax=mean+se),width=.5,
                  color="black",position = position_dodge(0.65))+
    geom_point(data=Mar_leaf.fung_decon_pr_petr_S.rar_overlap_Abun_sum_nan,
               aes(y=mean,shape=factor(rain_trt,levels = Rain_order_nan),
                   group=factor(rain_trt,levels = Rain_order_nan),
                   fill=factor(rain_trt,levels = Rain_order_nan),x=factor(variable,levels = c("overlap_seed_sum_prop","overlap_rain_sum_prop"))), size=10,
               color="black",position = position_dodge(0.65))+facet_nested(.~(factor(collect_date,labels = c("Round1","Round2"))))+
    scale_shape_manual(values = c(25,24),name=NULL)+scale_x_discrete(labels=c("Seed","Rain"))+
    scale_fill_manual(values = c("#6A3D9A","#1F78B4"),name=NULL)+
    ylab("Proportion reads overlap")+theme_classic()+
    theme(axis.title.x = element_blank(),axis.title.y = element_text(size = 28),axis.text = element_text(size = 30),
          panel.border =  element_rect(fill = NA,linetype="solid"),
          legend.position = "none",strip.text = element_text(size = 36),strip.placement = "outside",strip.background = element_rect(linetype="blank"),
          panel.spacing = unit(0, "lines")))


Petri_reads_over_1000_nan_mod=lmer((value)~collect_date*rain_trt*variable+(1|sampleID_bact)+(1|gh_block),data=data.frame(Mar_leaf.fung_decon_pr_petr_S.rar_overlap_Abun_M_nan))
plot(Petri_reads_over_1000_nan_mod)
hist(resid(Petri_reads_over_1000_nan_mod))
qqPlot(resid(Petri_reads_over_1000_nan_mod))
shapiro.test(resid(Petri_reads_over_1000_nan_mod))
#W = 0.97003, p-value = 0.2539

#####Table S3 Read Overlap####
anova(Petri_reads_over_1000_nan_mod)
#rain_trt                       0.112136 0.112136     1 13.180  6.3233 0.025653 * 
#variable                       0.187000 0.187000     1 20.000 10.5448 0.004036 **   


emmeans(Petri_reads_over_1000_nan_mod,pairwise~rain_trt*variable*collect_date,adjust="fdr")

emmeans(Petri_reads_over_1000_nan_mod,pairwise~rain_trt|variable)
emmeans(Petri_reads_over_1000_nan_mod,pairwise~variable|collect_date)
emmeans(Petri_reads_over_1000_nan_mod,pairwise~rain_trt|variable|collect_date)

#####POSTHOC TEST: Figure S2d####
emmeans(Petri_reads_over_1000_nan_mod,pairwise~rain_trt*variable|collect_date,adjust="fdr")



#####Petri Indicator analyses for importance of seed versus rain####

ntaxa(Mar_leaf.fung_decon_pr_petr_S_nan.rar)
#944
nsamples(Mar_leaf.fung_decon_pr_petr_S_nan.rar)
#35

#We Need to extract Rain and Seed OTUs
Mar_leaf.fung_decon_pr_petr_S.rar_Rain_seed=subset_samples(Mar_leaf.fung_decon_pr_petr_S_nan.rar, plant_type=="Seed"|plant_type=="Rain")
nsamples(Mar_leaf.fung_decon_pr_petr_S.rar_Rain_seed)
#11

Mar_leaf.fung_decon_pr_petr_S.rar_Rain_seed=prune_taxa(taxa_sums(Mar_leaf.fung_decon_pr_petr_S.rar_Rain_seed)>0,Mar_leaf.fung_decon_pr_petr_S.rar_Rain_seed)


ntaxa(Mar_leaf.fung_decon_pr_petr_S.rar_Rain_seed)
#713
Mar_leaf.fung_decon_pr_petr_S.rar_Rain_seed_otu=data.frame(as(t(otu_table(Mar_leaf.fung_decon_pr_petr_S.rar_Rain_seed)),"matrix"))
P_rain_seed_order=sample_data(Mar_leaf.fung_decon_pr_petr_S.rar_Rain_seed)$plant_type

P_rain_seed.ind = multipatt(Mar_leaf.fung_decon_pr_petr_S.rar_Rain_seed_otu, 
                            P_rain_seed_order,duleg=TRUE, 
                          control = how(nperm=9999), func="IndVal.g")
summary(P_rain_seed.ind, indvalcomp=TRUE)

head(P_rain_seed.ind$str)
#Component `A' is the probability that the surveyed
#site belongs to the target site group given the fact that the species has
#been found. This conditional probability is called the specificity or positive
#predictive value of the species as indicator of the site group.

head(P_rain_seed.ind$A)
#Component `B' is the probability of finding the species in sites belonging to the site group.
#This second conditional probability is called the fidelity or sensitivity of the
#species as indicator of the target site group.
head(P_rain_seed.ind$B)
head(P_rain_seed.ind$sign)


#Let's create our own indicator values based off of A and B since multipatt square root transforms and groups automatically 

P_rain_seed_Specif=P_rain_seed.ind$A
colnames(P_rain_seed_Specif)=c("Rain_spec","Seed_spec")
P_rain_seed_Fidel=P_rain_seed.ind$B
colnames(P_rain_seed_Fidel)=c("Rain_fid","Seed_fid")
P_rain_seed_ind_comb=merge(P_rain_seed_Specif,P_rain_seed_Fidel, by="row.names")
colnames(P_rain_seed_ind_comb)[colnames(P_rain_seed_ind_comb)=="Row.names"]="OTUs"
P_rain_seed_ind_comb$Rain_IndV=P_rain_seed_ind_comb$Rain_spec*P_rain_seed_ind_comb$Rain_fid
P_rain_seed_ind_comb$Seed_IndV=P_rain_seed_ind_comb$Seed_spec*P_rain_seed_ind_comb$Seed_fid
head(P_rain_seed_ind_comb)

#Let's add in the significance levels

P_rain_seed_ind_comb=merge(P_rain_seed_ind_comb,P_rain_seed.ind$sign,by.x = "OTUs",by.y = "row.names")


write.csv(P_rain_seed_ind_comb, here::here("R_file","IndVal_Petri_rain_seed_comb_taxa_class.csv"),row.names = F)

P_rain_seed_ind_comb=read.csv(here::here("R_file","IndVal_Petri_rain_seed_comb_taxa_class.csv"))

ntaxa(Mar_leaf.fung_decon_pr_petr_S_nan.rar)
#944
nsamples(Mar_leaf.fung_decon_pr_petr_S_nan.rar)
#35
#We need OTUs from plants only
Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS=subset_samples(Mar_leaf.fung_decon_pr_petr_S_nan.rar, plant_type!="Seed"&plant_type!="Rain")
nsamples(Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS)
#24

Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS=prune_taxa(taxa_sums(Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS)>0,Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS)


ntaxa(Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS)
#380

Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_otu=data.frame(otu_table(Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS))
head(Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_otu)
Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_otu$OTUs=row.names(Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_otu)
Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_otu_M=melt(Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_otu)

head(Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_otu_M)
colnames(Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_otu_M)=c("OTUs","sampleID","abund")
nrow(Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_otu_M)
#9120

#let's remove zeros to reduce the computational load
Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_otu_M_NZ=subset(Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_otu_M,abund>0)
nrow(Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_otu_M_NZ)
#1127

Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_otu_M_IndV=merge(Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_otu_M_NZ,
                                                                  P_rain_seed_ind_comb,by="OTUs",all.x = T)


head(Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_otu_M_IndV)

#Now let's calculate the community weighted trait mean 
Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_otu_M_IndV$seed_IndV_abun=Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_otu_M_IndV$Seed_IndV*Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_otu_M_IndV$abund
Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_otu_M_IndV$rain_IndV_abun=Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_otu_M_IndV$Rain_IndV*Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_otu_M_IndV$abund
summary(Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_otu_M_IndV)



Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_otu_M_IndV_sum=Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_otu_M_IndV%>%group_by(sampleID)%>%summarise_at(vars(rain_IndV_abun,seed_IndV_abun),~sum(.,na.rm=TRUE))
Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_IndVal=merge(Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_otu_M_IndV_sum,
                                                              sample_data(Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS),by.x = "sampleID",
                                                              by.y = "row.names")

head(Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_IndVal)



#Make the data long so we can compare rain versus seed importance and average
Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_IndVal$rain_IndV_abun_prop=Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_IndVal$rain_IndV_abun/1000
Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_IndVal$seed_IndV_abun_prop=Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_IndVal$seed_IndV_abun/1000
head(Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_IndVal)

Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_IndVal_M_prop=melt(Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_IndVal[,c("sampleID","seed_IndV_abun_prop",
                                                                                                                "rain_IndV_abun_prop","rain_trt",
                                                                                                                "gh_block",
                                                                                                                "gh_block","collect_date")])
head(Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_IndVal_M_prop)
nrow(Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_IndVal_M_prop)
#48

Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_IndVal_M_prop_sum=Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_IndVal_M_prop%>%
  group_by(collect_date,variable,rain_trt)%>%
  summarise_at("value",list(~mean(.),se=~sd(.)/sqrt(n()),~n(),~sd(.)))
Rain_order=c("Sterile_rain","Live_Rain")

round_order=c("Round1","Round2")

rain_trt_date_order=c("Sterile_rain.8/14/2018","Live_Rain.8/14/2018","Sterile_rain.9/18/2018","Live_Rain.9/18/2018")
brewer.pal(n = 12, name = "Paired")


(IndVal_petri_nan_p=ggplot(Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_IndVal_M_prop)+
    geom_point(aes(x=factor(variable,levels=c("seed_IndV_abun_prop","rain_IndV_abun_prop")),y=value,
                   shape=factor(rain_trt,levels = Rain_order),fill=factor(rain_trt,levels = Rain_order),
                   group=factor(interaction(rain_trt,collect_date),levels = rain_trt_date_order)), 
               color="black",size=4,position = position_dodge(0.65),alpha=.15)+
    geom_errorbar(data=Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_IndVal_M_prop_sum, 
                  aes(x=factor(variable,levels = c("seed_IndV_abun_prop","rain_IndV_abun_prop")),
                      group=factor(interaction(rain_trt,collect_date),levels = rain_trt_date_order),
                      y=mean,ymin=mean-se,ymax=mean+se),width=.5,
                  color="black",position = position_dodge(0.65))+
    geom_point(data=Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_IndVal_M_prop_sum,
               aes(shape=factor(rain_trt,levels = Rain_order),group=factor(interaction(rain_trt,collect_date),levels = rain_trt_date_order),
                   fill=factor(rain_trt,levels = Rain_order),x=factor(variable,levels = c("seed_IndV_abun_prop","rain_IndV_abun_prop")),y=mean), size=10,
               color="black",position = position_dodge(0.65))+facet_wrap(vars(factor(collect_date,labels = c("Round1","Round2"))),ncol = 2)+
    scale_shape_manual(values = c(25,24),name=NULL)+scale_x_discrete(labels=c("Seed","Rain"))+
    scale_fill_manual(values = c("#6A3D9A","#1F78B4"),name=NULL)+
    ylab("Average indicator value")+theme_classic()+
    theme(axis.title.x = element_blank(),axis.title.y = element_text(size = 28),axis.text = element_text(size = 30),
          panel.border =  element_rect(fill = NA,linetype="solid"),
          legend.position = "none",strip.text = element_text(size = 36),strip.placement = "outside",strip.background = element_rect(linetype="blank"),
          panel.spacing = unit(0, "lines")))
#Petri_IndVal_RS_LFE
#1200*600



Petri_IndVal_SR_mod=lmer((value)~rain_trt*collect_date*variable+
                           (1|sampleID)+(1|gh_block),data=Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_IndVal_M_prop)
plot(Petri_IndVal_SR_mod)
hist(resid(Petri_IndVal_SR_mod))
qqPlot(resid(Petri_IndVal_SR_mod))
shapiro.test(resid(Petri_IndVal_SR_mod))
#p-value =  0.8326

####Table S4 Average Indicator Value####
anova(Petri_IndVal_SR_mod)

#rain_trt                       0.15671 0.15671     1    40  8.1587 0.0067659 ** 
#rain_trt:variable              0.33891 0.33891     1    40 17.6438 0.0001447 ***


#####POSTHOC TEST: Figure 3a####
emmeans(Petri_IndVal_SR_mod,pairwise~variable*rain_trt|collect_date, adjust="fdr")
#$contrasts
#collect_date = 8/14/2018:
#  contrast             
#rain_IndV_abun_prop Live_Rain - rain_IndV_abun_prop Sterile_rain     0.24439 0.0800 37.1  3.054  0.0250 

#collect_date = 9/18/2018:
#contrast
#seed_IndV_abun_prop Live_Rain - rain_IndV_abun_prop Sterile_rain     0.21322 0.0809 38.1  2.635  0.0264 
#rain_IndV_abun_prop Live_Rain - rain_IndV_abun_prop Sterile_rain     0.32027 0.0809 38.1  3.959  0.0019 
#seed_IndV_abun_prop Sterile_rain - rain_IndV_abun_prop Sterile_rain  0.21763 0.0800 20.0  2.720  0.0264 



#####Significant Indicator Petri analyses#####


ntaxa(Mar_leaf.fung_decon_pr_petr_S_nan.rar)
#944
nsamples(Mar_leaf.fung_decon_pr_petr_S_nan.rar)
#35
#We need OTUs from plants only
Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS=subset_samples(Mar_leaf.fung_decon_pr_petr_S_nan.rar, plant_type!="Seed"&plant_type!="Rain")
nsamples(Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS)
#24

Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS=prune_taxa(taxa_sums(Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS)>0,Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS)


ntaxa(Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS)
#380

Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_otu=data.frame(otu_table(Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS))
head(Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_otu)
Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_otu$OTUs=row.names(Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_otu)
Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_otu_M=melt(Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_otu)

head(Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_otu_M)
colnames(Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_otu_M)=c("OTUs","sampleID","abund")
nrow(Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_otu_M)
#9120

#let's remove zeros to reduce the computational load
Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_otu_M_NZ=subset(Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_otu_M,abund>0)
nrow(Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_otu_M_NZ)
#1127

#Let's load in the indicator matrix and combine that with the below matrices

P_rain_seed_ind_comb=read.csv(here::here("R_file","IndVal_Petri_rain_seed_comb_taxa_class.csv"))

Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_otu_M_IndV=merge(Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_otu_M_NZ,
                                                             P_rain_seed_ind_comb,by="OTUs",all.x = T)


head(Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_otu_M_IndV)
summary(Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_otu_M_IndV)
nrow(Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_otu_M_IndV)
#1127

Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_otu_M_IndV_sig=subset(Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_otu_M_IndV, p.value<0.05)

Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_otu_M_IndV_sig$ind_grp=
  ifelse(Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_otu_M_IndV_sig$index==1,"Rain","Seed")
summary(Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_otu_M_IndV_sig)
Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_otu_M_IndV_sig_sum=
  Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_otu_M_IndV_sig%>%group_by(sampleID,ind_grp)%>%summarise_at("abund",~sum(.))

Mar_leaf.fung_decon_pr_petr_S.rar_IndV_sig_sum_map=merge(Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_otu_M_IndV_sig_sum,
                                                     sample_data(Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS),
                                                     by.x = "sampleID",by.y = "row.names")
head(Mar_leaf.fung_decon_pr_petr_S.rar_IndV_sig_sum_map)

Mar_leaf.fung_decon_pr_petr_S.rar_IndV_sig_sum_map$prop_abund=Mar_leaf.fung_decon_pr_petr_S.rar_IndV_sig_sum_map$abund/1000

summary(Mar_leaf.fung_decon_pr_petr_S.rar_IndV_sig_sum_map)
nrow(Mar_leaf.fung_decon_pr_petr_S.rar_IndV_sig_sum_map)
#48
#Graph the frequency 

Mar_leaf.fung_decon_pr_petr_S.rar_IndV_sig_sum_map_sum=Mar_leaf.fung_decon_pr_petr_S.rar_IndV_sig_sum_map%>%
  group_by(collect_date,ind_grp,rain_trt)%>%
  summarise_at("prop_abund",list(~mean(.),se=~sd(.)/sqrt(n()),~n(),~sd(.)))
Rain_order=c("Sterile_rain","Live_Rain")

round_order=c("Round1","Round2")
with(Mar_leaf.fung_decon_pr_petr_S.rar_IndV_sig_sum_map_sum,interaction(rain_trt,collect_date))
rain_trt_date_order=c("Sterile_rain.8/14/2018","Live_Rain.8/14/2018","Sterile_rain.9/18/2018","Live_Rain.9/18/2018")
brewer.pal(n = 12, name = "Paired")


(IndVal_sig_petri_nan_p=ggplot(Mar_leaf.fung_decon_pr_petr_S.rar_IndV_sig_sum_map)+
    geom_point(aes(x=factor(ind_grp,levels=c("Seed","Rain")),y=prop_abund,
                   shape=factor(rain_trt,levels = Rain_order),fill=factor(rain_trt,levels = Rain_order),
                   group=factor(interaction(rain_trt,collect_date),levels = rain_trt_date_order)), 
               color="black",size=4,position = position_dodge(0.65),alpha=.15)+
    geom_errorbar(data=Mar_leaf.fung_decon_pr_petr_S.rar_IndV_sig_sum_map_sum, 
                  aes(x=factor(ind_grp,levels = c("Seed","Rain")),
                      group=factor(interaction(rain_trt,collect_date),levels = rain_trt_date_order),
                      y=mean,ymin=mean-se,ymax=mean+se),width=.5,
                  color="black",position = position_dodge(0.65))+
    geom_point(data=Mar_leaf.fung_decon_pr_petr_S.rar_IndV_sig_sum_map_sum,
               aes(shape=factor(rain_trt,levels = Rain_order),group=factor(interaction(rain_trt,collect_date),levels = rain_trt_date_order),
                   fill=factor(rain_trt,levels = Rain_order),x=factor(ind_grp,levels = c("Seed","Rain")),y=mean), size=10,
               color="black",position = position_dodge(0.65))+facet_wrap(vars(factor(collect_date,labels = c("Round1","Round2"))),ncol = 2)+
    scale_shape_manual(values = c(25,24),name=NULL)+scale_x_discrete(labels=c("Seed","Rain"))+
    scale_fill_manual(values = c("#6A3D9A","#1F78B4"),name=NULL)+
    ylab("Proportion of reads indicator")+theme_classic()+
    theme(axis.title.x = element_blank(),axis.title.y = element_text(size = 28),axis.text = element_text(size = 30),panel.border =  element_rect(fill = NA,linetype="solid"),
          legend.position = "none",strip.text = element_text(size = 36),strip.placement = "outside",strip.background = element_rect(linetype="blank"),
          panel.spacing = unit(0, "lines")))
#
#1200*600


Mar_leaf.fung_decon_pr_petr_S.rar_IndV_sig_sum_map_sum
#Round2 increase in rain taxa 

(0.319-0.0287)/0.0287

#Analyses

Petri_IndVal_sig_mod=lmer((prop_abund)^(1/3)~rain_trt*collect_date*ind_grp+
                           (1|sampleID)+(1|gh_block),data=Mar_leaf.fung_decon_pr_petr_S.rar_IndV_sig_sum_map)
plot(Petri_IndVal_sig_mod)
hist(resid(Petri_IndVal_sig_mod))
qqPlot(resid(Petri_IndVal_sig_mod))
shapiro.test(resid(Petri_IndVal_sig_mod))
#W = 0.98776, p-value = 0.8931

####Table S4 Abundance of sig Indicator Taxa####
anova(Petri_IndVal_sig_mod)

#rain_trt                      0.199696 0.199696     1 17.972  5.6887 0.028297 * 
#rain_trt:collect_date         0.118225 0.118225     1 17.972  3.3678 0.083088 . 
#rain_trt:ind_grp              0.302852 0.302852     1 20.000  8.6272 0.008147 **

#####POSTHOC TEST: Figure S3a####
emmeans(Petri_IndVal_sig_mod,pairwise~ind_grp*rain_trt|collect_date, adjust="fdr")



#####Jaccard Petri Pariwise####
Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt=read.csv(here::here("R_file","Pairwise_turnover_Jaccard_fung_v4.2.2020_decon_pruned_ZOTU_full_rar_Marshall.csv"))
head(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt)
nrow(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt)
#6786

#Subset to only petri dishes between rounds within trt

Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_dishes=subset(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt,
                                                                   s1_s2_plant_type=="Petri.Petri")
nrow(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_dishes)
#496
unique(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_dishes$s1_s2_collect_date)

Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_dishes_rounds=subset(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_dishes,
                                                                        s1_s2_collect_date=="9/18/2018.8/14/2018"|
                                                                          s1_s2_collect_date=="8/14/2018.9/18/2018")
nrow(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_dishes_rounds)
#256

Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_dishes_rounds_w_trt=subset(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_dishes_rounds,
                                                                                  s1_s2_rain_trt=="Live_Rain.Live_Rain"|s1_s2_rain_trt=="Nano.Nano"|
                                                                                    s1_s2_rain_trt=="Sterile_rain.Sterile_rain")

unique(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_dishes_rounds_w_trt$s1_s2_rain_trt)

Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_dishes_rounds_w_trt$treatment=rep("Between_Rounds")


#Subset Seed to trt and round

unique(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt$s1_s2_plant_type)

Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_seed=subset(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt, 
                                                              s1_s2_plant_type=="Petri.Seed"|s1_s2_plant_type=="Seed.Petri")
nrow(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_seed)
#128

unique(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_seed$s1_s2_rain_trt)


Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_seed$treatment=rep("Seed")


#Subset Rain to trt and round

unique(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt$s1_s2_plant_type)

Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_rain0=subset(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt, 
                                                              s1_s2_plant_type=="Petri.Rain"|s1_s2_plant_type=="Rain.Petri")


nrow(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_rain0)
#1792

unique(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_rain0$s1_s2_collect_date)

#Petri 
#start_date collect_dates 
#8/21/2018	9/18/2018
#7/21/2018	8/14/2018

Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_rain=subset(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_rain0,s1_s2_collect_date=="8/21/2018.9/18/2018"|
                                                                  s1_s2_collect_date=="9/18/2018.8/21/2018"|s1_s2_collect_date=="7/21/2018.8/14/2018"|
                                                                  s1_s2_collect_date=="8/14/2018.7/21/2018")


unique(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_rain$s1_s2_collect_date)
nrow(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_rain)
#112

Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_rain$treatment=rep("Rain")


#combine the matrix 
Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pairwise=rbind(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_dishes_rounds_w_trt,
                                                                      Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_seed,
                                                               Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_rain)


Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pairwise$leaf_sample=ifelse(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pairwise$
                                                                                     s1_plant_type=="Petri",as.character(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pairwise$
                                                                                                                           sample1),as.character(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pairwise$
                                                                                                                                                   sample2))
unique(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pairwise$leaf_sample)
Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pairwise$unique_rain_trt=ifelse(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pairwise$
         s1_plant_type=="Petri",as.character(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pairwise$
                                               s1_rain_trt),as.character(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pairwise$
                                                                           s2_rain_trt))

unique(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pairwise$unique_rain_trt)

Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pairwise$round=ifelse(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pairwise$
                                                                               s1_collect_date=="9/18/2018"|Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pairwise$
                                                                               s2_collect_date=="9/18/2018","Round2","Round1")



head(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pairwise)
nrow(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pairwise)
#328

#Blocking term

Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pairwise$gh_block_comb=
  ifelse(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pairwise$s1_plant_type=="Petri",
         as.character(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pairwise$s1_gh_block),
         as.character(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pairwise$s2_gh_block))

head(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pairwise)
nrow(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pairwise)
#328
unique(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pairwise$gh_block_comb)


#Drop the between the rounds analyses
Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub=subset(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pairwise,treatment!="Between_Rounds")
nrow(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub)
#240


#Jaccard

Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub_sum=Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub%>%
  group_by(leaf_sample,round,unique_rain_trt,treatment,gh_block_comb)%>%
  summarise_at("jaccard",list(~mean(.),se=~sd(.)/sqrt(n()),~n()))

Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub_sum1=data.frame(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub_sum
                                                                                     [,c("round","unique_rain_trt","treatment","mean")])

colnames(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub_sum1)=c("round","unique_rain_trt","treatment","jaccard")
Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub_sum_sum=Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub_sum1%>%
  group_by(round,unique_rain_trt,treatment)%>%
  summarise_all(list(~mean(.),se=~sd(.)/sqrt(n()),~n(),~sd(.)))

#Nanopure samples need to be removed
Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub_sum_nan=subset(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub_sum, unique_rain_trt!="Nano")
Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub_sum_sum_nan=subset(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub_sum_sum,
                                                                                   unique_rain_trt!="Nano")
(jacc1000_petri_nan_p=ggplot(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub_sum_nan,aes(x=factor(treatment,levels = c("Seed","Rain")),y=mean))+
    geom_point(aes(shape=factor(unique_rain_trt,levels = Rain_order),fill=factor(unique_rain_trt,levels = Rain_order),
                   group=factor(interaction(unique_rain_trt,round),levels = rain_trt_round_nan)), 
               color="black",size=4,position = position_dodge(0.65),alpha=.15)+
    geom_errorbar(data=Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub_sum_sum_nan, 
                  aes(x=factor(treatment,levels = c("Seed","Rain")),group=factor(interaction(unique_rain_trt,round),levels = rain_trt_round_nan),
                      y=mean,ymin=mean-se,ymax=mean+se),width=.5,
                  color="black",position = position_dodge(0.65))+
    geom_point(data=Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub_sum_sum_nan,
               aes(shape=factor(unique_rain_trt,levels = Rain_order),group=factor(interaction(unique_rain_trt,round),levels = rain_trt_round_nan),
                   fill=factor(unique_rain_trt,levels = Rain_order),x=factor(treatment,levels = c("Seed","Rain"))), size=10,
               color="black",position = position_dodge(0.65))+facet_nested(.~round)+
    scale_shape_manual(values = c(25,24),name=NULL)+scale_x_discrete(labels=c("Seed","Rain"))+
    scale_fill_manual(values = c("#6A3D9A","#1F78B4"),name=NULL)+
    ylab("Jaccard distance")+theme_classic()+
    theme(axis.title.x = element_blank(),axis.title.y = element_text(size = 28),axis.text = element_text(size = 30),panel.border =  element_rect(fill = NA,linetype="solid"),
          legend.position = "none",strip.text = element_text(size = 36),strip.placement = "outside",strip.background = element_rect(linetype="blank"),
          panel.spacing = unit(0, "lines")))
#Dist_Jacc_petri_dot_error_nan
#1200*600



Petri_pair_dist_J1000_nan_mod=lmer((mean)~round*unique_rain_trt*treatment+(1|leaf_sample)+(1|gh_block_comb),data=Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub_sum_nan)
plot(Petri_pair_dist_J1000_nan_mod)
hist(resid(Petri_pair_dist_J1000_nan_mod))
qqPlot(resid(Petri_pair_dist_J1000_nan_mod))
shapiro.test(resid(Petri_pair_dist_J1000_nan_mod))
#p-value = 0.1747

#####Table S3 Jaccard####
anova(Petri_pair_dist_J1000_nan_mod)
#treatment                       0.134385 0.134385     1 20.000 209.8911 4.56e-12 ***


emmeans(Petri_pair_dist_J1000_nan_mod, pairwise~treatment|round)
emmeans(Petri_pair_dist_J1000_nan_mod, pairwise~treatment|round|unique_rain_trt)
emmeans(Petri_pair_dist_J1000_nan_mod, pairwise~treatment*round)
emmeans(Petri_pair_dist_J1000_nan_mod, pairwise~treatment*round*unique_rain_trt,adjust="fdr")

#####POSTHOC TEST: Figure 2b####
emmeans(Petri_pair_dist_J1000_nan_mod, pairwise~treatment*unique_rain_trt|round,adjust="fdr")


#Nestedness
head(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub)
Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub_nest_sum=Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub%>%
  group_by(leaf_sample,round,unique_rain_trt,treatment,gh_block_comb)%>%
  summarise_at("nestedness",list(~mean(.),se=~sd(.)/sqrt(n()),~n()))


Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub_nest_sum1=data.frame(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub_nest_sum
                                                                                [,c("round","unique_rain_trt","treatment","mean")])

colnames(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub_nest_sum1)=c("round","unique_rain_trt","treatment","nestedness")
Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub_nest_sum_sum=Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub_nest_sum1%>%
  group_by(round,unique_rain_trt,treatment)%>%
  summarise_all(list(~mean(.),se=~sd(.)/sqrt(n()),~n(),~sd(.)))
Rain_order=c("Nano","Sterile_rain","Live_Rain")

round_order=c("Round1","Round2")


brewer.pal(n = 12, name = "Paired")


#Remove Nanopure samples


Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub_nest_sum_nan=subset(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub_nest_sum, 
                                                                                    unique_rain_trt!="Nano")
Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub_nest_sum_sum_nan=subset(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub_nest_sum_sum,
                                                                                   unique_rain_trt!="Nano")
(nest1000_petri_nan_p=ggplot(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub_nest_sum_nan,aes(x=factor(treatment,levels = c("Seed","Rain")),y=mean))+
    geom_point(aes(shape=factor(unique_rain_trt,levels = Rain_order),fill=factor(unique_rain_trt,levels = Rain_order),
                   group=factor(interaction(unique_rain_trt,round),levels = rain_trt_round_nan)), 
               color="black",size=4,position = position_dodge(0.65),alpha=.15)+
    geom_errorbar(data=Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub_nest_sum_sum_nan, 
                  aes(x=factor(treatment,levels = c("Seed","Rain")),group=factor(interaction(unique_rain_trt,round),levels = rain_trt_round_nan),
                      y=mean,ymin=mean-se,ymax=mean+se),width=.5,
                  color="black",position = position_dodge(0.65))+
    geom_point(data=Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub_nest_sum_sum_nan,
               aes(shape=factor(unique_rain_trt,levels = Rain_order),group=factor(interaction(unique_rain_trt,round),levels = rain_trt_round_nan),
                   fill=factor(unique_rain_trt,levels = Rain_order),x=factor(treatment,levels = c("Seed","Rain"))), size=10,
               color="black",position = position_dodge(0.65))+facet_nested(.~round)+
    scale_shape_manual(values = c(25,24),name=NULL)+scale_x_discrete(labels=c("Seed","Rain"))+
    scale_fill_manual(values = c("#6A3D9A","#1F78B4"),name=NULL)+
    ylab("Nestedness distance")+theme_classic()+
    theme(axis.title.x = element_blank(),axis.title.y = element_text(size = 28),axis.text = element_text(size = 30),panel.border =  element_rect(fill = NA,linetype="solid"),
          legend.position = "none",strip.text = element_text(size = 36),strip.placement = "outside",strip.background = element_rect(linetype="blank"),
          panel.spacing = unit(0, "lines")))

#1200*600


Petri_pair_dist_nest1000_nan_mod=lmer((mean)~round*unique_rain_trt*treatment+(1|leaf_sample)+(1|gh_block_comb),data=Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub_nest_sum_nan)
plot(Petri_pair_dist_nest1000_nan_mod)
hist(resid(Petri_pair_dist_nest1000_nan_mod))
qqPlot(resid(Petri_pair_dist_nest1000_nan_mod))
shapiro.test(resid(Petri_pair_dist_nest1000_nan_mod))
#W = 0.97145, p-value = 0.2881


#####Table S3 Nestedness####
anova(Petri_pair_dist_nest1000_nan_mod)
#treatment                       0.231011 0.231011     1    20 92.8984 5.856e-09 *** 

#changed
emmeans(Petri_pair_dist_nest1000_nan_mod, pairwise~treatment*round*unique_rain_trt,adjust="fdr")
emmeans(Petri_pair_dist_nest1000_nan_mod, pairwise~treatment|round)

#####POSTHOC TEST: Figure S2a####
emmeans(Petri_pair_dist_nest1000_nan_mod, pairwise~treatment*unique_rain_trt|round,adjust="fdr")





#Turnover
head(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub)
Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub_turn_sum=Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub%>%
  group_by(leaf_sample,round,unique_rain_trt,treatment,gh_block_comb)%>%
  summarise_at("turnover",list(~mean(.),se=~sd(.)/sqrt(n()),~n()))



Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub_turn_sum1=data.frame(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub_turn_sum
                                                                                     [,c("round","unique_rain_trt","treatment","mean")])

colnames(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub_turn_sum1)=c("round","unique_rain_trt","treatment","turnover")
Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub_turn_sum_sum=Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub_turn_sum1%>%
  group_by(round,unique_rain_trt,treatment)%>%
  summarise_all(list(~mean(.),se=~sd(.)/sqrt(n()),~n(),~sd(.)))
Rain_order=c("Nano","Sterile_rain","Live_Rain")
round_order=c("Round1","Round2")




#Remove Nanopure samples


Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub_turn_sum_nan=subset(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub_turn_sum, 
                                                                                    unique_rain_trt!="Nano")
Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub_turn_sum_sum_nan=subset(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub_turn_sum_sum,
                                                                                        unique_rain_trt!="Nano")
(turn1000_petri_nan_p=ggplot(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub_turn_sum_nan,aes(x=factor(treatment,levels = c("Seed","Rain")),y=mean))+
    geom_point(aes(shape=factor(unique_rain_trt,levels = Rain_order),fill=factor(unique_rain_trt,levels = Rain_order),
                   group=factor(interaction(unique_rain_trt,round),levels = rain_trt_round_nan)), 
               color="black",size=4,position = position_dodge(0.65),alpha=.15)+
    geom_errorbar(data=Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub_turn_sum_sum_nan, 
                  aes(x=factor(treatment,levels = c("Seed","Rain")),group=factor(interaction(unique_rain_trt,round),levels = rain_trt_round_nan),
                      y=mean,ymin=mean-se,ymax=mean+se),width=.5,
                  color="black",position = position_dodge(0.65))+
    geom_point(data=Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub_turn_sum_sum_nan,
               aes(shape=factor(unique_rain_trt,levels = Rain_order),group=factor(interaction(unique_rain_trt,round),levels = rain_trt_round_nan),
                   fill=factor(unique_rain_trt,levels = Rain_order),x=factor(treatment,levels = c("Seed","Rain"))), size=10,
               color="black",position = position_dodge(0.65))+facet_nested(.~round)+
    scale_shape_manual(values = c(25,24),name=NULL)+scale_x_discrete(labels=c("Seed","Rain"))+
    scale_fill_manual(values = c("#6A3D9A","#1F78B4"),name=NULL)+
    ylab("Turnover distance")+theme_classic()+
    theme(axis.title.x = element_blank(),axis.title.y = element_text(size = 28),axis.text = element_text(size = 30),panel.border =  element_rect(fill = NA,linetype="solid"),
          legend.position = "none",strip.text = element_text(size = 36),strip.placement = "outside",strip.background = element_rect(linetype="blank"),
          panel.spacing = unit(0, "lines")))
#Dist_Jacc_petri_dot_error_nan
#1200*600



Petri_pair_dist_turn1000_nan_mod=lmer((mean)^3~round*unique_rain_trt*treatment+(1|leaf_sample)+(1|gh_block_comb),data=Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub_turn_sum_nan)
plot(Petri_pair_dist_turn1000_nan_mod)
hist(resid(Petri_pair_dist_turn1000_nan_mod))
qqPlot(resid(Petri_pair_dist_turn1000_nan_mod))
shapiro.test(resid(Petri_pair_dist_turn1000_nan_mod))
#p-value =  0.96

#####Table S3 Turnover####
anova(Petri_pair_dist_turn1000_nan_mod)
#treatment                       0.0226145 0.0226145     1 20.000  5.3377 0.03166 *
#round:treatment                 0.0193902 0.0193902     1 20.000  4.5767 0.04493 *


emmeans(Petri_pair_dist_turn1000_nan_mod,pairwise~treatment|round)
emmeans(Petri_pair_dist_turn1000_nan_mod,pairwise~treatment|unique_rain_trt)
emmeans(Petri_pair_dist_turn1000_nan_mod, pairwise~treatment*round*unique_rain_trt,adjust="fdr")

#####POSTHOC TEST: Figure S2c####
emmeans(Petri_pair_dist_turn1000_nan_mod, pairwise~treatment*unique_rain_trt|round,adjust="fdr")

#####PLOT: Figure S2####
#####Plot Turnover Nestedness Overlap####

plot_grid(nest1000_petri_nan_p,rich_over_petri_nan_p,turn1000_petri_nan_p,abun_over_petri_nan_p,nrow = 2, align="v")
#2000*1200
#Nest_Turn_OverLap_Petri_dot_error_comb_nan_raw


#####Petri Ratio of Nestedness to Turnover####
Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub$ratio_nest_turn=Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub$nestedness/Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub$turnover

head(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub)
Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub_rat_NT_sum=Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub%>%
  group_by(leaf_sample,round,unique_rain_trt,treatment,gh_block_comb)%>%
  summarise_at("ratio_nest_turn",list(~mean(.),se=~sd(.)/sqrt(n()),~n()))



Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub_rat_NT_sum1=data.frame(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub_rat_NT_sum
                                                                                     [,c("round","unique_rain_trt","treatment","mean")])

colnames(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub_rat_NT_sum1)=c("round","unique_rain_trt","treatment","ratio_nest_turn")
Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub_rat_NT_sum_sum=Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub_rat_NT_sum1%>%
  group_by(round,unique_rain_trt,treatment)%>%
  summarise_all(list(~mean(.),se=~sd(.)/sqrt(n()),~n(),~sd(.)))
Rain_order=c("Nano","Sterile_rain","Live_Rain")
with(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub_sum,interaction(round,treatment))
round_order=c("Round1","Round2")
rain_trt_round_nan=c("Sterile_rain.Round1","Live_Rain.Round1","Sterile_rain.Round2","Live_Rain.Round2")
with(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub_sum,interaction(unique_rain_trt,round))

brewer.pal(n = 12, name = "Paired")


#Remove Nanopure samples


Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub_rat_NT_sum_nan=subset(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub_rat_NT_sum, 
                                                                                    unique_rain_trt!="Nano")
Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub_rat_NT_sum_sum_nan=subset(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub_rat_NT_sum_sum,
                                                                                        unique_rain_trt!="Nano")
(rat_NT_1000_petri_nan_p=ggplot(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub_rat_NT_sum_nan,aes(x=factor(treatment,levels = c("Seed","Rain")),y=mean))+
    geom_point(aes(shape=factor(unique_rain_trt,levels = Rain_order),fill=factor(unique_rain_trt,levels = Rain_order),
                   group=factor(interaction(unique_rain_trt,round),levels = rain_trt_round_nan)), 
               color="black",size=4,position = position_dodge(0.65),alpha=.15)+
    geom_errorbar(data=Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub_rat_NT_sum_sum_nan, 
                  aes(x=factor(treatment,levels = c("Seed","Rain")),group=factor(interaction(unique_rain_trt,round),levels = rain_trt_round_nan),
                      y=mean,ymin=mean-se,ymax=mean+se),width=.5,
                  color="black",position = position_dodge(0.65))+
    geom_point(data=Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub_rat_NT_sum_sum_nan,
               aes(shape=factor(unique_rain_trt,levels = Rain_order),group=factor(interaction(unique_rain_trt,round),levels = rain_trt_round_nan),
                   fill=factor(unique_rain_trt,levels = Rain_order),x=factor(treatment,levels = c("Seed","Rain"))), size=10,
               color="black",position = position_dodge(0.65))+facet_nested(.~round)+
    scale_shape_manual(values = c(25,24),name=NULL)+scale_x_discrete(labels=c("Seed","Rain"))+
    scale_fill_manual(values = c("#6A3D9A","#1F78B4"),name=NULL)+
    ylab("Ratio of\nNestedness to Turnover")+theme_classic()+
    theme(axis.title.x = element_blank(),axis.title.y = element_text(size = 28),axis.text = element_text(size = 30),panel.border =  element_rect(fill = NA,linetype="solid"),
          legend.position = "none",strip.text = element_text(size = 36),strip.placement = "outside",strip.background = element_rect(linetype="blank"),
          panel.spacing = unit(0, "lines")))
#dot_Ratio_Nest_Turn_Petri_dot_error_raw
#1200*600



Petri_pair_dist_rat_NT_1000_nan_mod=lmer(sqrt(mean)~round*unique_rain_trt*treatment+(1|leaf_sample)+(1|gh_block_comb),data=Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub_rat_NT_sum_nan)
plot(Petri_pair_dist_rat_NT_1000_nan_mod)
hist(resid(Petri_pair_dist_rat_NT_1000_nan_mod))
qqPlot(resid(Petri_pair_dist_rat_NT_1000_nan_mod))
shapiro.test(resid(Petri_pair_dist_rat_NT_1000_nan_mod))
#W = 0.97267, p-value = 0.3206

#####Table S3 Ratio Nestedness:Turnover####
anova(Petri_pair_dist_rat_NT_1000_nan_mod)
#treatment                       0.63979 0.63979     1 20.000 79.2153 2.16e-08 ***


emmeans(Petri_pair_dist_rat_NT_1000_nan_mod, pairwise~treatment*round*unique_rain_trt,adjust="fdr")
emmeans(Petri_pair_dist_rat_NT_1000_nan_mod, pairwise~treatment|round)


#####POSTHOC TEST: Figure 2c####
emmeans(Petri_pair_dist_rat_NT_1000_nan_mod, pairwise~treatment*unique_rain_trt|round,adjust="fdr")




#####PLOT: Figure 2####
#####Plot Petri Bray Jaccard and Ratio#####



(bray1000_petri_nan_p=ggplot(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_pair_sub_sum_nan,aes(x=factor(treatment,levels = c("Seed","Rain")),y=mean))+
    geom_point(aes(shape=factor(unique_rain_trt,levels = Rain_order),fill=factor(unique_rain_trt,levels = Rain_order),
                   group=factor(interaction(unique_rain_trt,round),levels = rain_trt_round_nan)), 
               color="black",size=4,position = position_dodge(0.65),alpha=.15)+
    geom_errorbar(data=Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_pair_sub_sum_sum_nan, 
                  aes(x=factor(treatment,levels = c("Seed","Rain")),group=factor(interaction(unique_rain_trt,round),levels = rain_trt_round_nan),
                      y=mean,ymin=mean-se,ymax=mean+se),width=.5,
                  color="black",position = position_dodge(0.65))+
    geom_point(data=Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_petri_pair_sub_sum_sum_nan,
               aes(shape=factor(unique_rain_trt,levels = Rain_order),group=factor(interaction(unique_rain_trt,round),levels = rain_trt_round_nan),
                   fill=factor(unique_rain_trt,levels = Rain_order),x=factor(treatment,levels = c("Seed","Rain"))), size=10,
               color="black",position = position_dodge(0.65))+facet_nested(.~round)+
    scale_shape_manual(values = c(25,24),name=NULL)+scale_x_discrete(labels=c("Seed","Rain"))+
    scale_fill_manual(values = c("#6A3D9A","#1F78B4"),name=NULL)+
    ylab("Bray-Curtis distance")+theme_classic()+
    theme(axis.title.x = element_blank(),axis.title.y = element_text(size = 28),axis.text.y = element_text(size = 30),axis.text.x = element_blank(),
          panel.border =  element_rect(fill = NA,linetype="solid"),
          legend.position = "none",strip.text = element_text(size = 36),strip.placement = "outside",strip.background = element_rect(linetype="blank"),
          panel.spacing = unit(0, "lines")))


(jacc1000_petri_nan_p=ggplot(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub_sum_nan,aes(x=factor(treatment,levels = c("Seed","Rain")),y=mean))+
    geom_point(aes(shape=factor(unique_rain_trt,levels = Rain_order),fill=factor(unique_rain_trt,levels = Rain_order),
                   group=factor(interaction(unique_rain_trt,round),levels = rain_trt_round_nan)), 
               color="black",size=4,position = position_dodge(0.65),alpha=.15)+
    geom_errorbar(data=Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub_sum_sum_nan, 
                  aes(x=factor(treatment,levels = c("Seed","Rain")),group=factor(interaction(unique_rain_trt,round),levels = rain_trt_round_nan),
                      y=mean,ymin=mean-se,ymax=mean+se),width=.5,
                  color="black",position = position_dodge(0.65))+
    geom_point(data=Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub_sum_sum_nan,
               aes(shape=factor(unique_rain_trt,levels = Rain_order),group=factor(interaction(unique_rain_trt,round),levels = rain_trt_round_nan),
                   fill=factor(unique_rain_trt,levels = Rain_order),x=factor(treatment,levels = c("Seed","Rain"))), size=10,
               color="black",position = position_dodge(0.65))+facet_nested(.~round)+
    scale_shape_manual(values = c(25,24),name=NULL)+scale_x_discrete(labels=c("Seed","Rain"))+
    scale_fill_manual(values = c("#6A3D9A","#1F78B4"),name=NULL)+
    ylab("Jaccard distance")+theme_classic()+
    theme(axis.title.x = element_blank(),axis.title.y = element_text(size = 28),axis.text.y = element_text(size = 30),axis.text.x = element_blank(),
            panel.border =  element_rect(fill = NA,linetype="solid"),
          legend.position = "none",strip.text = element_blank(),strip.placement = "none",
          panel.spacing = unit(0, "lines")))
##900x800


(rat_NT_1000_petri_nan_p=ggplot(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub_rat_NT_sum_nan,aes(x=factor(treatment,levels = c("Seed","Rain")),y=mean))+
    geom_point(aes(shape=factor(unique_rain_trt,levels = Rain_order),fill=factor(unique_rain_trt,levels = Rain_order),
                   group=factor(interaction(unique_rain_trt,round),levels = rain_trt_round_nan)), 
               color="black",size=4,position = position_dodge(0.65),alpha=.15)+
    geom_errorbar(data=Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub_rat_NT_sum_sum_nan, 
                  aes(x=factor(treatment,levels = c("Seed","Rain")),group=factor(interaction(unique_rain_trt,round),levels = rain_trt_round_nan),
                      y=mean,ymin=mean-se,ymax=mean+se),width=.5,
                  color="black",position = position_dodge(0.65))+
    geom_point(data=Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_petri_pair_sub_rat_NT_sum_sum_nan,
               aes(shape=factor(unique_rain_trt,levels = Rain_order),group=factor(interaction(unique_rain_trt,round),levels = rain_trt_round_nan),
                   fill=factor(unique_rain_trt,levels = Rain_order),x=factor(treatment,levels = c("Seed","Rain"))), size=10,
               color="black",position = position_dodge(0.65))+facet_nested(.~round)+
    scale_shape_manual(values = c(25,24),name=NULL)+scale_x_discrete(labels=c("Seed","Rain"))+
    scale_fill_manual(values = c("#6A3D9A","#1F78B4"),name=NULL)+
    scale_y_continuous(name = "Ratio\nNestedness:Turnover",breaks = c(seq(0,1,by=0.3)))+theme_classic()+
    theme(axis.title.x = element_blank(),axis.title.y = element_text(size = 28),axis.text = element_text(size = 30),panel.border =  element_rect(fill = NA,linetype="solid"),
          legend.position = "none",strip.text = element_blank(),strip.placement = "none",
          panel.spacing = unit(0, "lines")))


plot_grid(bray1000_petri_nan_p,jacc1000_petri_nan_p,rat_NT_1000_petri_nan_p,nrow = 3, align="v",rel_heights = c(0.95,0.8,0.9))
#1100*1300

#Petri_Bray_Jac_Ratio_dot_raw



#####Field experiment####


load(here::here("R_file","Mar_leaf_rain_ZOTU_full_v4.2.2020.fung_decon_rar_pruned_phyloseq_obj.RData"))
head(sample_data(Mar_leaf_rain.fung_decon_rar))

Mar_leaf_rain.fung_decon_rar_Mar=subset_samples(Mar_leaf_rain.fung_decon_rar, sub_proj=="Marshall")
nsamples(Mar_leaf_rain.fung_decon_rar_Mar)
#117

Mar_leaf_rain.fung_decon_rar_Mar=prune_taxa(taxa_sums(Mar_leaf_rain.fung_decon_rar_Mar) > 0, Mar_leaf_rain.fung_decon_rar_Mar)
ntaxa(Mar_leaf_rain.fung_decon_rar_Mar)
#2484
sum(taxa_sums(Mar_leaf_rain.fung_decon_rar_Mar))
#117000
min(sample_sums(Mar_leaf_rain.fung_decon_rar_Mar))
#1000
max(sample_sums(Mar_leaf_rain.fung_decon_rar_Mar))
#1000
#Let's look at Rain and the field experiment
head(sample_data(Mar_leaf_rain.fung_decon_rar))
Mar_leaf_rain.fung_decon_rar_Mar_field=subset_samples(Mar_leaf_rain.fung_decon_rar_Mar, plant_type!="Petri")
nsamples(Mar_leaf_rain.fung_decon_rar_Mar_field)
#85

unique(sample_data(Mar_leaf_rain.fung_decon_rar_Mar_field)$sub_proj)

Mar_leaf_rain.fung_decon_rar_Mar_field=prune_taxa(taxa_sums(Mar_leaf_rain.fung_decon_rar_Mar_field) > 0, Mar_leaf_rain.fung_decon_rar_Mar_field)
ntaxa(Mar_leaf_rain.fung_decon_rar_Mar_field)
#2243
sum(taxa_sums(Mar_leaf_rain.fung_decon_rar_Mar_field))
#85000
min(sample_sums(Mar_leaf_rain.fung_decon_rar_Mar_field))
#1000
max(sample_sums(Mar_leaf_rain.fung_decon_rar_Mar_field))
#1000


#final number of samples in each group of seedlings

sample_data(subset_samples(Mar_leaf_rain.fung_decon_rar_Mar_field,plant_type=="Seedling")) %>% group_by(plant_type, collect_date, field_plot) %>% summarise_at("sample_type",~n())


#####Field experiment Bray analyses####


Mar_leaf_rain.fung_decon_rar_Mar_field.ord=ordinate(Mar_leaf_rain.fung_decon_rar_Mar_field,method = "NMDS",distance = "bray")
#*** No convergence -- monoMDS stopping criteria:
#19: stress ratio > sratmax
#1: scale factor of the gradient < sfgrmin
#0.1492751  

#Create files for PERMANOVA analyses of commmunity composition in PRIMER v6
Mar_leaf_rain.fung_decon_rar_Mar_field_dis=phyloseq::distance(Mar_leaf_rain.fung_decon_rar_Mar_field,method = "bray")
Mar_leaf_rain.fung_decon_rar_Mar_field_map=sample_data(Mar_leaf_rain.fung_decon_rar_Mar_field)

write.csv(as.matrix(Mar_leaf_rain.fung_decon_rar_Mar_field_dis),file = here::here("R_file","Bray_Mar_leaf_rain.fung_decon_rar_Mar_field_dis.csv"))
write.csv(Mar_leaf_rain.fung_decon_rar_Mar_field_map,file = here::here("R_file","Mar_leaf_rain.fung_decon_rar_Mar_field_map.csv"))


#####Field Bray Beta-dipsersion####

Mar_leaf_rain.fung_decon_rar_Mar_field_dis_betamod=betadisper(Mar_leaf_rain.fung_decon_rar_Mar_field_dis,
                                                              with(Mar_leaf_rain.fung_decon_rar_Mar_field_map,interaction(collect_date,plant_type)))

Mar_leaf_rain.fung_decon_rar_Mar_field_dis_betadisp=data.frame(Mar_leaf_rain.fung_decon_rar_Mar_field_dis_betamod$distances)

colnames(Mar_leaf_rain.fung_decon_rar_Mar_field_dis_betadisp)="betadisp"

Mar_leaf_rain.fung_decon_rar_Mar_field_dis_betadisp=merge(Mar_leaf_rain.fung_decon_rar_Mar_field_dis_betadisp,Mar_leaf_rain.fung_decon_rar_Mar_field_map,
                                                          by="row.names")


#####Graphing Bray Beta-dipsersion####

Mar_leaf_rain.fung_decon_rar_Mar_field_dis_betadisp_sum=Mar_leaf_rain.fung_decon_rar_Mar_field_dis_betadisp%>%group_by(collect_date,plant_type)%>%
  summarise_at(vars("betadisp"),list(~mean(.),~n(),se=~sd(.)/sqrt(n())))

with(Mar_leaf_rain.fung_decon_rar_Mar_field_dis_betadisp_sum,interaction(collect_date,plant_type))

plant_date_order=c("7/10/2018.Seed","7/16/2018.Adult","9/5/2018.Adult","7/16/2018.Seedling","9/5/2018.Seedling")

#Bray Beta disp
#plants
(B_field_betadisp=ggplot(subset(Mar_leaf_rain.fung_decon_rar_Mar_field_dis_betadisp,plant_type!="Rain"))+
    geom_point(aes(x=factor(interaction(collect_date,plant_type),levels = plant_date_order,
                            labels = c("Seed", "Start","End","Start","End")),y=betadisp, shape=plant_type,
                   fill=factor(interaction(collect_date,plant_type),
                               levels = plant_date_order),
                   color=factor(interaction(collect_date,plant_type),
                                levels = plant_date_order)),size=2, alpha = 0.5)+
    geom_errorbar(data = subset(Mar_leaf_rain.fung_decon_rar_Mar_field_dis_betadisp_sum,plant_type!="Rain"),
                  aes(x=factor(interaction(collect_date,plant_type),levels = plant_date_order,
                               labels = c("Seed", "Start","End","Start","End")),
                      y=mean, ymax=mean+se,
                      ymin=mean-se), width=0.3)+
    geom_point(data = subset(Mar_leaf_rain.fung_decon_rar_Mar_field_dis_betadisp_sum,plant_type!="Rain"),
               aes(x=factor(interaction(collect_date,plant_type),
                            levels = plant_date_order,
                            labels = c("Seed", "Start","End","Start","End")),
                   y=mean, shape=plant_type,
                   fill=factor(interaction(collect_date,plant_type),
                               levels = plant_date_order),
                   color=factor(interaction(collect_date,plant_type),
                                levels = plant_date_order)),size=10)+scale_x_discrete(labels=)+
    scale_shape_manual(values = c(22,13,24),name=NULL)+ylab("Beta-dispersion\n(Bray-Curtis)")+
    facet_grid(cols = vars(factor(plant_type,levels = c("Seed","Adult","Seedling"))),scales = "free_x",space = "free_x")+
    scale_color_manual(values = c("#E31A1C","black", "black","black","black"),name=NULL)+
    scale_fill_manual(values = c("#E31A1C","#A6CEE3","#1F78B4","#B2DF8A","#33A02C"),name=NULL)+
    theme_cowplot()+theme(axis.title.y = element_text(size = 34),axis.title.x = element_blank(),axis.text =  element_text(size = 30),
                          legend.position = "none",strip.background = element_rect(fill = NA),strip.text = element_text(size = 34)))

#1000*800




#Rain
(B_rain_betadisp=ggplot(subset(Mar_leaf_rain.fung_decon_rar_Mar_field_dis_betadisp,plant_type=="Rain"))+
    geom_point(aes(x=as.Date(collect_date,format="%m/%d/%Y"),y=betadisp),size=2, alpha = 0.5, shape=21, fill="black")+
    geom_errorbar(data = subset(Mar_leaf_rain.fung_decon_rar_Mar_field_dis_betadisp_sum,plant_type=="Rain"),aes(x=as.Date(collect_date,format="%m/%d/%Y"),
                                                                                                                y=mean, ymax=mean+se,
                                                                                                                ymin=mean-se), width=1.5)+
    geom_point(data = subset(Mar_leaf_rain.fung_decon_rar_Mar_field_dis_betadisp_sum,plant_type=="Rain"),aes(x=as.Date(collect_date,format="%m/%d/%Y"),
                                                                                                             y=mean),size=10,shape=21, fill="black")+
    theme_cowplot()+scale_x_date(date_labels = "%b %d",date_breaks = "11 days")+ylab("Beta-dispersion\n(Bray-Curtis)")+
    facet_grid(cols = vars(factor(plant_type,labels = "")))+
    theme(axis.title.y = element_text(size = 34),axis.title.x = element_blank(),axis.text =  element_text(size = 30),
          legend.position = "none",strip.background = element_rect(fill = NA),strip.text = element_text(size = 34)))

#1300*800
#



#####Field Betadisp Stats####




Mar_leaf_rain.fung_decon_rar_Mar_field_dis_betadisp_plant=subset(Mar_leaf_rain.fung_decon_rar_Mar_field_dis_betadisp,plant_type!="Rain")
Mar_leaf_rain.fung_decon_rar_Mar_field_dis_betadisp_plant$plant_date=with(Mar_leaf_rain.fung_decon_rar_Mar_field_dis_betadisp_plant,
                                                                          interaction(plant_type,collect_date))

#Beta disp all plants
Betdisp_field_r1000_mod=lmer(sqrt(betadisp)~plant_date+(1|gh_block),data = Mar_leaf_rain.fung_decon_rar_Mar_field_dis_betadisp_plant)
hist(resid(Betdisp_field_r1000_mod))
qqPlot(resid(Betdisp_field_r1000_mod))
shapiro.test(resid(Betdisp_field_r1000_mod))
#p-value = 0.3225
plot(Betdisp_field_r1000_mod)

anova(Betdisp_field_r1000_mod)
#plant_date 0.05566 0.013915     4 8.2672  1.9186 0.1981

emmeans(Betdisp_field_r1000_mod,pairwise~plant_date, adjust="fdr")


#Seedling and Adult only 
Mar_leaf_rain.fung_decon_rar_Mar_field_dis_betadisp_LFE=subset(Mar_leaf_rain.fung_decon_rar_Mar_field_dis_betadisp_plant,plant_type!="Seed")

Betdisp_field_LFE_r1000_mod=lmer(sqrt(betadisp)~plant_type*collect_date+(1|gh_block),data = Mar_leaf_rain.fung_decon_rar_Mar_field_dis_betadisp_LFE)
hist(resid(Betdisp_field_LFE_r1000_mod))
qqPlot(resid(Betdisp_field_LFE_r1000_mod))
shapiro.test(resid(Betdisp_field_LFE_r1000_mod))
#p-value = 0.3298
plot(Betdisp_field_LFE_r1000_mod)

####Table S10 Beta-dispersion (Bray)####
anova(Betdisp_field_LFE_r1000_mod)
#plant_type              0.0033512 0.0033512     1 18.515  0.4968 0.48970  
#collect_date            0.0227004 0.0227004     1 18.286  3.3649 0.08292 .
#plant_type:collect_date 0.0004254 0.0004254     1 18.515  0.0631 0.80448  


####Graph Field experiment Bray NMDS####

#Creation of Hulls for the graphing of Bray NMDS
Mar_leaf_rain.fung_decon_rar_Mar_field_map=sample_data(Mar_leaf_rain.fung_decon_rar_Mar_field)
Mar_leaf_rain.fung_decon_rar_Mar_field_coordinates=merge(Mar_leaf_rain.fung_decon_rar_Mar_field.ord$points,Mar_leaf_rain.fung_decon_rar_Mar_field_map, by="row.names")
Mar_leaf_rain.fung_decon_rar_Mar_field_coordinates$plant_data_grp=ifelse(Mar_leaf_rain.fung_decon_rar_Mar_field_coordinates$plant_type=="Rain",
                                                                         "Rain", ifelse(Mar_leaf_rain.fung_decon_rar_Mar_field_coordinates$plant_type=="Seedling",
                                                                                        ifelse(Mar_leaf_rain.fung_decon_rar_Mar_field_coordinates$collect_date=="9/5/2018","Sept_Seedling","Jul_Seedling"),
                                                                                        ifelse(Mar_leaf_rain.fung_decon_rar_Mar_field_coordinates$plant_type=="Adult",
                                                                                               ifelse(Mar_leaf_rain.fung_decon_rar_Mar_field_coordinates$collect_date=="9/5/2018","Sept_Adult","Jul_Adult"),"Seed")))
grp1000.Adult_Jul <- Mar_leaf_rain.fung_decon_rar_Mar_field_coordinates[Mar_leaf_rain.fung_decon_rar_Mar_field_coordinates$plant_type == "Adult"&
                                                                          Mar_leaf_rain.fung_decon_rar_Mar_field_coordinates$collect_date == "7/16/2018", 
][chull(Mar_leaf_rain.fung_decon_rar_Mar_field_coordinates[
  Mar_leaf_rain.fung_decon_rar_Mar_field_coordinates$
    plant_type =="Adult"&Mar_leaf_rain.fung_decon_rar_Mar_field_coordinates$
    collect_date == "7/16/2018", c("MDS1", "MDS2")]), ] 
nrow(grp1000.Adult_Jul)
#3
grp1000.Adult_Sept <- Mar_leaf_rain.fung_decon_rar_Mar_field_coordinates[Mar_leaf_rain.fung_decon_rar_Mar_field_coordinates$plant_type == "Adult"&
                                                                           Mar_leaf_rain.fung_decon_rar_Mar_field_coordinates$collect_date == "9/5/2018", 
][chull(Mar_leaf_rain.fung_decon_rar_Mar_field_coordinates[
  Mar_leaf_rain.fung_decon_rar_Mar_field_coordinates$
    plant_type =="Adult"&Mar_leaf_rain.fung_decon_rar_Mar_field_coordinates$
    collect_date == "9/5/2018", c("MDS1", "MDS2")]), ] 

nrow(grp1000.Adult_Sept)
#3


grp1000.Seedling_Jul <- Mar_leaf_rain.fung_decon_rar_Mar_field_coordinates[Mar_leaf_rain.fung_decon_rar_Mar_field_coordinates$plant_type == "Seedling"&
                                                                             Mar_leaf_rain.fung_decon_rar_Mar_field_coordinates$collect_date == "7/16/2018", 
][chull(Mar_leaf_rain.fung_decon_rar_Mar_field_coordinates[
  Mar_leaf_rain.fung_decon_rar_Mar_field_coordinates$
    plant_type =="Seedling"&Mar_leaf_rain.fung_decon_rar_Mar_field_coordinates$
    collect_date == "7/16/2018", c("MDS1", "MDS2")]), ] 
nrow(grp1000.Seedling_Jul)
#3
grp1000.Seedling_Sept <- Mar_leaf_rain.fung_decon_rar_Mar_field_coordinates[Mar_leaf_rain.fung_decon_rar_Mar_field_coordinates$plant_type == "Seedling"&
                                                                              Mar_leaf_rain.fung_decon_rar_Mar_field_coordinates$collect_date == "9/5/2018", 
][chull(Mar_leaf_rain.fung_decon_rar_Mar_field_coordinates[
  Mar_leaf_rain.fung_decon_rar_Mar_field_coordinates$
    plant_type =="Seedling"&Mar_leaf_rain.fung_decon_rar_Mar_field_coordinates$
    collect_date == "9/5/2018", c("MDS1", "MDS2")]), ] 

nrow(grp1000.Seedling_Sept)
#5

grp1000.Rain <- Mar_leaf_rain.fung_decon_rar_Mar_field_coordinates[Mar_leaf_rain.fung_decon_rar_Mar_field_coordinates$plant_type == "Rain", 
][chull(Mar_leaf_rain.fung_decon_rar_Mar_field_coordinates[
  Mar_leaf_rain.fung_decon_rar_Mar_field_coordinates$
    plant_type =="Rain", c("MDS1", "MDS2")]), ] 

nrow(grp1000.Rain)
#9

grp1000.Seed <- Mar_leaf_rain.fung_decon_rar_Mar_field_coordinates[Mar_leaf_rain.fung_decon_rar_Mar_field_coordinates$plant_type == "Seed", 
][chull(Mar_leaf_rain.fung_decon_rar_Mar_field_coordinates[
  Mar_leaf_rain.fung_decon_rar_Mar_field_coordinates$
    plant_type =="Seed", c("MDS1", "MDS2")]), ] 

nrow(grp1000.Seed)
#4


hull1000_Field=rbind(grp1000.Adult_Jul,grp1000.Adult_Sept,grp1000.Seedling_Jul,grp1000.Seedling_Sept,grp1000.Rain, grp1000.Seed)
nrow(hull1000_Field)
#27
Mar_leaf_rain.fung_decon_rar_Mar_field_coor_rain=subset(Mar_leaf_rain.fung_decon_rar_Mar_field_coordinates, plant_type=="Rain")
nrow(Mar_leaf_rain.fung_decon_rar_Mar_field_coor_rain)
#56
Mar_leaf_rain.fung_decon_rar_Mar_field_coor_plant=subset(Mar_leaf_rain.fung_decon_rar_Mar_field_coordinates, plant_type!="Rain"&
                                                           plant_type!="Seed")
nrow(Mar_leaf_rain.fung_decon_rar_Mar_field_coor_plant)
#25
Mar_leaf_rain.fung_decon_rar_Mar_field_coor_seed=subset(Mar_leaf_rain.fung_decon_rar_Mar_field_coordinates,plant_type=="Seed")
nrow(Mar_leaf_rain.fung_decon_rar_Mar_field_coor_seed)
#4


plant_data_grp_order=c("Seed","Jul_Seedling","Sept_Seedling","Jul_Adult","Sept_Adult","Rain")


ggplot(Mar_leaf_rain.fung_decon_rar_Mar_field_coordinates,aes(x= MDS1, y=MDS2))+
  geom_polygon(data=hull1000_Field,aes(x=MDS1,y=MDS2,fill=factor(plant_data_grp,levels = plant_data_grp_order),
                                       group=factor(plant_data_grp,levels = plant_data_grp_order)), alpha=0.5)+
  geom_point(size=4,aes(fill=factor(plant_data_grp,levels = plant_data_grp_order),
                        shape=factor(plant_type,levels = c("Adult","Seed","Seedling","Rain")), color=factor(plant_data_grp,levels = plant_data_grp_order)),stroke=2)+
  scale_shape_manual(values = c(22,13,24,21),name=NULL)+xlab("NMDS1")+ylab("NMDS2")+
  scale_color_manual(values = c("#E31A1C","black", "black","black","black","black"),name=NULL)+
  scale_fill_manual(values = c("#E31A1C","#B2DF8A","#33A02C","#A6CEE3","#1F78B4","black"),name=NULL)+
  theme_bw()+theme(axis.title = element_text(size = 30),axis.text =  element_text(size = 30),legend.position = "none")


#Leaf55 and Leaf25 are outliers
#900x800

#####Stacked bar graph for the Taxonomy####

#####Plant Stack Graphs####


#Remove the Rain for the this graph
Mar_leaf_rain.fung_decon_rar_Mar_field_NR=subset_samples( Mar_leaf_rain.fung_decon_rar_Mar_field,plant_type!="Rain")
nsamples(Mar_leaf_rain.fung_decon_rar_Mar_field_NR)
#29
Mar_leaf_rain.fung_decon_rar_Mar_field_NR=prune_taxa(taxa_sums(Mar_leaf_rain.fung_decon_rar_Mar_field_NR) > 0, Mar_leaf_rain.fung_decon_rar_Mar_field_NR)
#Stack bar graphs

sample_data(Mar_leaf_rain.fung_decon_rar_Mar_field_NR)$date_plant=droplevels(with(sample_data(Mar_leaf_rain.fung_decon_rar_Mar_field_NR), interaction(collect_date, plant_type)))


unique(sample_data(Mar_leaf_rain.fung_decon_rar_Mar_field_NR)$date_plant)

ntaxa(Mar_leaf_rain.fung_decon_rar_Mar_field_NR)
#653
#merge OTUs by the soil and precipitation treatment
Mar_leaf_rain.fung_decon_rar_Mar_field_NR_fact=merge_samples(Mar_leaf_rain.fung_decon_rar_Mar_field_NR, "date_plant")
sample_names(Mar_leaf_rain.fung_decon_rar_Mar_field_NR_fact)     

#combine the reads at Phylum level
get_taxa_unique(Mar_leaf_rain.fung_decon_rar_Mar_field_NR_fact, taxonomic.rank="Phylum")
#"Basidiomycota" "Unknown"       "Ascomycota" 
(Mar_leaf_rain.fung_decon_rar_Mar_field_NR_fact.phylum<-tax_glom(Mar_leaf_rain.fung_decon_rar_Mar_field_NR_fact, taxrank="Phylum"))



#subset so there is only the top ten most abundant phyla

Mar_leaf_rain.fung_decon_rar_Mar_field_NR_Names=(data.frame(tax_table(Mar_leaf_rain.fung_decon_rar_Mar_field_NR_fact.phylum))$Phylum)

#Transform the read counts to prop of total reads

Mar_leaf_rain.fung_decon_rar_Mar_field_NR_fact.phylum.prop=transform_sample_counts(Mar_leaf_rain.fung_decon_rar_Mar_field_NR_fact.phylum, function(x)x/sum(x))

taxon_positions=c("7/10/2018.Seed", "7/16/2018.Seedling","9/5/2018.Seedling","7/16/2018.Adult", "9/5/2018.Adult")


Mar_leaf_rain.fung_decon_rar_Mar_field_NR_fact.phylum.prop_otu=as.data.frame(t(otu_table(Mar_leaf_rain.fung_decon_rar_Mar_field_NR_fact.phylum.prop)))
Mar_leaf_rain.fung_decon_rar_Mar_field_NR_fact.phylum.prop_otu[,"Phylum"]=Mar_leaf_rain.fung_decon_rar_Mar_field_NR_Names
Mar_leaf_rain.fung_decon_rar_Mar_field_NR_fact.phylum.prop_otu_M=melt(Mar_leaf_rain.fung_decon_rar_Mar_field_NR_fact.phylum.prop_otu,id="Phylum")



(p_Field_plants_color=ggplot(Mar_leaf_rain.fung_decon_rar_Mar_field_NR_fact.phylum.prop_otu_M,aes(x=variable,y=value,fill=Phylum))+
    geom_bar(aes( fill=Phylum), stat="identity", position="stack",color="black")+theme_bw()+
    theme(axis.text.y=element_text(size=18),axis.text.x=element_text(size=18),
          axis.title=element_text(size=20),panel.grid.major=element_blank(),legend.text = element_text(size=16), legend.title = element_text(size=20),
          panel.grid.minor=element_blank())+xlab(NULL)+ylab("Proportion")+
    scale_x_discrete(limits = taxon_positions,labels=c("Seed","Start\nSeedling","End\nSeedling", "Start\nAdult","End\nAdult"))+
    scale_fill_brewer(palette="Paired")+
    guides(fill=guide_legend(title="Phyla")))

#1000x700




#class Level
#combine the reads at Class level
get_taxa_unique(Mar_leaf_rain.fung_decon_rar_Mar_field_NR_fact, taxonomic.rank="Class")
# [1] "Tremellomycetes"      "Unknown"              "Dothideomycetes"      "Taphrinomycetes"      "Sordariomycetes"      "Exobasidiomycetes"    "Microbotryomycetes"   "Cystobasidiomycetes" 
#[9] "Agaricomycetes"       "Eurotiomycetes"       "Ustilaginomycetes"    "Spiculogloeomycetes"  "Wallemiomycetes"      "Pucciniomycetes"      "Orbiliomycetes"       "Leotiomycetes"       
#[17] "Saccharomycetes"      "Agaricostilbomycetes"
(Mar_leaf_rain.fung_decon_rar_Mar_field_NR_fact.class<-tax_glom(Mar_leaf_rain.fung_decon_rar_Mar_field_NR_fact, taxrank="Class"))
#20 taxa

#How many 

Mar_leaf_rain.fung_decon_rar_Mar_field_NR_class_Names_raw=
  data.frame("Phylum"=data.frame(tax_table(Mar_leaf_rain.fung_decon_rar_Mar_field_NR_fact.class))$Phylum,
             "Class"=data.frame(tax_table(Mar_leaf_rain.fung_decon_rar_Mar_field_NR_fact.class))$Class,
             "OTU"=row.names(data.frame(tax_table(Mar_leaf_rain.fung_decon_rar_Mar_field_NR_fact.class))))

Mar_leaf_rain.fung_decon_rar_Mar_field_NR_class_Names_raw$P_C=with(Mar_leaf_rain.fung_decon_rar_Mar_field_NR_class_Names_raw,interaction(Phylum,Class))
#Let's take the top 10 taxa 

TopCLASS_Field = names(sort(taxa_sums(Mar_leaf_rain.fung_decon_rar_Mar_field_NR_fact.class), TRUE)[1:10])
TOP10_Class_Names=data.frame("Phylum"=data.frame(tax_table(prune_taxa(TopCLASS_Field, Mar_leaf_rain.fung_decon_rar_Mar_field_NR_fact.class)))$Phylum,
                             "Class"=data.frame(tax_table(prune_taxa(TopCLASS_Field, Mar_leaf_rain.fung_decon_rar_Mar_field_NR_fact.class)))$Class,
                             "OTU"=row.names(data.frame(tax_table(prune_taxa(TopCLASS_Field, Mar_leaf_rain.fung_decon_rar_Mar_field_NR_fact.class)))))



#Need to make an other section

Mar_leaf_rain.fung_decon_rar_Mar_field_NR_class_Names=merge(Mar_leaf_rain.fung_decon_rar_Mar_field_NR_class_Names_raw,TOP10_Class_Names, by="OTU", all.x = T)

Mar_leaf_rain.fung_decon_rar_Mar_field_NR_class_Names$Real_P_C=ifelse(is.na(Mar_leaf_rain.fung_decon_rar_Mar_field_NR_class_Names$Class.y),
                                                                      paste(Mar_leaf_rain.fung_decon_rar_Mar_field_NR_class_Names$Phylum.x,
                                                                            "Other_spp",sep="."),
                                                                      ifelse(Mar_leaf_rain.fung_decon_rar_Mar_field_NR_class_Names$Class.x=="Unknown",
                                                                             paste(Mar_leaf_rain.fung_decon_rar_Mar_field_NR_class_Names$Phylum.x,
                                                                                   "Other_spp",sep="."),                              
                                                                             as.character(Mar_leaf_rain.fung_decon_rar_Mar_field_NR_class_Names$P_C)
                                                                      ))

#Transform the read counts to prop of total reads

Mar_leaf_rain.fung_decon_rar_Mar_field_NR_fact.class.prop=transform_sample_counts(Mar_leaf_rain.fung_decon_rar_Mar_field_NR_fact.class, function(x)x/sum(x))




Mar_leaf_rain.fung_decon_rar_Mar_field_NR_fact.class.prop_otu=as.data.frame(t(otu_table(Mar_leaf_rain.fung_decon_rar_Mar_field_NR_fact.class.prop)))
Mar_leaf_rain.fung_decon_rar_Mar_field_NR_fact.class.prop_otu2=merge(Mar_leaf_rain.fung_decon_rar_Mar_field_NR_class_Names[,c("Real_P_C","OTU")],
                                                                     Mar_leaf_rain.fung_decon_rar_Mar_field_NR_fact.class.prop_otu,by.y = "row.names",by.x = "OTU")
summary(Mar_leaf_rain.fung_decon_rar_Mar_field_NR_fact.class.prop_otu2)
Mar_leaf_rain.fung_decon_rar_Mar_field_NR_fact.class.prop_otu2$OTU=NULL
Mar_leaf_rain.fung_decon_rar_Mar_field_NR_fact.class.prop_otu_sum=Mar_leaf_rain.fung_decon_rar_Mar_field_NR_fact.class.prop_otu2%>%group_by(Real_P_C)%>%
  summarise_all(~sum(.))

Mar_leaf_rain.fung_decon_rar_Mar_field_NR_fact.class.prop_otu_sum_M=melt(Mar_leaf_rain.fung_decon_rar_Mar_field_NR_fact.class.prop_otu_sum,id="Real_P_C")
unique(Mar_leaf_rain.fung_decon_rar_Mar_field_NR_fact.class.prop_otu_sum_M$Real_P_C)
taxon_positions=c("7/10/2018.Seed", "7/16/2018.Seedling","9/5/2018.Seedling","7/16/2018.Adult", "9/5/2018.Adult")
fill_order=c("Ascomycota.Dothideomycetes","Ascomycota.Eurotiomycetes","Ascomycota.Sordariomycetes","Ascomycota.Other_spp","Basidiomycota.Agaricomycetes",
             "Basidiomycota.Cystobasidiomycetes","Basidiomycota.Exobasidiomycetes","Basidiomycota.Microbotryomycetes","Basidiomycota.Tremellomycetes",
             "Basidiomycota.Ustilaginomycetes","Basidiomycota.Other_spp","Unknown.Other_spp" )

(p_class_Field_plants_color=ggplot(Mar_leaf_rain.fung_decon_rar_Mar_field_NR_fact.class.prop_otu_sum_M,aes(x=factor(variable,levels = taxon_positions),
                                                                                                           y=value,fill=factor(Real_P_C,levels = fill_order)))+
    geom_bar(aes( fill=factor(Real_P_C,levels = fill_order)), stat="identity", position="stack",color="black")+theme_bw()+
    theme(axis.text.y=element_text(size=18),axis.text.x=element_text(size=18),
          axis.title=element_text(size=20),panel.grid.major=element_blank(),legend.text = element_text(size=16), legend.title = element_text(size=20),
          panel.grid.minor=element_blank())+xlab(NULL)+ylab("Proportion")+
    scale_x_discrete(limits = taxon_positions,labels=c("Seed","Start\nSeedling","End\nSeedling", "Start\nAdult","End\nAdult"))+scale_fill_brewer(palette="Paired")+
    guides(fill=guide_legend(title="Phyla")))

#1200x700



#####PLOT: Figure S12####
#####Plot Field Stacked Taxa####

(p_Field_plants_color=ggplot(Mar_leaf_rain.fung_decon_rar_Mar_field_NR_fact.phylum.prop_otu_M,aes(x=variable,y=value,fill=Phylum))+
   geom_bar(aes( fill=Phylum), stat="identity", position="stack",color="black")+theme_cowplot()+
   theme(axis.text.y=element_text(size=24),axis.text.x=element_blank(),
         axis.title=element_text(size=40),panel.grid.major=element_blank(),legend.position = "none",
         panel.grid.minor=element_blank())+xlab(NULL)+scale_y_continuous(name = "Proportion",limits = c(0,1.1),breaks = c(seq(0,1,by=0.2)))+
   scale_x_discrete(limits = taxon_positions,labels=c("Seed","Start\nSeedling","End\nSeedling", "Start\nAdult","End\nAdult"))+
   scale_fill_brewer(palette="Paired")+
   guides(fill=guide_legend(title="Phyla")))

(p_class_Field_plants_color=ggplot(Mar_leaf_rain.fung_decon_rar_Mar_field_NR_fact.class.prop_otu_sum_M,aes(x=factor(variable,levels = taxon_positions),
                                                                                                           y=value,fill=factor(Real_P_C,levels = fill_order)))+
    geom_bar(aes( fill=factor(Real_P_C,levels = fill_order)), stat="identity", position="stack",color="black")+theme_cowplot()+
    theme(axis.text=element_text(size=24),
          axis.title=element_text(size=40),panel.grid.major=element_blank(),legend.position = "none",
          panel.grid.minor=element_blank())+xlab(NULL)+scale_y_continuous(name = "Proportion",limits = c(0,1.1),breaks = c(seq(0,1,by=0.2)))+
    scale_x_discrete(limits = taxon_positions,labels=c("Seed","Start\nSeedling","End\nSeedling", "Start\nAdult","End\nAdult"))+
    scale_fill_brewer(palette="Paired")+
    guides(fill=guide_legend(title="Phyla")))


plot_grid(p_Field_plants_color,p_class_Field_plants_color,nrow = 2,rel_heights = c(0.8,1))
#700x1000
#####Rain Stack Graphs####
#Let's graph rain 

#Remove the Rain for the this graph
Mar_rain.fung_decon_rar_Mar_field=subset_samples( Mar_leaf_rain.fung_decon_rar_Mar_field,plant_type=="Rain")
nsamples(Mar_rain.fung_decon_rar_Mar_field)
#56
Mar_rain.fung_decon_rar_Mar_field=prune_taxa(taxa_sums(Mar_rain.fung_decon_rar_Mar_field) > 0, Mar_rain.fung_decon_rar_Mar_field)
#Stack bar graphs
length(unique(sample_data(Mar_rain.fung_decon_rar_Mar_field)$collect_date))
#15
ntaxa(Mar_rain.fung_decon_rar_Mar_field)
#1907
#merge OTUs by the soil and precipitation treatment
Mar_rain.fung_decon_rar_Mar_field_fact=merge_samples(Mar_rain.fung_decon_rar_Mar_field, "collect_date")
sample_names(Mar_rain.fung_decon_rar_Mar_field_fact)     

#combine the reads at Phylum level
get_taxa_unique(Mar_rain.fung_decon_rar_Mar_field_fact, taxonomic.rank="Phylum")
#"Basidiomycota"   "Ascomycota"      "Unknown"         "Chytridiomycota" "Mucoromycota"    "Rozellomycota" 
(Mar_rain.fung_decon_rar_Mar_field_fact.phylum<-tax_glom(Mar_rain.fung_decon_rar_Mar_field_fact, taxrank="Phylum"))



#subset so there is only the top ten most abundant phyla

Mar_rain.fung_decon_rar_Mar_field_Names=(data.frame(tax_table(Mar_rain.fung_decon_rar_Mar_field_fact.phylum))$Phylum)

#Transform the read counts to prop of total reads

Mar_rain.fung_decon_rar_Mar_field_fact.phylum.prop=transform_sample_counts(Mar_rain.fung_decon_rar_Mar_field_fact.phylum, function(x)x/sum(x))




Mar_rain.fung_decon_rar_Mar_field_fact.phylum.prop_otu=as.data.frame(t(otu_table(Mar_rain.fung_decon_rar_Mar_field_fact.phylum.prop)))
Mar_rain.fung_decon_rar_Mar_field_fact.phylum.prop_otu[,"Phylum"]=Mar_rain.fung_decon_rar_Mar_field_Names
Mar_rain.fung_decon_rar_Mar_field_fact.phylum.prop_otu_M=melt(Mar_rain.fung_decon_rar_Mar_field_fact.phylum.prop_otu,id="Phylum")

#####PLOT: Figure S14####

(p_Rain_color=ggplot(Mar_rain.fung_decon_rar_Mar_field_fact.phylum.prop_otu_M,aes(x=as.Date(variable,format="%m/%d/%Y"),y=value))+
    geom_area(aes( fill=Phylum))+theme_bw()+
    theme(axis.text.y=element_text(size=18),axis.text.x=element_text(size=18),
          axis.title=element_text(size=20),panel.grid.major=element_blank(),legend.text = element_text(size=16), legend.title = element_text(size=20),
          panel.grid.minor=element_blank())+xlab(NULL)+ylab("Proportion")+scale_fill_brewer(palette="Paired")+
    guides(fill=guide_legend(title="Phyla"))+ 
    scale_x_date(date_labels = "%b %d",date_breaks = "1 week"))

#1000x700

(p_Rain_color=ggplot(Mar_rain.fung_decon_rar_Mar_field_fact.phylum.prop_otu_M,aes(x=as.Date(variable,format="%m/%d/%Y"),y=value))+
    geom_area(aes( fill=Phylum))+theme_bw()+
    theme(axis.text.y=element_text(size=18),axis.text.x=element_text(size=18),
          axis.title=element_text(size=20),panel.grid.major=element_blank(),legend.position = "none",
          panel.grid.minor=element_blank())+xlab(NULL)+ylab("Proportion")+scale_fill_brewer(palette="Paired")+
    guides(fill=guide_legend(title="Phyla"))+ 
    scale_x_date(date_labels = "%b %d",date_breaks = "1 week"))

#1000x700

#class Level
#combine the reads at Class level
get_taxa_unique(Mar_rain.fung_decon_rar_Mar_field_fact, taxonomic.rank="Class")
#[1] "Tremellomycetes"      "Dothideomycetes"      "Taphrinomycetes"      "Sordariomycetes"      "Exobasidiomycetes"    "Microbotryomycetes"   "Cystobasidiomycetes"  "Unknown"             
#[9] "Agaricomycetes"       "Lecanoromycetes"      "Ustilaginomycetes"    "Leotiomycetes"        "Agaricostilbomycetes" "Eurotiomycetes"       "Spiculogloeomycetes"  "Classiculomycetes"   
#[17] "Saccharomycetes"      "Pucciniomycetes"      "Glomeromycetes"       "Orbiliomycetes"       "Mucoromycetes"        "Mortierellomycetes"   "Atractiellomycetes"   "Dacrymycetes"        
#[25] "Laboulbeniomycetes" 
(Mar_rain.fung_decon_rar_Mar_field_fact.class<-tax_glom(Mar_rain.fung_decon_rar_Mar_field_fact, taxrank="Class"))
#27 taxa



Mar_rain.fung_decon_rar_Mar_field_fact.class_Names_raw=data.frame("Phylum"=data.frame(tax_table(Mar_rain.fung_decon_rar_Mar_field_fact.class))$Phylum,
                                                                  "Class"=data.frame(tax_table(Mar_rain.fung_decon_rar_Mar_field_fact.class))$Class,
                                                                  "OTU"=taxa_names(Mar_rain.fung_decon_rar_Mar_field_fact.class))

Mar_rain.fung_decon_rar_Mar_field_fact.class_Names_raw$P_C=with(Mar_rain.fung_decon_rar_Mar_field_fact.class_Names_raw,interaction(Phylum,Class))
#Let's take the top 10 taxa 

R_TopCLASS_Field = names(sort(taxa_sums(Mar_rain.fung_decon_rar_Mar_field_fact.class), TRUE)[1:9])
R_TOP10_Class_Names=data.frame("Phylum"=data.frame(tax_table(prune_taxa(R_TopCLASS_Field, Mar_rain.fung_decon_rar_Mar_field_fact.class)))$Phylum,
                               "Class"=data.frame(tax_table(prune_taxa(R_TopCLASS_Field, Mar_rain.fung_decon_rar_Mar_field_fact.class)))$Class,
                               "OTU"=R_TopCLASS_Field)



#Need to make an other section

Mar_rain.fung_decon_rar_Mar_field_fact.class_Names=merge(Mar_rain.fung_decon_rar_Mar_field_fact.class_Names_raw,R_TOP10_Class_Names, by="OTU", all.x = T)

Mar_rain.fung_decon_rar_Mar_field_fact.class_Names$Real_P_C=ifelse(is.na(Mar_rain.fung_decon_rar_Mar_field_fact.class_Names$Class.y),
                                                                   paste(Mar_rain.fung_decon_rar_Mar_field_fact.class_Names$Phylum.x,
                                                                         "Other_spp",sep="."),
                                                                   ifelse(Mar_rain.fung_decon_rar_Mar_field_fact.class_Names$Class.x=="Unknown",
                                                                          paste(Mar_rain.fung_decon_rar_Mar_field_fact.class_Names$Phylum.x,
                                                                                "Other_spp",sep="."),                              
                                                                          as.character(Mar_rain.fung_decon_rar_Mar_field_fact.class_Names$P_C)
                                                                   ))

#Transform the read counts to prop of total reads

Mar_rain.fung_decon_rar_Mar_field_fact.class.prop=transform_sample_counts(Mar_rain.fung_decon_rar_Mar_field_fact.class, function(x)x/sum(x))




Mar_rain.fung_decon_rar_Mar_field_fact.class.prop_otu=as.data.frame(t(otu_table(Mar_rain.fung_decon_rar_Mar_field_fact.class.prop)))
Mar_rain.fung_decon_rar_Mar_field_fact.class.prop_otu2=merge(Mar_rain.fung_decon_rar_Mar_field_fact.class_Names[,c("Real_P_C","OTU")],
                                                             Mar_rain.fung_decon_rar_Mar_field_fact.class.prop_otu,by.y = "row.names",by.x = "OTU")

Mar_rain.fung_decon_rar_Mar_field_fact.class.prop_otu2$OTU=NULL
Mar_rain.fung_decon_rar_Mar_field_fact.class.prop_otu2_sum=Mar_rain.fung_decon_rar_Mar_field_fact.class.prop_otu2%>%group_by(Real_P_C)%>%
  summarise_all(~sum(.))
rowSums(Mar_rain.fung_decon_rar_Mar_field_fact.class.prop_otu2_sum[,2:ncol(Mar_rain.fung_decon_rar_Mar_field_fact.class.prop_otu2_sum)])
Mar_rain.fung_decon_rar_Mar_field_fact.class.prop_otu2_sum_trunc=
  Mar_rain.fung_decon_rar_Mar_field_fact.class.prop_otu2_sum[order(-rowSums(
    Mar_rain.fung_decon_rar_Mar_field_fact.class.prop_otu2_sum[,2:ncol(Mar_rain.fung_decon_rar_Mar_field_fact.class.prop_otu2_sum)])),][1:11,]

Mar_rain.fung_decon_rar_Mar_field_fact.class.prop_otu2_sum_trunc1=Mar_rain.fung_decon_rar_Mar_field_fact.class.prop_otu2_sum_trunc %>% 
  bind_rows(summarise_all(., list(~if(is.numeric(.)) 1-sum(.) else "Other_Phyla")))
summary(Mar_rain.fung_decon_rar_Mar_field_fact.class.prop_otu2_sum_trunc1)
Mar_rain.fung_decon_rar_Mar_field_fact.class.prop_otu2_sum_M=melt(Mar_rain.fung_decon_rar_Mar_field_fact.class.prop_otu2_sum_trunc1,id="Real_P_C")
unique(Mar_rain.fung_decon_rar_Mar_field_fact.class.prop_otu2_sum_M$Real_P_C)

fill_order_R=c("Ascomycota.Dothideomycetes","Ascomycota.Sordariomycetes","Ascomycota.Other_spp","Basidiomycota.Agaricomycetes","Basidiomycota.Cystobasidiomycetes",
               "Basidiomycota.Exobasidiomycetes","Basidiomycota.Microbotryomycetes","Basidiomycota.Tremellomycetes","Basidiomycota.Ustilaginomycetes","Basidiomycota.Other_spp",
               "Other_Phyla","Unknown.Other_spp")
(p_CLASS_Rain_color=ggplot(Mar_rain.fung_decon_rar_Mar_field_fact.class.prop_otu2_sum_M,aes(x=as.Date(variable,format="%m/%d/%Y"),y=value))+
    geom_area(aes( fill=factor(Real_P_C,levels = fill_order_R)))+theme_bw()+
    theme(axis.text.y=element_text(size=18),axis.text.x=element_text(size=18),
          axis.title=element_text(size=20),panel.grid.major=element_blank(),legend.text = element_text(size=16), legend.title = element_text(size=20),
          panel.grid.minor=element_blank())+xlab(NULL)+ylab("Proportion")+scale_fill_brewer(palette="Paired")+
    guides(fill=guide_legend(title="Phyla"))+ 
    scale_x_date(date_labels = "%b %d",date_breaks = "1 week"))

#1300x700


(p_CLASS_Rain_color=ggplot(Mar_rain.fung_decon_rar_Mar_field_fact.class.prop_otu2_sum_M,aes(x=as.Date(variable,format="%m/%d/%Y"),y=value))+
    geom_area(aes( fill=factor(Real_P_C,levels = fill_order_R)))+theme_bw()+
    theme(axis.text.y=element_text(size=18),axis.text.x=element_text(size=18),
          axis.title=element_text(size=20),panel.grid.major=element_blank(),legend.position = "none",
          panel.grid.minor=element_blank())+xlab(NULL)+ylab("Proportion")+scale_fill_brewer(palette="Paired")+
    guides(fill=guide_legend(title="Phyla"))+ 
    scale_x_date(date_labels = "%b %d",date_breaks = "1 week"))

#1000x700





#####Field and Rain  Diversity#####
load(here::here("R_file","Mar_leaf_rain_ZOTU_full_v4.2.2020.fung_decon_rar_pruned_phyloseq_obj.RData"))
Mar_leaf_rain.fung_decon_rar_Mar=subset_samples(Mar_leaf_rain.fung_decon_rar, sub_proj=="Marshall")
nsamples(Mar_leaf_rain.fung_decon_rar_Mar)
#117
Mar_leaf_rain.fung_decon_rar_Mar_field=subset_samples(Mar_leaf_rain.fung_decon_rar_Mar, plant_type!="Petri")
nsamples(Mar_leaf_rain.fung_decon_rar_Mar_field)
#85

unique(sample_data(Mar_leaf_rain.fung_decon_rar_Mar_field)$sub_proj)

Mar_leaf_rain.fung_decon_rar_Mar_field=prune_taxa(taxa_sums(Mar_leaf_rain.fung_decon_rar_Mar_field) > 0, Mar_leaf_rain.fung_decon_rar_Mar_field)
ntaxa(Mar_leaf_rain.fung_decon_rar_Mar_field)
#2243
sum(taxa_sums(Mar_leaf_rain.fung_decon_rar_Mar_field))
#85000
min(sample_sums(Mar_leaf_rain.fung_decon_rar_Mar_field))
#1000
max(sample_sums(Mar_leaf_rain.fung_decon_rar_Mar_field))
#1000

Mar_leaf_rain.fung_decon_rar_Mar_field_map=sample_data(Mar_leaf_rain.fung_decon_rar_Mar_field)
Mar_leaf_rain.fung_decon_rar_Mar_field_div=estimate_richness(Mar_leaf_rain.fung_decon_rar_Mar_field,measures=c("Observed","InvSimpson","Shannon"))

Mar_leaf_rain.fung_decon_rar_Mar_field_div=merge(Mar_leaf_rain.fung_decon_rar_Mar_field_div,Mar_leaf_rain.fung_decon_rar_Mar_field_map, by="row.names")
head(Mar_leaf_rain.fung_decon_rar_Mar_field_div)

#####PLOT: Figure S11####
#####Graphing Field Diversity####

Mar_leaf_rain.fung_decon_rar_Mar_field_div_sum=Mar_leaf_rain.fung_decon_rar_Mar_field_div%>%group_by(collect_date,plant_type)%>%
  summarise_at(vars("Observed","InvSimpson"),list(~mean(.),~n(),se=~sd(.)/sqrt(n())))

with(Mar_leaf_rain.fung_decon_rar_Mar_field_div_sum,interaction(collect_date,plant_type))

plant_date_order=c("7/10/2018.Seed","7/16/2018.Adult","9/5/2018.Adult","7/16/2018.Seedling","9/5/2018.Seedling")

#Richness
#plants
(field_rich=ggplot(subset(Mar_leaf_rain.fung_decon_rar_Mar_field_div,plant_type!="Rain"))+
    geom_point(aes(x=factor(interaction(collect_date,plant_type),levels = plant_date_order,
                            labels = c("Seed", "Start","End","Start","End")),y=Observed, shape=plant_type,
                   fill=factor(interaction(collect_date,plant_type),
                               levels = plant_date_order),
                   color=factor(interaction(collect_date,plant_type),
                                levels = plant_date_order)),size=2, alpha = 0.5)+
    geom_errorbar(data = subset(Mar_leaf_rain.fung_decon_rar_Mar_field_div_sum,plant_type!="Rain"),
                  aes(x=factor(interaction(collect_date,plant_type),levels = plant_date_order,
                               labels = c("Seed", "Start","End","Start","End")),
                      y=Observed_mean, ymax=Observed_mean+Observed_se,
                      ymin=Observed_mean-Observed_se), width=0.3)+
    geom_point(data = subset(Mar_leaf_rain.fung_decon_rar_Mar_field_div_sum,plant_type!="Rain"),
               aes(x=factor(interaction(collect_date,plant_type),
                            levels = plant_date_order,
                            labels = c("Seed", "Start","End","Start","End")),
                   y=Observed_mean, shape=plant_type,
                   fill=factor(interaction(collect_date,plant_type),
                               levels = plant_date_order),
                   color=factor(interaction(collect_date,plant_type),
                                levels = plant_date_order)),size=10)+scale_x_discrete(labels=)+
    scale_shape_manual(values = c(22,13,24),name=NULL)+ylab("Richness (ZOTUs)")+
    facet_grid(cols = vars(factor(plant_type,levels = c("Seed","Adult","Seedling"))),scales = "free_x",space = "free_x")+
    scale_color_manual(values = c("#E31A1C","black", "black","black","black"),name=NULL)+
    scale_fill_manual(values = c("#E31A1C","#A6CEE3","#1F78B4","#B2DF8A","#33A02C"),name=NULL)+
    theme_cowplot()+theme(axis.title.y = element_text(size = 34),axis.title.x = element_blank(),axis.text =  element_text(size = 30),
                          legend.position = "none",strip.background = element_rect(fill = NA),strip.text = element_text(size = 34)))

#1000*800

#Increase in seedling richness

subset(Mar_leaf_rain.fung_decon_rar_Mar_field_div_sum,plant_type=="Seedling")[,c("collect_date","Observed_mean")]
(103-50)/50


#Rain
(rain_rich=ggplot(subset(Mar_leaf_rain.fung_decon_rar_Mar_field_div,plant_type=="Rain"))+geom_point(aes(x=as.Date(collect_date,format="%m/%d/%Y"),
                                                                                                        y=Observed),size=2, alpha = 0.5, shape=21, fill="black")+
    geom_errorbar(data = subset(Mar_leaf_rain.fung_decon_rar_Mar_field_div_sum,plant_type=="Rain"),aes(x=as.Date(collect_date,format="%m/%d/%Y"),
                                                                                                       y=Observed_mean, ymax=Observed_mean+Observed_se,
                                                                                                       ymin=Observed_mean-Observed_se), width=1.5)+
    geom_point(data = subset(Mar_leaf_rain.fung_decon_rar_Mar_field_div_sum,plant_type=="Rain"),aes(x=as.Date(collect_date,format="%m/%d/%Y"),
                                                                                                    y=Observed_mean),size=10,shape=21, fill="black")+
    theme_cowplot()+scale_x_date(date_labels = "%b %d",date_breaks = "11 days")+ylab("Richness (ZOTUs)")+facet_grid(cols = vars(factor(plant_type,labels = "")))+
    theme(axis.title.y = element_text(size = 34),axis.title.x = element_blank(),axis.text =  element_text(size = 30),
          legend.position = "none",strip.background = element_rect(fill = NA),strip.text = element_text(size = 34)))

#1300*800
#Rich_Rain_Field_rar1000_dot_error

#InvSimpson
#Plant
(field_invSimp=ggplot(subset(Mar_leaf_rain.fung_decon_rar_Mar_field_div,plant_type!="Rain"))+
    geom_point(aes(x=factor(interaction(collect_date,plant_type),levels = plant_date_order,
                            labels = c("Seed", "Start","End","Start","End")),y=InvSimpson, shape=plant_type,
                   fill=factor(interaction(collect_date,plant_type),
                               levels = plant_date_order),
                   color=factor(interaction(collect_date,plant_type),
                                levels = plant_date_order)),size=2, alpha = 0.5)+
    geom_errorbar(data = subset(Mar_leaf_rain.fung_decon_rar_Mar_field_div_sum,plant_type!="Rain"),
                  aes(x=factor(interaction(collect_date,plant_type),levels = plant_date_order,
                               labels = c("Seed", "Start","End","Start","End")),
                      y=InvSimpson_mean, ymax=InvSimpson_mean+InvSimpson_se,
                      ymin=InvSimpson_mean-InvSimpson_se), width=0.3)+
    geom_point(data = subset(Mar_leaf_rain.fung_decon_rar_Mar_field_div_sum,plant_type!="Rain"),
               aes(x=factor(interaction(collect_date,plant_type),
                            levels = plant_date_order,
                            labels = c("Seed", "Start","End","Start","End")),
                   y=InvSimpson_mean, shape=plant_type,
                   fill=factor(interaction(collect_date,plant_type),
                               levels = plant_date_order),
                   color=factor(interaction(collect_date,plant_type),
                                levels = plant_date_order)),size=10)+scale_x_discrete(labels=)+
    scale_shape_manual(values = c(22,13,24),name=NULL)+ylab("inverse Simpson")+
    facet_grid(cols = vars(factor(plant_type,levels = c("Seed","Adult","Seedling"))),scales = "free_x",space = "free_x")+
    scale_color_manual(values = c("#E31A1C","black", "black","black","black"),name=NULL)+
    scale_fill_manual(values = c("#E31A1C","#A6CEE3","#1F78B4","#B2DF8A","#33A02C"),name=NULL)+
    theme_cowplot()+theme(axis.title.y = element_text(size = 34),axis.title.x = element_blank(),axis.text =  element_text(size = 30),
                          legend.position = "none",strip.background = element_rect(fill = NA),strip.text = element_text(size = 34)))

#1000x800
#Rain
(rain_invSimp=ggplot(subset(Mar_leaf_rain.fung_decon_rar_Mar_field_div,plant_type=="Rain"))+geom_point(aes(x=as.Date(collect_date,format="%m/%d/%Y"),
                                                                                                           y=InvSimpson),size=2, alpha = 0.5, shape=21, fill="black")+
    geom_errorbar(data = subset(Mar_leaf_rain.fung_decon_rar_Mar_field_div_sum,plant_type=="Rain"),aes(x=as.Date(collect_date,format="%m/%d/%Y"),
                                                                                                       y=InvSimpson_mean, ymax=InvSimpson_mean+InvSimpson_se,
                                                                                                       ymin=InvSimpson_mean-InvSimpson_se), width=1.5)+
    geom_point(data = subset(Mar_leaf_rain.fung_decon_rar_Mar_field_div_sum,plant_type=="Rain"),aes(x=as.Date(collect_date,format="%m/%d/%Y"),
                                                                                                    y=InvSimpson_mean),size=10,shape=21, fill="black")+
    theme_cowplot()+scale_x_date(date_labels = "%b %d",date_breaks = "11 days")+ylab("inverse Simpson")+facet_grid(cols = vars(factor(plant_type,labels = "")))+
    theme(axis.title.y = element_text(size = 34),axis.title.x = element_blank(),axis.text =  element_text(size = 30),
          legend.position = "none",strip.background = element_rect(fill = NA),strip.text = element_text(size = 34)))

#1300*800
#invSimp_Rain_Field_rar1000_dot_error

#####Plot Combined Field diversity####

plot_grid(field_rich,rain_rich,field_invSimp,rain_invSimp,nrow = 2,ncol = 2, align="v",axis = "l")
#Dot_error_Field_Rain_Diversity_raw

#####Diversity Analyses####




Mar_leaf_rain.fung_decon_rar_Mar_field_div_plant=subset(Mar_leaf_rain.fung_decon_rar_Mar_field_div,plant_type!="Rain")
Mar_leaf_rain.fung_decon_rar_Mar_field_div_plant$plant_date=with(Mar_leaf_rain.fung_decon_rar_Mar_field_div_plant,interaction(plant_type,collect_date))

#Richness
Rich_field_r1000_mod=lmer(Observed~plant_date+(1|gh_block),data = Mar_leaf_rain.fung_decon_rar_Mar_field_div_plant)
hist(resid(Rich_field_r1000_mod))
qqPlot(resid(Rich_field_r1000_mod))
shapiro.test(resid(Rich_field_r1000_mod))
#p-value = 0.2455
plot(Rich_field_r1000_mod)

anova(Rich_field_r1000_mod)
#plant_date  10712  2677.9     4 8.3702  4.9571 0.02445 *

#####POSTHOC TEST: Figure S11a####
emmeans(Rich_field_r1000_mod,pairwise~plant_date, adjust="fdr")

#invSimpson
invSimp_field_r1000_mod=lmer(sqrt(InvSimpson)~plant_date+(1|gh_block),data = Mar_leaf_rain.fung_decon_rar_Mar_field_div_plant)
hist(resid(invSimp_field_r1000_mod))
qqPlot(resid(invSimp_field_r1000_mod))
shapiro.test(resid(invSimp_field_r1000_mod))
#p-value = 0.2078
plot(invSimp_field_r1000_mod)

anova(invSimp_field_r1000_mod)
#plant_date   4.55  1.1375     4    24  0.6162 0.6552

#####POSTHOC TEST: Figure S11c####
emmeans(invSimp_field_r1000_mod,pairwise~plant_date, adjust="fdr")

#No Seeds
Mar_leaf_rain.fung_decon_rar_Mar_field_div_plant_NS=subset(Mar_leaf_rain.fung_decon_rar_Mar_field_div_plant,plant_type!="Seed")


#Richness
Rich_field_NS_r1000_mod=lmer((Observed)^2~plant_type*collect_date+(1|gh_block),data = Mar_leaf_rain.fung_decon_rar_Mar_field_div_plant_NS)
hist(resid(Rich_field_NS_r1000_mod))
qqPlot(resid(Rich_field_NS_r1000_mod))
shapiro.test(resid(Rich_field_NS_r1000_mod))
#p-value = 0.5266
plot(Rich_field_NS_r1000_mod)


####Table S10 Richness####
anova(Rich_field_NS_r1000_mod)
#collect_date            121295821 121295821     1 18.339  8.9777 0.007633 **
#plant_type:collect_date  58806443  58806443     1 18.532  4.3525 0.051027 . 

emmeans(Rich_field_NS_r1000_mod,pairwise~collect_date|plant_type)

emmeans(Rich_field_NS_r1000_mod,pairwise~collect_date*plant_type, adjust="fdr")

#invSimpson
invSimp_field_NS_r1000_mod=lmer(sqrt(InvSimpson)~plant_type*collect_date+(1|gh_block),data = Mar_leaf_rain.fung_decon_rar_Mar_field_div_plant_NS)
hist(resid(invSimp_field_NS_r1000_mod))
qqPlot(resid(invSimp_field_NS_r1000_mod))
shapiro.test(resid(invSimp_field_NS_r1000_mod))
#p-value = 0.231
plot(invSimp_field_NS_r1000_mod)

####Table S10 Inverse Simpson####
anova(invSimp_field_NS_r1000_mod)
#nada



#####Jaccard Field Analyses####


Mar_leaf_rain.fung_decon_rar_Mar_field_J.ord=ordinate(Mar_leaf_rain.fung_decon_rar_Mar_field,method = "NMDS",distance =  "jaccard",binary = TRUE)
#*** No convergence -- monoMDS stopping criteria:
#20: stress ratio > sratmax
#0.1536061     

#Creation of files for PERMANOVA analyses in Primer v6

Mar_leaf_rain.fung_decon_rar_Mar_field_J_dis=phyloseq::distance(Mar_leaf_rain.fung_decon_rar_Mar_field,method = "jaccard",binary = TRUE)
Mar_leaf_rain.fung_decon_rar_Mar_field_map=sample_data(Mar_leaf_rain.fung_decon_rar_Mar_field)
write.csv(as.matrix(Mar_leaf_rain.fung_decon_rar_Mar_field_J_dis),here::here("R_file","Jaccard_Mar_leaf_rain.fung_decon_rar_Mar_field_J_dis.csv"))



#####Field Jaccard Beta-dipsersion####

Mar_leaf_rain.fung_decon_rar_Mar_field_J_dis_betamod=betadisper(Mar_leaf_rain.fung_decon_rar_Mar_field_J_dis,
                                                                with(Mar_leaf_rain.fung_decon_rar_Mar_field_map,interaction(collect_date,plant_type)))

Mar_leaf_rain.fung_decon_rar_Mar_field_J_dis_betadisp=data.frame(Mar_leaf_rain.fung_decon_rar_Mar_field_J_dis_betamod$distances)

colnames(Mar_leaf_rain.fung_decon_rar_Mar_field_J_dis_betadisp)="betadisp"

Mar_leaf_rain.fung_decon_rar_Mar_field_J_dis_betadisp=merge(Mar_leaf_rain.fung_decon_rar_Mar_field_J_dis_betadisp,Mar_leaf_rain.fung_decon_rar_Mar_field_map,
                                                            by="row.names")


#####Graphing Jaccard Beta-dipsersion####

Mar_leaf_rain.fung_decon_rar_Mar_field_J_dis_betadisp_sum=Mar_leaf_rain.fung_decon_rar_Mar_field_J_dis_betadisp%>%group_by(collect_date,plant_type)%>%
  summarise_at(vars("betadisp"),list(~mean(.),~n(),se=~sd(.)/sqrt(n())))

with(Mar_leaf_rain.fung_decon_rar_Mar_field_J_dis_betadisp_sum,interaction(collect_date,plant_type))

plant_date_order=c("7/10/2018.Seed","7/16/2018.Adult","9/5/2018.Adult","7/16/2018.Seedling","9/5/2018.Seedling")

#Bray betadisp
#plants
(J_field_betadisp=ggplot(subset(Mar_leaf_rain.fung_decon_rar_Mar_field_J_dis_betadisp,plant_type!="Rain"))+
    geom_point(aes(x=factor(interaction(collect_date,plant_type),levels = plant_date_order,
                            labels = c("Seed", "Start","End","Start","End")),y=betadisp, shape=plant_type,
                   fill=factor(interaction(collect_date,plant_type),
                               levels = plant_date_order),
                   color=factor(interaction(collect_date,plant_type),
                                levels = plant_date_order)),size=2, alpha = 0.5)+
    geom_errorbar(data = subset(Mar_leaf_rain.fung_decon_rar_Mar_field_J_dis_betadisp_sum,plant_type!="Rain"),
                  aes(x=factor(interaction(collect_date,plant_type),levels = plant_date_order,
                               labels = c("Seed", "Start","End","Start","End")),
                      y=mean, ymax=mean+se,
                      ymin=mean-se), width=0.3)+
    geom_point(data = subset(Mar_leaf_rain.fung_decon_rar_Mar_field_J_dis_betadisp_sum,plant_type!="Rain"),
               aes(x=factor(interaction(collect_date,plant_type),
                            levels = plant_date_order,
                            labels = c("Seed", "Start","End","Start","End")),
                   y=mean, shape=plant_type,
                   fill=factor(interaction(collect_date,plant_type),
                               levels = plant_date_order),
                   color=factor(interaction(collect_date,plant_type),
                                levels = plant_date_order)),size=10)+scale_x_discrete(labels=)+
    scale_shape_manual(values = c(22,13,24),name=NULL)+ylab("Beta-dispersion\n(Jaccard)")+
    facet_grid(cols = vars(factor(plant_type,levels = c("Seed","Adult","Seedling"))),scales = "free_x",space = "free_x")+
    scale_color_manual(values = c("#E31A1C","black", "black","black","black"),name=NULL)+
    scale_fill_manual(values = c("#E31A1C","#A6CEE3","#1F78B4","#B2DF8A","#33A02C"),name=NULL)+
    theme_cowplot()+theme(axis.title.y = element_text(size = 34),axis.title.x = element_blank(),axis.text =  element_text(size = 30),
                          legend.position = "none",strip.background = element_rect(fill = NA),strip.text = element_text(size = 34)))

#1000*800




#Rain
(J_rain_betadisp=ggplot(subset(Mar_leaf_rain.fung_decon_rar_Mar_field_J_dis_betadisp,plant_type=="Rain"))+
    geom_point(aes(x=as.Date(collect_date,format="%m/%d/%Y"),y=betadisp),size=2, alpha = 0.5, shape=21, fill="black")+
    geom_errorbar(data = subset(Mar_leaf_rain.fung_decon_rar_Mar_field_J_dis_betadisp_sum,plant_type=="Rain"),aes(x=as.Date(collect_date,format="%m/%d/%Y"),
                                                                                                                  y=mean, ymax=mean+se,
                                                                                                                  ymin=mean-se), width=1.5)+
    geom_point(data = subset(Mar_leaf_rain.fung_decon_rar_Mar_field_J_dis_betadisp_sum,plant_type=="Rain"),aes(x=as.Date(collect_date,format="%m/%d/%Y"),
                                                                                                               y=mean),size=10,shape=21, fill="black")+
    theme_cowplot()+scale_x_date(date_labels = "%b %d",date_breaks = "11 days")+ylab("Beta-dispersion\n(Jaccard)")+
    facet_grid(cols = vars(factor(plant_type,labels = "")))+
    theme(axis.title.y = element_text(size = 34),axis.title.x = element_blank(),axis.text =  element_text(size = 30),
          legend.position = "none",strip.background = element_rect(fill = NA),strip.text = element_text(size = 34)))

#1300*800
#




#####Field Jaccard Betadisp Stats####




Mar_leaf_rain.fung_decon_rar_Mar_field_J_dis_betadisp_plant=subset(Mar_leaf_rain.fung_decon_rar_Mar_field_J_dis_betadisp,plant_type!="Rain")
Mar_leaf_rain.fung_decon_rar_Mar_field_J_dis_betadisp_plant$plant_date=with(Mar_leaf_rain.fung_decon_rar_Mar_field_J_dis_betadisp_plant,
                                                                            interaction(plant_type,collect_date))

#all plant material
Betdisp_field_J_r1000_mod=lmer((betadisp)^(-1)~plant_date+(1|gh_block),data = Mar_leaf_rain.fung_decon_rar_Mar_field_J_dis_betadisp_plant)
hist(resid(Betdisp_field_J_r1000_mod))
qqPlot(resid(Betdisp_field_J_r1000_mod))
shapiro.test(resid(Betdisp_field_J_r1000_mod))
#p-value = 0.04292
plot(Betdisp_field_J_r1000_mod)

anova(Betdisp_field_J_r1000_mod)
#plant_date 0.95369 0.23842     4 8.9793  1.6789 0.2381

emmeans(Betdisp_field_r1000_mod,pairwise~plant_date, adjust="fdr")

#Seedling and Adult only 
Mar_leaf_rain.fung_decon_rar_Mar_field_J_dis_betadisp_LFE=subset(Mar_leaf_rain.fung_decon_rar_Mar_field_J_dis_betadisp_plant,plant_type!="Seed")

Betdisp_field_J_LFE_r1000_mod=lmer((betadisp)^(-1)~plant_type*collect_date+(1|gh_block),data = Mar_leaf_rain.fung_decon_rar_Mar_field_J_dis_betadisp_LFE)
hist(resid(Betdisp_field_J_LFE_r1000_mod))
qqPlot(resid(Betdisp_field_J_LFE_r1000_mod))
shapiro.test(resid(Betdisp_field_J_LFE_r1000_mod))
#p-value = 0.07337
plot(Betdisp_field_J_LFE_r1000_mod)

####Table S10 Beta-dispersion (Jaccard)####
anova(Betdisp_field_J_LFE_r1000_mod)
#plant_type              0.08345 0.08345     1 18.656  0.5836 0.4545
#collect_date            0.36076 0.36076     1 18.342  2.5228 0.1293
#plant_type:collect_date 0.40222 0.40222     1 18.656  2.8127 0.1102 

#emmeans(Betdisp_field_J_LFE_r1000_mod,pairwise~plant_type*collect_date, adjust="fdr")


#####PLOT: Figure S13####
#####Plot Betadisp Field Bray and Jaccard####

#Bray betadisp
#plants
(B_field_betadisp2=ggplot(subset(Mar_leaf_rain.fung_decon_rar_Mar_field_dis_betadisp,plant_type!="Rain"))+
   geom_point(aes(x=factor(interaction(collect_date,plant_type),levels = plant_date_order,
                           labels = c("Seed", "Start","End","Start","End")),y=betadisp, shape=plant_type,
                  fill=factor(interaction(collect_date,plant_type),
                              levels = plant_date_order),
                  color=factor(interaction(collect_date,plant_type),
                               levels = plant_date_order)),size=2, alpha = 0.5)+
   geom_errorbar(data = subset(Mar_leaf_rain.fung_decon_rar_Mar_field_dis_betadisp_sum,plant_type!="Rain"),
                 aes(x=factor(interaction(collect_date,plant_type),levels = plant_date_order,
                              labels = c("Seed", "Start","End","Start","End")),
                     y=mean, ymax=mean+se,
                     ymin=mean-se), width=0.3)+
   geom_point(data = subset(Mar_leaf_rain.fung_decon_rar_Mar_field_dis_betadisp_sum,plant_type!="Rain"),
              aes(x=factor(interaction(collect_date,plant_type),
                           levels = plant_date_order,
                           labels = c("Seed", "Start","End","Start","End")),
                  y=mean, shape=plant_type,
                  fill=factor(interaction(collect_date,plant_type),
                              levels = plant_date_order),
                  color=factor(interaction(collect_date,plant_type),
                               levels = plant_date_order)),size=10)+scale_x_discrete(labels=)+
   scale_shape_manual(values = c(22,13,24),name=NULL)+
   scale_y_continuous(name = "Beta-dispersion\n(Bray-Curtis)",limits = c(0.12,0.75))+
   facet_grid(cols = vars(factor(plant_type,levels = c("Seed","Adult","Seedling"))),scales = "free_x",space = "free_x")+
   scale_color_manual(values = c("#E31A1C","black", "black","black","black"),name=NULL)+
   scale_fill_manual(values = c("#E31A1C","#A6CEE3","#1F78B4","#B2DF8A","#33A02C"),name=NULL)+
   theme_cowplot()+theme(axis.title.y = element_text(size = 34),axis.title.x = element_blank(),axis.text =  element_text(size = 30),
                         legend.position = "none",strip.background = element_rect(fill = NA),strip.text = element_text(size = 34)))

#1000*800

(J_field_betadisp2=ggplot(subset(Mar_leaf_rain.fung_decon_rar_Mar_field_J_dis_betadisp,plant_type!="Rain"))+
    geom_point(aes(x=factor(interaction(collect_date,plant_type),levels = plant_date_order,
                            labels = c("Seed", "Start","End","Start","End")),y=betadisp, shape=plant_type,
                   fill=factor(interaction(collect_date,plant_type),
                               levels = plant_date_order),
                   color=factor(interaction(collect_date,plant_type),
                                levels = plant_date_order)),size=2, alpha = 0.5)+
    geom_errorbar(data = subset(Mar_leaf_rain.fung_decon_rar_Mar_field_J_dis_betadisp_sum,plant_type!="Rain"),
                  aes(x=factor(interaction(collect_date,plant_type),levels = plant_date_order,
                               labels = c("Seed", "Start","End","Start","End")),
                      y=mean, ymax=mean+se,
                      ymin=mean-se), width=0.3)+
    geom_point(data = subset(Mar_leaf_rain.fung_decon_rar_Mar_field_J_dis_betadisp_sum,plant_type!="Rain"),
               aes(x=factor(interaction(collect_date,plant_type),
                            levels = plant_date_order,
                            labels = c("Seed", "Start","End","Start","End")),
                   y=mean, shape=plant_type,
                   fill=factor(interaction(collect_date,plant_type),
                               levels = plant_date_order),
                   color=factor(interaction(collect_date,plant_type),
                                levels = plant_date_order)),size=10)+scale_x_discrete(labels=)+
    scale_shape_manual(values = c(22,13,24),name=NULL)+
    scale_y_continuous(name = "Beta-dispersion\n(Jaccard)",limits = c(0.25,0.73))+
    facet_grid(cols = vars(factor(plant_type,levels = c("Seed","Adult","Seedling"))),scales = "free_x",space = "free_x")+
    scale_color_manual(values = c("#E31A1C","black", "black","black","black"),name=NULL)+
    scale_fill_manual(values = c("#E31A1C","#A6CEE3","#1F78B4","#B2DF8A","#33A02C"),name=NULL)+
    theme_cowplot()+theme(axis.title.y = element_text(size = 34),axis.title.x = element_blank(),axis.text =  element_text(size = 30),
                          legend.position = "none",strip.background = element_rect(fill = NA),strip.text = element_text(size = 34)))


#Rain
(B_rain_betadisp2=ggplot(subset(Mar_leaf_rain.fung_decon_rar_Mar_field_dis_betadisp,plant_type=="Rain"))+
    geom_point(aes(x=as.Date(collect_date,format="%m/%d/%Y"),y=betadisp),size=2, alpha = 0.5, shape=21, fill="black")+
    geom_errorbar(data = subset(Mar_leaf_rain.fung_decon_rar_Mar_field_dis_betadisp_sum,plant_type=="Rain"),aes(x=as.Date(collect_date,format="%m/%d/%Y"),
                                                                                                                y=mean, ymax=mean+se,
                                                                                                                ymin=mean-se), width=1.5)+
    geom_point(data = subset(Mar_leaf_rain.fung_decon_rar_Mar_field_dis_betadisp_sum,plant_type=="Rain"),aes(x=as.Date(collect_date,format="%m/%d/%Y"),
                                                                                                             y=mean),size=10,shape=21, fill="black")+
    theme_cowplot()+scale_x_date(date_labels = "%b %d",date_breaks = "11 days")+
    scale_y_continuous(name = "Beta-dispersion\n(Bray-Curtis)",limits = c(0.12,0.75))+
    facet_grid(cols = vars(factor(plant_type,labels = "")))+
    theme(axis.title= element_blank(),axis.text =  element_text(size = 30),axis.text.y = element_blank(),
          legend.position = "none",strip.background = element_rect(fill = NA),strip.text = element_text(size = 34)))





#Rain
(J_rain_betadisp2=ggplot(subset(Mar_leaf_rain.fung_decon_rar_Mar_field_J_dis_betadisp,plant_type=="Rain"))+
    geom_point(aes(x=as.Date(collect_date,format="%m/%d/%Y"),y=betadisp),size=2, alpha = 0.5, shape=21, fill="black")+
    geom_errorbar(data = subset(Mar_leaf_rain.fung_decon_rar_Mar_field_J_dis_betadisp_sum,plant_type=="Rain"),aes(x=as.Date(collect_date,format="%m/%d/%Y"),
                                                                                                                  y=mean, ymax=mean+se,
                                                                                                                  ymin=mean-se), width=1.5)+
    geom_point(data = subset(Mar_leaf_rain.fung_decon_rar_Mar_field_J_dis_betadisp_sum,plant_type=="Rain"),aes(x=as.Date(collect_date,format="%m/%d/%Y"),
                                                                                                               y=mean),size=10,shape=21, fill="black")+
    theme_cowplot()+scale_x_date(date_labels = "%b %d",date_breaks = "11 days")+
    scale_y_continuous(name = "Beta-dispersion\n(Jaccard)",limits = c(0.25,0.73))+
    facet_grid(cols = vars(factor(plant_type,labels = "")))+
    theme(axis.title = element_blank(),axis.text =  element_text(size = 30),axis.text.y = element_blank(),
          legend.position = "none",strip.background = element_rect(fill = NA),strip.text = element_text(size = 34)))

plot_grid(B_field_betadisp2,B_rain_betadisp2,J_field_betadisp2,J_rain_betadisp2,nrow = 2,ncol = 2, align="h",axis = "l",
          rel_widths = c(1.15,1))
#Field_betadsp_Bray_Jacc_raw


####Graph Field experiment Jaccard NMDS####

#Hull for the Jaccard NMDS
Mar_leaf_rain.fung_decon_rar_Mar_field_map=sample_data(Mar_leaf_rain.fung_decon_rar_Mar_field)
Mar_leaf_rain.fung_decon_rar_Mar_field_J_coordinates=merge(Mar_leaf_rain.fung_decon_rar_Mar_field_J.ord$points,Mar_leaf_rain.fung_decon_rar_Mar_field_map, by="row.names")
Mar_leaf_rain.fung_decon_rar_Mar_field_J_coordinates$plant_data_grp=ifelse(Mar_leaf_rain.fung_decon_rar_Mar_field_J_coordinates$plant_type=="Rain",
                                                                           "Rain", ifelse(Mar_leaf_rain.fung_decon_rar_Mar_field_J_coordinates$plant_type=="Seedling",
                                                                                          ifelse(Mar_leaf_rain.fung_decon_rar_Mar_field_J_coordinates$collect_date=="9/5/2018","Sept_Seedling","Jul_Seedling"),
                                                                                          ifelse(Mar_leaf_rain.fung_decon_rar_Mar_field_J_coordinates$plant_type=="Adult",
                                                                                                 ifelse(Mar_leaf_rain.fung_decon_rar_Mar_field_J_coordinates$collect_date=="9/5/2018","Sept_Adult","Jul_Adult"),"Seed")))
grp1000_J.Adult_Jul <- Mar_leaf_rain.fung_decon_rar_Mar_field_J_coordinates[Mar_leaf_rain.fung_decon_rar_Mar_field_J_coordinates$plant_type == "Adult"&
                                                                              Mar_leaf_rain.fung_decon_rar_Mar_field_J_coordinates$collect_date == "7/16/2018", 
][chull(Mar_leaf_rain.fung_decon_rar_Mar_field_J_coordinates[
  Mar_leaf_rain.fung_decon_rar_Mar_field_J_coordinates$
    plant_type =="Adult"&Mar_leaf_rain.fung_decon_rar_Mar_field_J_coordinates$
    collect_date == "7/16/2018", c("MDS1", "MDS2")]), ] 
nrow(grp1000_J.Adult_Jul)
#3
grp1000_J.Adult_Sept <- Mar_leaf_rain.fung_decon_rar_Mar_field_J_coordinates[Mar_leaf_rain.fung_decon_rar_Mar_field_J_coordinates$plant_type == "Adult"&
                                                                               Mar_leaf_rain.fung_decon_rar_Mar_field_J_coordinates$collect_date == "9/5/2018", 
][chull(Mar_leaf_rain.fung_decon_rar_Mar_field_J_coordinates[
  Mar_leaf_rain.fung_decon_rar_Mar_field_J_coordinates$
    plant_type =="Adult"&Mar_leaf_rain.fung_decon_rar_Mar_field_J_coordinates$
    collect_date == "9/5/2018", c("MDS1", "MDS2")]), ] 

nrow(grp1000_J.Adult_Sept)
#4


grp1000_J.Seedling_Jul <- Mar_leaf_rain.fung_decon_rar_Mar_field_J_coordinates[Mar_leaf_rain.fung_decon_rar_Mar_field_J_coordinates$plant_type == "Seedling"&
                                                                                 Mar_leaf_rain.fung_decon_rar_Mar_field_J_coordinates$collect_date == "7/16/2018", 
][chull(Mar_leaf_rain.fung_decon_rar_Mar_field_J_coordinates[
  Mar_leaf_rain.fung_decon_rar_Mar_field_J_coordinates$
    plant_type =="Seedling"&Mar_leaf_rain.fung_decon_rar_Mar_field_J_coordinates$
    collect_date == "7/16/2018", c("MDS1", "MDS2")]), ] 
nrow(grp1000_J.Seedling_Jul)
#3
grp1000_J.Seedling_Sept <- Mar_leaf_rain.fung_decon_rar_Mar_field_J_coordinates[Mar_leaf_rain.fung_decon_rar_Mar_field_J_coordinates$plant_type == "Seedling"&
                                                                                  Mar_leaf_rain.fung_decon_rar_Mar_field_J_coordinates$collect_date == "9/5/2018", 
][chull(Mar_leaf_rain.fung_decon_rar_Mar_field_J_coordinates[
  Mar_leaf_rain.fung_decon_rar_Mar_field_J_coordinates$
    plant_type =="Seedling"&Mar_leaf_rain.fung_decon_rar_Mar_field_J_coordinates$
    collect_date == "9/5/2018", c("MDS1", "MDS2")]), ] 

nrow(grp1000_J.Seedling_Sept)
#4

grp1000_J.Rain <- Mar_leaf_rain.fung_decon_rar_Mar_field_J_coordinates[Mar_leaf_rain.fung_decon_rar_Mar_field_J_coordinates$plant_type == "Rain", 
][chull(Mar_leaf_rain.fung_decon_rar_Mar_field_J_coordinates[
  Mar_leaf_rain.fung_decon_rar_Mar_field_J_coordinates$
    plant_type =="Rain", c("MDS1", "MDS2")]), ] 

nrow(grp1000_J.Rain)
#8

grp1000_J.Seed <- Mar_leaf_rain.fung_decon_rar_Mar_field_J_coordinates[Mar_leaf_rain.fung_decon_rar_Mar_field_J_coordinates$plant_type == "Seed", 
][chull(Mar_leaf_rain.fung_decon_rar_Mar_field_J_coordinates[
  Mar_leaf_rain.fung_decon_rar_Mar_field_J_coordinates$
    plant_type =="Seed", c("MDS1", "MDS2")]), ] 

nrow(grp1000_J.Seed)
#4


hull1000_field_J=rbind(grp1000_J.Adult_Jul,grp1000_J.Adult_Sept,grp1000_J.Seedling_Jul,grp1000_J.Seedling_Sept,grp1000_J.Rain, grp1000_J.Seed)
nrow(hull1000_field_J)
#26
Mar_leaf_rain.fung_decon_rar_Mar_field_J_coor_rain=subset(Mar_leaf_rain.fung_decon_rar_Mar_field_J_coordinates, plant_type=="Rain")
nrow(Mar_leaf_rain.fung_decon_rar_Mar_field_J_coor_rain)
#56
Mar_leaf_rain.fung_decon_rar_Mar_field_J_coor_plant=subset(Mar_leaf_rain.fung_decon_rar_Mar_field_J_coordinates, plant_type!="Rain"&
                                                             plant_type!="Seed")
nrow(Mar_leaf_rain.fung_decon_rar_Mar_field_J_coor_plant)
#25
Mar_leaf_rain.fung_decon_rar_Mar_field_J_coor_seed=subset(Mar_leaf_rain.fung_decon_rar_Mar_field_J_coordinates,plant_type=="Seed")
nrow(Mar_leaf_rain.fung_decon_rar_Mar_field_J_coor_seed)
#4
library(RColorBrewer)
brewer.pal(n = 8, name = "Paired")
with(Mar_leaf_rain.fung_decon_rar_Mar_field_J_coordinates,interaction(collect_date,plant_type))


unique(Mar_leaf_rain.fung_decon_rar_Mar_field_J_coordinates$plant_data_grp)

plant_data_grp_order=c("Seed","Jul_Seedling","Sept_Seedling","Jul_Adult","Sept_Adult","Rain")
unique(Mar_leaf_rain.fung_decon_rar_Mar_field_J_coordinates$plant_type)


ggplot(Mar_leaf_rain.fung_decon_rar_Mar_field_J_coordinates,aes(x= MDS1, y=MDS2))+
  geom_polygon(data=hull1000_field_J,aes(x=MDS1,y=MDS2,fill=factor(plant_data_grp,levels = plant_data_grp_order),
                                         group=factor(plant_data_grp,levels = plant_data_grp_order)), alpha=0.5)+
  geom_point(size=4,aes(fill=factor(plant_data_grp,levels = plant_data_grp_order),
                        shape=factor(plant_type,levels = c("Adult","Seed","Seedling","Rain")), color=factor(plant_data_grp,levels = plant_data_grp_order)),stroke=2)+
  scale_shape_manual(values = c(22,13,24,21),name=NULL)+
  scale_color_manual(values = c("#E31A1C","black", "black","black","black","black"),name=NULL)+
  scale_fill_manual(values = c("#E31A1C","#B2DF8A","#33A02C","#A6CEE3","#1F78B4","black"),name=NULL)+
  theme_bw()+theme(axis.title = element_text(size = 30),axis.text =  element_text(size = 30),legend.position = "none")


#900x800
#####PLOT: Figure S9####
#####Plot field NMDS####

(NMDS_Field_Bray=ggplot(Mar_leaf_rain.fung_decon_rar_Mar_field_coordinates,aes(x= MDS1, y=MDS2))+
   geom_polygon(data=hull1000_Field,aes(x=MDS1,y=MDS2,fill=factor(plant_data_grp,levels = plant_data_grp_order),
                                        group=factor(plant_data_grp,levels = plant_data_grp_order)), alpha=0.5)+
   geom_point(size=4,aes(fill=factor(plant_data_grp,levels = plant_data_grp_order),
                         shape=factor(plant_type,levels = c("Adult","Seed","Seedling","Rain")), 
                         color=factor(plant_data_grp,levels = plant_data_grp_order)),stroke=2)+
   scale_shape_manual(values = c(22,13,24,21),name=NULL)+xlab("NMDS1")+scale_y_continuous(name = "NMDS2",limits = c(-1.3,0.7))+
   scale_color_manual(values = c("#E31A1C","black", "black","black","black","black"),name=NULL)+
   scale_fill_manual(values = c("#E31A1C","#B2DF8A","#33A02C","#A6CEE3","#1F78B4","black"),name=NULL)+
   theme_cowplot()+theme(axis.title.y = element_text(size = 36),axis.text.y =  element_text(size = 30),
                         axis.title.x = element_blank(),axis.text.x =  element_blank(),legend.position = "none",
                         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_rect(colour = "black")))






(NMDS_Field_Jacc=ggplot(Mar_leaf_rain.fung_decon_rar_Mar_field_J_coordinates,aes(x= MDS1, y=MDS2))+
    geom_polygon(data=hull1000_field_J,aes(x=MDS1,y=MDS2,fill=factor(plant_data_grp,levels = plant_data_grp_order),
                                           group=factor(plant_data_grp,levels = plant_data_grp_order)), alpha=0.5)+
    geom_point(size=4,aes(fill=factor(plant_data_grp,levels = plant_data_grp_order),
                          shape=factor(plant_type,levels = c("Adult","Seed","Seedling","Rain")), color=factor(plant_data_grp,levels = plant_data_grp_order)),stroke=2)+
    scale_shape_manual(values = c(22,13,24,21),name=NULL)+xlab("NMDS1")+scale_y_continuous(name = "NMDS2",limits = c(-1.3,0.7))+
    scale_color_manual(values = c("#E31A1C","black", "black","black","black","black"),name=NULL)+
    scale_fill_manual(values = c("#E31A1C","#B2DF8A","#33A02C","#A6CEE3","#1F78B4","black"),name=NULL)+
    theme_cowplot()+theme(axis.title = element_text(size = 30),axis.text =  element_text(size = 30),legend.position = "none",
                          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_rect(colour = "black")))



plot_grid(NMDS_Field_Bray,NMDS_Field_Jacc,nrow = 2, align="hv")


#####Field Published leaf endophyte sequence analyses####
#Run using Usearch 

#cd HardDrive/FungiRainLeaf2019/R_file/
#The program needs for the codons to be uppercase
#awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' rep_set.FunRainLeaf2019_v4.2.2020_ZOTU_rar1000.fna > CAP_rep_set.FunRainLeaf2019_v4.2.2020_ZOTU_rar1000.fna
#~/HardDrive/Sciencey_Program/usearch11 -usearch_global CAP_rep_set.FunRainLeaf2019_v4.2.2020_ZOTU_rar1000.fna -db switchgrass_compiled_ITS_database_v2.fasta -id 0.97 -strand both -maxaccepts 0 -maxhits 10 -matched MATCHED_pub_seq_rep_set.FunRainLeaf2019_ZOTU_v4.2.2020_rar1000.fa -notmatched NOT_matched_pub_seq_rep_set.FunRainLeaf2019_ZOTU_v4.2.2020_rar1000.fa -userout TBL_pub_seq_rep_set.FunRainLeaf2019_ZOTU_v4.2.2020_rar1000.txt -userfields query+target+id+mid+bits+evalue+ql+ts+qlor+qhir+tlor+thir
#00:02 22Mb    100.0% Searching, 5.2% matched

publ_matched_rar_OTUs=read.delim(here::here("R_file","TBL_pub_seq_rep_set.FunRainLeaf2019_ZOTU_v4.2.2020_rar1000.txt"), header = F)
head(publ_matched_rar_OTUs)

colnames(publ_matched_rar_OTUs)=c("query","target","id","mid","bits","evalue","ql","ts","qlor","qhir","tlor","thir")

nrow(publ_matched_rar_OTUs)
#468
#Subset to have only the top hiting taxa 
publ_matched_rar_OTUs_top_hits=publ_matched_rar_OTUs %>% group_by(query) %>% slice(which(id == max(id)))

nrow(publ_matched_rar_OTUs_top_hits)
#258






load(here::here("R_file","Mar_leaf_rain_ZOTU_full_v4.2.2020.fung_decon_rar_pruned_phyloseq_obj.RData"))

Mar_leaf_rain.fung_decon_rar_Mar=subset_samples(Mar_leaf_rain.fung_decon_rar, sub_proj=="Marshall")
Mar_leaf_rain.fung_decon_rar_Mar_field=subset_samples(Mar_leaf_rain.fung_decon_rar_Mar, plant_type!="Petri")
nsamples(Mar_leaf_rain.fung_decon_rar_Mar_field)
#85

unique(sample_data(Mar_leaf_rain.fung_decon_rar_Mar_field)$sub_proj)

Mar_leaf_rain.fung_decon_rar_Mar_field=prune_taxa(taxa_sums(Mar_leaf_rain.fung_decon_rar_Mar_field) > 0, Mar_leaf_rain.fung_decon_rar_Mar_field)
ntaxa(Mar_leaf_rain.fung_decon_rar_Mar_field)
#2243
sum(taxa_sums(Mar_leaf_rain.fung_decon_rar_Mar_field))
#85000
min(sample_sums(Mar_leaf_rain.fung_decon_rar_Mar_field))
#1000
max(sample_sums(Mar_leaf_rain.fung_decon_rar_Mar_field))
#1000
head(publ_matched_rar_OTUs_top_hits)

Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo=merge(otu_table(Mar_leaf_rain.fung_decon_rar_Mar_field),
                                                      publ_matched_rar_OTUs_top_hits,
                                                      by.x="row.names", by.y="query",all.x=T)
head(Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo)
colnames(Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo)
Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo$OTU_sum=rowSums(Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo[,2:86])
Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_sort=Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo[order(-Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo$OTU_sum),]

unique(Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo$target)
#40

#Replace NAs with not matched 

Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo[is.na(Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo)] <- "Not_matched"


#Let's look at the grouped Taxa

colnames(Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo)
summary(Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo)
nrow(Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo)
#2267
#Group by interaction and indicator type

#Load in the names table
publ_matchednames=read.csv(here::here("R_file","Matched_pub_seq_names_TBL.csv"), header = T)
colnames(publ_matchednames)
Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_names=merge(Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo,publ_matchednames,by = "target")
nrow(Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_names)
#2267
colnames(Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_names)

#First, let' look at the groups

Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_names_derup=
  Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_names[!duplicated(Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_names[,c("Row.names","Group")]),]

nrow(Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_names_derup)
#2258

Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_names_derup_sum=
  Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_names_derup[,c(3:87,102)]%>%group_by(Group)%>%summarise_all(~sum(.))




Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_names_derup_sum_M=melt(Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_names_derup_sum)

#Merge with the metadata
Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_names_derup_sum_M_trt=merge(Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_names_derup_sum_M,sample_data(Mar_leaf_rain.fung_decon_rar_Mar_field),by.x = "variable",
                                                                            by.y = "row.names")


Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_names_derup_sum_M_trt2=merge(Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_names_derup_sum_M_trt,
                                                                             sample_sums(Mar_leaf_rain.fung_decon_rar_Mar_field),by.x = "variable",
                                                                             by.y = "row.names")

Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_names_derup_sum_M_trt2$prop_abun=Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_names_derup_sum_M_trt2$value/
  Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_names_derup_sum_M_trt2$y

unique(Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_names_derup_sum_M_trt2$Group)


#Second, let' look at the Interactions

Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_names_derup_int=
  Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_names[!duplicated(Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_names[,c("Row.names","Interaction")]),]
unique(Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_names_derup_int$Interaction)
nrow(Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_names_derup_int)
#2258

Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_names_derup_int_sum=
  Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_names_derup_int[,c(3:87,103)]%>%group_by(Interaction)%>%summarise_all(~sum(.))


unique(Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_names_derup_int_sum$Interaction)

Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_names_derup_int_sum_M=melt(Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_names_derup_int_sum)

#Merge with the metadata
Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_names_derup_int_sum_trt=merge(Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_names_derup_int_sum_M,sample_data(Mar_leaf_rain.fung_decon_rar_Mar_field),by.x = "variable",
                                                                              by.y = "row.names")


Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_names_derup_int_sum_trt2=merge(Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_names_derup_int_sum_trt,
                                                                               sample_sums(Mar_leaf_rain.fung_decon_rar_Mar_field),by.x = "variable",
                                                                               by.y = "row.names")

Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_names_derup_int_sum_trt2$prop_abun=Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_names_derup_int_sum_trt2$value/
  Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_names_derup_int_sum_trt2$y


#Let's combine the interaction groups with the indicators
colnames(Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_names_derup_int_sum_trt2)[colnames(Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_names_derup_int_sum_trt2)=="Interaction"]="grp_inter"
colnames(Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_names_derup_sum_M_trt2)[colnames(Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_names_derup_sum_M_trt2)=="Group"]="grp_inter"

Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_comb=rbind(subset(Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_names_derup_int_sum_trt2,
                                                                  grp_inter!="Unknown"),
                                                           subset(Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_names_derup_sum_M_trt2,
                                                                  grp_inter=="Cave-in-Rock Indicator"|grp_inter=="Prairie Indicator"))

Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_comb_not_match=subset(Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_comb,grp_inter=="Not_matched")
Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_comb_not_match$grp_inter=rep("Endophyte")
Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_comb_not_match$value=1000-Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_comb_not_match$value
Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_comb_not_match$prop_abun=1-Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_comb_not_match$prop_abun
Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_comb=rbind(Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_comb_not_match,
                                                           Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_comb)

Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_comb$plant_data_grp=ifelse(Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_comb$plant_type=="Rain",
                                                                           "Rain", ifelse(Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_comb$plant_type=="Seedling",
                                                                                          ifelse(Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_comb$collect_date=="9/5/2018","Sept_Seedling","Jul_Seedling"),
                                                                                          ifelse(Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_comb$plant_type=="Adult",
                                                                                                 ifelse(Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_comb$collect_date=="9/5/2018","Sept_Adult","Jul_Adult"),"Seed")))




Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_comb_sum=Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_comb%>%group_by(grp_inter,plant_data_grp)%>%
  summarise_at("prop_abun",list(~mean(.),se=~sd(.)/sqrt(n())))
unique(Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_comb_sum$grp_inter)

#####PLOT: Figure 5####
#Plants and Rain


grp_inter_order=c("Endophyte","Pathogen","Context Mutualist","Mutualist","Prairie Indicator","Cave-in-Rock Indicator")
data_plant=c("Seed","Jul_Seedling","Sept_Seedling","Jul_Adult","Sept_Adult","Rain")

ggplot(subset(Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_comb,
              grp_inter!="Not_matched")) + 
  geom_point(data=subset(Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_comb,
                         grp_inter!="Not_matched"),
             aes(x=factor(plant_data_grp,levels = data_plant),fill=factor(plant_data_grp,levels = data_plant),
                 color=factor(plant_data_grp,levels = data_plant),y=prop_abun,shape=factor(plant_data_grp,levels = data_plant),
                 group=variable),size=3,alpha=0.15,position = position_dodge2(0.8))+
  geom_errorbar(data=subset(Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_comb_sum,
                            grp_inter!="Not_matched"),
                aes(x=factor(plant_data_grp,levels = data_plant),y=mean,ymax=mean+se,ymin=mean-se),width=0.5)+
  geom_point(data=subset(Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_comb_sum,
                         grp_inter!="Not_matched"),
             size=10,aes(x=factor(plant_data_grp,levels = data_plant),fill=factor(plant_data_grp,levels = data_plant),
                         color=factor(plant_data_grp,levels = data_plant),y=mean,shape=factor(plant_data_grp,levels = data_plant)))+
  facet_wrap(vars(factor(grp_inter,levels = grp_inter_order)),scales="free_y")+
  scale_shape_manual(values = c(13,24,24,22,22,21))+
  scale_x_discrete(labels=c("Seed", "Start","End","Start","End","Rain"),name=NULL)+theme_cowplot()+
  ylab("Proportion of read\nper sample")+
  scale_color_manual(values = c("#E31A1C","black","black","black","black","black"),name=NULL)+
  scale_fill_manual(values = c("#E31A1C","#B2DF8A","#33A02C","#A6CEE3","#1F78B4","black"),name=NULL)+
  theme(axis.title = element_text(size = 30),axis.text = element_text(size = 24),legend.position = "none",strip.text = element_text(size = 28),
        strip.background = element_rect(fill = NA,colour = "black"))
#2000*900
#Field_pub_seq_grp_prop_dot

#Mean Endophytes in Rain 
subset(Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_comb_sum,plant_data_grp=="Rain")


#####PLOT: Figure S14####
#Rain

Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_comb_R=subset(Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_comb,plant_data_grp=="Rain")
Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_comb_R_sum=Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_comb_R%>%
  group_by(grp_inter,collect_date)%>%summarise_at(vars(prop_abun),list(~mean(.),se=~sd(.)/sqrt(n())))

ggplot(subset(Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_comb_R,grp_inter!="Not_matched")) +
  geom_point(data = subset(Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_comb_R,grp_inter!="Not_matched"),aes(x=as.Date(collect_date,format="%m/%d/%Y"),y=prop_abun),
             size=2,shape=21, fill="black",alpha=0.15)+
  geom_point(data = subset(Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_comb_R_sum,grp_inter!="Not_matched"),aes(x=as.Date(collect_date,format="%m/%d/%Y"),y=mean),
             size=7,shape=21, fill="black")+
  geom_errorbar(data = subset(Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_comb_R_sum,grp_inter!="Not_matched"),
                aes(x=as.Date(collect_date,format="%m/%d/%Y"),y=mean,ymax=mean+se,ymin=mean-se),width=1)+
  theme_cowplot()+scale_x_date(date_labels = "%b %d",date_breaks = "2 week",name = NULL)+
  ylab("Proportion of read\nper sample")+facet_wrap(vars(factor(grp_inter,levels = grp_inter_order)),scales="free_y")+
  theme(axis.title = element_text(size = 30),axis.text = element_text(size = 24),legend.position = "none",strip.text = element_text(size = 28),
        strip.background = element_rect(fill = NA,colour = "black"))
#2000*900
#Rain_pub_seq_grp_prop_dot


#####Field Statistical analyses of Pub seqs####


#Endophytes


Endop_field_r1000_mod=lmer(sqrt(value)~plant_data_grp+(1|gh_block),data = subset(Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_comb,grp_inter=="Endophyte"))
hist(resid(Endop_field_r1000_mod))
qqPlot(resid(Endop_field_r1000_mod))
shapiro.test(resid(Endop_field_r1000_mod))
#p-value = 0.536
plot(Endop_field_r1000_mod)

anova(Endop_field_r1000_mod)
#plant_data_grp 754.68  150.94     5    79  9.3723 4.844e-07 ***

emmeans(Endop_field_r1000_mod,pairwise~plant_data_grp, adjust="fdr")
#Jul_Adult - Sept_Seedling      -8.885 2.56 77.6 -3.472  0.0032
#Jul_Seedling - Sept_Seedling  -10.286 2.56 77.6 -4.022  0.0007 
#Rain - Sept_Seedling           -7.090 1.17 76.7 -6.059  <.0001
#Sept_Adult - Sept_Seedling     -9.898 2.26 76.1 -4.382  0.0003 


#Context Mutualist


C_Mut_field_r1000_mod=lmer(sqrt(value)~plant_data_grp+(1|gh_block),data = subset(Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_comb,grp_inter=="Context Mutualist"))
hist(resid(C_Mut_field_r1000_mod))
qqPlot(resid(C_Mut_field_r1000_mod))
shapiro.test(resid(C_Mut_field_r1000_mod))
#p-value = 0.1416
plot(C_Mut_field_r1000_mod)

anova(C_Mut_field_r1000_mod)
#plant_data_grp 105.58  21.117     5 58.704  8.3959 4.873e-06 ***


#####POSTHOC TEST: Figure 5c####
emmeans(C_Mut_field_r1000_mod,pairwise~plant_data_grp, adjust="fdr")
#Jul_Adult - Rain              -2.4576 0.947 77.2 -2.595  0.0242 
#Jul_Adult - Sept_Adult        -1.8069 1.217 76.7 -1.485  0.2124 
#Jul_Adult - Sept_Seedling     -4.2488 1.011 77.2 -4.204  0.0004 
#Jul_Seedling - Sept_Seedling  -4.1787 1.011 77.2 -4.134  0.0004 
#Rain - Seed                    2.4072 0.864 26.3  2.785  0.0242
#Rain - Sept_Seedling          -1.7912 0.462 76.5 -3.874  0.0008 
#Seed - Sept_Seedling          -4.1984 0.933 32.8 -4.501  0.0004 
#Sept_Adult - Sept_Seedling    -2.4419 0.893 76.1 -2.735  0.0233 



#Mutualist


Mut_field_r1000_mod=lmer(log(value+1)~plant_data_grp+(1|gh_block),data = subset(Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_comb,grp_inter=="Mutualist"))
hist(resid(Mut_field_r1000_mod))
qqPlot(resid(Mut_field_r1000_mod))
shapiro.test(resid(Mut_field_r1000_mod))
#p-value = 0.017
plot(Mut_field_r1000_mod)

anova(Mut_field_r1000_mod)
#plant_data_grp 9.2568  1.8514     5 66.948  2.8584 0.02126 *

emmeans(Mut_field_r1000_mod,pairwise~plant_data_grp, adjust="fdr")
#Nada


#Pathogen


Path_field_r1000_mod=lmer(sqrt(value)~plant_data_grp+(1|gh_block),data = subset(Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_comb,grp_inter=="Pathogen"))
hist(resid(Path_field_r1000_mod))
qqPlot(resid(Path_field_r1000_mod))
shapiro.test(resid(Path_field_r1000_mod))
#p-value = 0.7355
plot(Path_field_r1000_mod)

anova(Path_field_r1000_mod)
#plant_data_grp 202.86  40.572     5    79  3.5953 0.00556 **

#####POSTHOC TEST: Figure 5b####
emmeans(Path_field_r1000_mod,pairwise~plant_data_grp, adjust="fdr")
#Seed - Sept_Seedling           -6.064 1.89 46.6 -3.206  0.0364 


#Cave-in-Rock Indicator

Cave_field_r1000_mod=lmer(log(value+1)~plant_data_grp+(1|gh_block),data = subset(Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_comb,grp_inter=="Cave-in-Rock Indicator"))
hist(resid(Cave_field_r1000_mod))
qqPlot(resid(Cave_field_r1000_mod))
shapiro.test(resid(Cave_field_r1000_mod))
#p-value = 7.554e-08
plot(Cave_field_r1000_mod)

anova(Cave_field_r1000_mod)
#plant_data_grp 3.0051 0.60102     5 58.623  0.9671 0.4454

emmeans(Cave_field_r1000_mod,pairwise~plant_data_grp, adjust="fdr")
#nada


#Prairie Indicator

Prair_field_r1000_mod=lmer(sqrt(value+1)~plant_data_grp+(1|gh_block),data = subset(Mar_leaf_rain.fung_decon_rar_Mar_field_pub_Endo_comb,grp_inter=="Prairie Indicator"))
hist(resid(Prair_field_r1000_mod))
qqPlot(resid(Prair_field_r1000_mod))
shapiro.test(resid(Prair_field_r1000_mod))
#p-value = 0.3005
plot(Cave_field_r1000_mod)

anova(Prair_field_r1000_mod)
#plant_data_grp 636.89  127.38     5 58.679  21.459 3.834e-12 ***

#####POSTHOC TEST: Figure 5f####
emmeans(Prair_field_r1000_mod,pairwise~plant_data_grp, adjust="fdr")
#Jul_Adult - Sept_Seedling     -7.5980 1.55 77.0 -4.895  <.0001
#Jul_Seedling - Sept_Seedling  -8.7397 1.55 77.0 -5.629  <.0001
#Jul_Seedling - Sept_Seedling  -8.7397 1.55 77.0 -5.629  <.0001
#Rain - Sept_Seedling          -6.8172 0.71 76.4 -9.598  <.0001
#Seed - Sept_Seedling          -8.1082 1.47 27.4 -5.522  <.0001 
#Sept_Adult - Sept_Seedling    -8.1645 1.37 76.0 -5.953  <.0001


#####Bray Pariwise to Sept####
Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt=read.csv(here::here("R_file","Pairwise_turnover_grad_distance_fung_v4.2.2020_decon_pruned_ZOTU_full_rar_Marshall.csv"))
head(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt)
nrow(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt)
#6786

#Subset to only seedling seedling and adult

Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_life_stage=subset(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt,
                                                                s1_s2_plant_type=="Seedling.Seedling"|s1_s2_plant_type=="Adult.Adult")
nrow(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_life_stage)
#174
unique(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_life_stage$s1_s2_collect_date)

Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_life_stage_sept=subset(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_life_stage,
                                                                     s1_s2_collect_date=="9/5/2018.7/16/2018"|
                                                                       s1_s2_collect_date=="7/16/2018.9/5/2018")
nrow(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_life_stage_sept)
#57

Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_life_stage_sept$treatment=rep("Start")


#Subset Seed to Sept

unique(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt$s1_s2_plant_type)

Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_seed1=subset(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt, 
                                                           s1_s2_plant_type=="Seed.Seedling"|s1_s2_plant_type=="Seedling.Seed"|
                                                             s1_s2_plant_type=="Seed.Adult"|s1_s2_plant_type=="Adult.Seed")
nrow(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_seed1)
#100

unique(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_seed1$s1_s2_collect_date)

Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_seed_sept=subset(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_seed1, 
                                                               s1_s2_collect_date=="9/5/2018.7/10/2018"|
                                                                 s1_s2_collect_date=="7/10/2018.9/5/2018")


nrow(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_seed_sept)
#76

Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_seed_sept$treatment=rep("Seed")


#Subset Rain to Sept

unique(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt$s1_s2_plant_type)

Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_rain1=subset(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt, 
                                                           s1_s2_plant_type=="Rain.Seedling"|s1_s2_plant_type=="Seedling.Rain"|
                                                             s1_s2_plant_type=="Rain.Adult"|s1_s2_plant_type=="Adult.Rain")
nrow(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_rain1)
#1400

unique(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_rain1$s1_s2_collect_date)

Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_rain_sept=subset(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_rain1, 
                                                               s1_collect_date=="9/5/2018"|s2_collect_date=="9/5/2018")


nrow(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_rain_sept)
#1064
Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_rain_sept$treatment=rep("Rain")


#combine the matrix 
Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_sept_fd=rbind(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_life_stage_sept,
                                                            Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_seed_sept,
                                                            Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_rain_sept)


Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_sept_fd$leaf_sample=ifelse(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_sept_fd$
                                                                           s1_collect_date=="9/5/2018",as.character(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_sept_fd$
                                                                                                                      sample1),as.character(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_sept_fd$
                                                                                                                                              sample2))

unique(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_sept_fd$leaf_sample)


Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_sept_fd$gh_block_comb=ifelse(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_sept_fd$
                                                                             s1_collect_date=="9/5/2018",as.character(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_sept_fd$
                                                                                                                        s1_gh_block),as.character(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_sept_fd$
                                                                                                                                                    s2_gh_block))

Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_sept_fd$life_stage=ifelse(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_sept_fd$
                                                                          s1_plant_type=="Seedling"|Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_sept_fd$
                                                                          s2_plant_type=="Seedling","Seedling","Adult")



head(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_sept_fd)


Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_sept_fd_sum=Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_sept_fd%>%
  group_by(leaf_sample,treatment,life_stage,gh_block_comb)%>%
  summarise_at("bray",list(~mean(.),se=~sd(.)/sqrt(n()),~n()))


ggplot(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_sept_fd_sum,aes(x=treatment,y=mean,color=life_stage,group=interaction(treatment,life_stage)))+
  geom_boxplot()

Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_sept_fd_sum1=data.frame(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_sept_fd_sum
                                                                      [,c("treatment","life_stage","mean")])

colnames(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_sept_fd_sum1)=c("treatment","life_stage","bray")
Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_sept_fd_sum_sum=Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_sept_fd_sum1%>%
  group_by(treatment,life_stage)%>%
  summarise_all(list(~mean(.),se=~sd(.)/sqrt(n()),~n(),~sd(.)))
life_order=c("Seedling","Adult")
trt_order=c("Start","Seed","Rain")
with(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_sept_fd_sum_sum,interaction(treatment,life_stage))
#brewer.pal(n = 8, name = "Paired")

(bray1000_field_p=ggplot(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_sept_fd_sum,aes(x=factor(life_stage,
                                                                                                levels = life_order),y=mean))+
    geom_point(aes(shape=interaction(factor(life_stage,levels = life_order),factor(treatment,levels = trt_order)),
                   color=interaction(factor(life_stage,levels = life_order),factor(treatment,levels = trt_order))), 
               size=4,position = position_dodge(0.75),alpha=.15)+
    geom_errorbar(data=Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_sept_fd_sum_sum, 
                  aes(group=interaction(factor(life_stage,levels = life_order),factor(treatment,levels = trt_order)),
                      x=factor(life_stage,levels = life_order),y=mean,ymin=mean-se,ymax=mean+se),width=.5,
                  position = position_dodge(0.75),color="black")+
    geom_point(data=Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_sept_fd_sum_sum,
               aes(shape=interaction(factor(life_stage,levels = life_order),factor(treatment,levels = trt_order)),
                   color=interaction(factor(life_stage,levels = life_order),factor(treatment,levels = trt_order)),
                   x=factor(life_stage,levels = life_order),y=mean), size=10,
               position = position_dodge(0.75))+ylim(0.42,0.95)+scale_shape_manual(values = c(17,15,13,13,19,19),name=NULL)+
    scale_color_manual(values = c("#B2DF8A","#A6CEE3","#E31A1C","#E31A1C","black","black"))+ylab("Bray-Curtis distance")+theme_classic()+
    theme(axis.title.x = element_blank(),axis.title.y = element_text(size = 36),axis.text = element_text(size = 30),legend.position = "none"))
#900*700

Field_pair_dist_bray1000_mod=lmer(logit(mean)~life_stage*treatment+(1|leaf_sample)+(1|gh_block_comb),
                                  data=Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_sept_fd_sum)
plot(Field_pair_dist_bray1000_mod)
hist(resid(Field_pair_dist_bray1000_mod))
qqPlot(resid(Field_pair_dist_bray1000_mod))
shapiro.test(resid(Field_pair_dist_bray1000_mod))
#p-value = 0.5964

#####Table S9 Bray-Curtis####
anova(Field_pair_dist_bray1000_mod)
#life_stage           0.07533 0.07533     1 14.315  0.6695 0.426645   
#treatment            1.39584 0.69792     2 34.000  6.2024 0.005053 **
#life_stage:treatment 0.90724 0.45362     2 34.000  4.0313 0.026848 * 
#no change

emmeans(Field_pair_dist_bray1000_mod,pairwise~life_stage|treatment)


emmeans(Field_pair_dist_bray1000_mod,pairwise~treatment)

#####POSTHOC TEST: Figure 4a####
emmeans(Field_pair_dist_bray1000_mod,pairwise~life_stage*treatment, adjust="fdr")

#


#
#####Overlap in Taxa between Field Taxa####

ntaxa(Mar_leaf_rain.fung_decon_rar_Mar_field)
#2243
nsamples(Mar_leaf_rain.fung_decon_rar_Mar_field)
#85
#what OTUs are in Seeds 
Mar_leaf_rain.fung_decon_rar_Mar_field_SEED=subset_samples(Mar_leaf_rain.fung_decon_rar_Mar_field, plant_type=="Seed")
Mar_leaf_rain.fung_decon_rar_Mar_field_SEED=prune_taxa(taxa_sums(Mar_leaf_rain.fung_decon_rar_Mar_field_SEED)>0,Mar_leaf_rain.fung_decon_rar_Mar_field_SEED)

nsamples(Mar_leaf_rain.fung_decon_rar_Mar_field_SEED)
#4
ntaxa(Mar_leaf_rain.fung_decon_rar_Mar_field_SEED)
#152

Mar_field_SEED=taxa_names(Mar_leaf_rain.fung_decon_rar_Mar_field_SEED)

#what OTUs are in Rain  
Mar_leaf_rain.fung_decon_rar_Mar_field_RAIN=subset_samples(Mar_leaf_rain.fung_decon_rar_Mar_field, plant_type=="Rain")
Mar_leaf_rain.fung_decon_rar_Mar_field_RAIN=prune_taxa(taxa_sums(Mar_leaf_rain.fung_decon_rar_Mar_field_RAIN)>0,Mar_leaf_rain.fung_decon_rar_Mar_field_RAIN)

nsamples(Mar_leaf_rain.fung_decon_rar_Mar_field_RAIN)
#56
ntaxa(Mar_leaf_rain.fung_decon_rar_Mar_field_RAIN)
#1907

Mar_field_RAIN=taxa_names(Mar_leaf_rain.fung_decon_rar_Mar_field_RAIN)


#what OTUs are in Start Seedling
Mar_leaf_rain.fung_decon_rar_Mar_field_START_SEEDLING=subset_samples(Mar_leaf_rain.fung_decon_rar_Mar_field, plant_type=="Seedling"&collect_date=="7/16/2018")
Mar_leaf_rain.fung_decon_rar_Mar_field_START_SEEDLING=prune_taxa(taxa_sums(Mar_leaf_rain.fung_decon_rar_Mar_field_START_SEEDLING)>0,Mar_leaf_rain.fung_decon_rar_Mar_field_START_SEEDLING)

nsamples(Mar_leaf_rain.fung_decon_rar_Mar_field_START_SEEDLING)
#3
ntaxa(Mar_leaf_rain.fung_decon_rar_Mar_field_START_SEEDLING)
#93

Mar_field_START_SEEDLING=taxa_names(Mar_leaf_rain.fung_decon_rar_Mar_field_START_SEEDLING)

#what OTUs are in Start Adult
Mar_leaf_rain.fung_decon_rar_Mar_field_START_ADULT=subset_samples(Mar_leaf_rain.fung_decon_rar_Mar_field, plant_type=="Adult"&collect_date=="7/16/2018")
Mar_leaf_rain.fung_decon_rar_Mar_field_START_ADULT=prune_taxa(taxa_sums(Mar_leaf_rain.fung_decon_rar_Mar_field_START_ADULT)>0,Mar_leaf_rain.fung_decon_rar_Mar_field_START_ADULT)

nsamples(Mar_leaf_rain.fung_decon_rar_Mar_field_START_ADULT)
#3
ntaxa(Mar_leaf_rain.fung_decon_rar_Mar_field_START_ADULT)
#147

Mar_field_START_ADULT=taxa_names(Mar_leaf_rain.fung_decon_rar_Mar_field_START_ADULT)

#End Seedling
Mar_leaf_rain.fung_decon_rar_Mar_field_END_SEEDLING=subset_samples(Mar_leaf_rain.fung_decon_rar_Mar_field, plant_type=="Seedling"&collect_date=="9/5/2018")
Mar_leaf_rain.fung_decon_rar_Mar_field_END_SEEDLING=prune_taxa(taxa_sums(Mar_leaf_rain.fung_decon_rar_Mar_field_END_SEEDLING)>0,Mar_leaf_rain.fung_decon_rar_Mar_field_END_SEEDLING)

nsamples(Mar_leaf_rain.fung_decon_rar_Mar_field_END_SEEDLING)
#15
ntaxa(Mar_leaf_rain.fung_decon_rar_Mar_field_END_SEEDLING)
#461
Mar_leaf_rain.fung_decon_rar_Mar_field_END_SEEDLING_overlap=data.frame(estimate_richness(Mar_leaf_rain.fung_decon_rar_Mar_field_END_SEEDLING,measures ="Observed"),
                                                                       estimate_richness(prune_taxa(Mar_field_START_SEEDLING,Mar_leaf_rain.fung_decon_rar_Mar_field_END_SEEDLING),
                                                                                         measures ="Observed"),
                                                                       estimate_richness(prune_taxa(Mar_field_SEED,Mar_leaf_rain.fung_decon_rar_Mar_field_END_SEEDLING),
                                                                                         measures ="Observed"),
                                                                       estimate_richness(prune_taxa(Mar_field_RAIN,Mar_leaf_rain.fung_decon_rar_Mar_field_END_SEEDLING),
                                                                                         measures ="Observed"),
                                                                       "overlap_start_sum"=sample_sums(prune_taxa(Mar_field_START_SEEDLING,Mar_leaf_rain.fung_decon_rar_Mar_field_END_SEEDLING)),
                                                                       "overlap_seed_sum"=sample_sums(prune_taxa(Mar_field_SEED,Mar_leaf_rain.fung_decon_rar_Mar_field_END_SEEDLING)),
                                                                       "overlap_rain_sum"=sample_sums(prune_taxa(Mar_field_RAIN,Mar_leaf_rain.fung_decon_rar_Mar_field_END_SEEDLING)),
                                                                       sample_data(Mar_leaf_rain.fung_decon_rar_Mar_field_END_SEEDLING))
colnames(Mar_leaf_rain.fung_decon_rar_Mar_field_END_SEEDLING_overlap)[1:4]=c("tot_rich","overlap_start_rich","overlap_seed_rich","overlap_rain_rich")

#End Adult

Mar_leaf_rain.fung_decon_rar_Mar_field_END_ADULT=subset_samples(Mar_leaf_rain.fung_decon_rar_Mar_field, plant_type=="Adult"&collect_date=="9/5/2018")
Mar_leaf_rain.fung_decon_rar_Mar_field_END_ADULT=prune_taxa(taxa_sums(Mar_leaf_rain.fung_decon_rar_Mar_field_END_ADULT)>0,Mar_leaf_rain.fung_decon_rar_Mar_field_END_ADULT)

nsamples(Mar_leaf_rain.fung_decon_rar_Mar_field_END_ADULT)
#4
ntaxa(Mar_leaf_rain.fung_decon_rar_Mar_field_END_ADULT)
#180

Mar_leaf_rain.fung_decon_rar_Mar_field_END_ADULT_overlap=data.frame(estimate_richness(Mar_leaf_rain.fung_decon_rar_Mar_field_END_ADULT,measures ="Observed"),
                                                                    estimate_richness(prune_taxa(Mar_field_START_ADULT,Mar_leaf_rain.fung_decon_rar_Mar_field_END_ADULT),
                                                                                      measures ="Observed"),
                                                                    estimate_richness(prune_taxa(Mar_field_SEED,Mar_leaf_rain.fung_decon_rar_Mar_field_END_ADULT),
                                                                                      measures ="Observed"),
                                                                    estimate_richness(prune_taxa(Mar_field_RAIN,Mar_leaf_rain.fung_decon_rar_Mar_field_END_ADULT),
                                                                                      measures ="Observed"),
                                                                    "overlap_start_sum"=sample_sums(prune_taxa(Mar_field_START_ADULT,Mar_leaf_rain.fung_decon_rar_Mar_field_END_ADULT)),
                                                                    "overlap_seed_sum"=sample_sums(prune_taxa(Mar_field_SEED,Mar_leaf_rain.fung_decon_rar_Mar_field_END_ADULT)),
                                                                    "overlap_rain_sum"=sample_sums(prune_taxa(Mar_field_RAIN,Mar_leaf_rain.fung_decon_rar_Mar_field_END_ADULT)),
                                                                    sample_data(Mar_leaf_rain.fung_decon_rar_Mar_field_END_ADULT))
colnames(Mar_leaf_rain.fung_decon_rar_Mar_field_END_ADULT_overlap)[1:4]=c("tot_rich","overlap_start_rich","overlap_seed_rich","overlap_rain_rich")


Mar_leaf_rain.fung_decon_rar_Mar_field_END_overlap=rbind(Mar_leaf_rain.fung_decon_rar_Mar_field_END_ADULT_overlap,Mar_leaf_rain.fung_decon_rar_Mar_field_END_SEEDLING_overlap)


Mar_leaf_rain.fung_decon_rar_Mar_field_END_overlap$overlap_start_rich_prop=Mar_leaf_rain.fung_decon_rar_Mar_field_END_overlap$overlap_start_rich/
  Mar_leaf_rain.fung_decon_rar_Mar_field_END_overlap$tot_rich
Mar_leaf_rain.fung_decon_rar_Mar_field_END_overlap$overlap_seed_rich_prop=Mar_leaf_rain.fung_decon_rar_Mar_field_END_overlap$overlap_seed_rich/
  Mar_leaf_rain.fung_decon_rar_Mar_field_END_overlap$tot_rich
Mar_leaf_rain.fung_decon_rar_Mar_field_END_overlap$overlap_rain_rich_prop=Mar_leaf_rain.fung_decon_rar_Mar_field_END_overlap$overlap_rain_rich/
  Mar_leaf_rain.fung_decon_rar_Mar_field_END_overlap$tot_rich
Mar_leaf_rain.fung_decon_rar_Mar_field_END_overlap$overlap_start_sum_prop=Mar_leaf_rain.fung_decon_rar_Mar_field_END_overlap$overlap_start_sum/1000
Mar_leaf_rain.fung_decon_rar_Mar_field_END_overlap$overlap_seed_sum_prop=Mar_leaf_rain.fung_decon_rar_Mar_field_END_overlap$overlap_seed_sum/1000
Mar_leaf_rain.fung_decon_rar_Mar_field_END_overlap$overlap_rain_sum_prop=Mar_leaf_rain.fung_decon_rar_Mar_field_END_overlap$overlap_rain_sum/1000


#Richness overlap

Mar_leaf_rain.fung_decon_rar_Mar_field_END_overlap_Rich_M=melt(Mar_leaf_rain.fung_decon_rar_Mar_field_END_overlap[c("plant_type","sampleID_bact","overlap_rain_rich_prop",
                                                                                                                    "overlap_start_rich_prop","overlap_seed_rich_prop",
                                                                                                                    "gh_block")])




Mar_leaf_rain.fung_decon_rar_Mar_field_END_overlap_Rich_M1=data.frame(Mar_leaf_rain.fung_decon_rar_Mar_field_END_overlap_Rich_M
                                                                      [,c("plant_type","variable","value")])


Mar_leaf_rain.fung_decon_rar_Mar_field_END_overlap_Rich_M_sum=Mar_leaf_rain.fung_decon_rar_Mar_field_END_overlap_Rich_M1%>%
  group_by(plant_type,variable)%>%
  summarise_all(list(~mean(.),se=~sd(.)/sqrt(n()),~n(),~sd(.)))
life_order=c("Seedling","Adult")
rich_over_order=c("overlap_start_rich_prop","overlap_seed_rich_prop","overlap_rain_rich_prop")
with(Mar_leaf_rain.fung_decon_rar_Mar_field_END_overlap_Rich_M1,interaction(variable,plant_type))
#brewer.pal(n = 8, name = "Paired")

(overlap_1000_field_p=ggplot(Mar_leaf_rain.fung_decon_rar_Mar_field_END_overlap_Rich_M,aes(x=factor(plant_type,
                                                                                                    levels = life_order),y=value))+
    geom_point(aes(shape=interaction(factor(plant_type,levels = life_order),factor(variable,levels = rich_over_order)),
                   color=interaction(factor(plant_type,levels = life_order),factor(variable,levels = rich_over_order))), 
               size=4,position = position_dodge(0.75),alpha=.15)+
    geom_errorbar(data=Mar_leaf_rain.fung_decon_rar_Mar_field_END_overlap_Rich_M_sum, 
                  aes(group=interaction(factor(plant_type,levels = life_order),factor(variable,levels = rich_over_order)),
                      x=factor(plant_type,levels = life_order),y=mean,ymin=mean-se,ymax=mean+se),width=.5,
                  position = position_dodge(0.75),color="black")+
    geom_point(data=Mar_leaf_rain.fung_decon_rar_Mar_field_END_overlap_Rich_M_sum,
               aes(shape=interaction(factor(plant_type,levels = life_order),factor(variable,levels = rich_over_order)),
                   color=interaction(factor(plant_type,levels = life_order),factor(variable,levels = rich_over_order)),
                   x=factor(plant_type,levels = life_order),y=mean), size=10,
               position = position_dodge(0.75))+scale_shape_manual(values = c(17,15,13,13,19,19),name=NULL)+
    scale_color_manual(values = c("#B2DF8A","#A6CEE3","#E31A1C","#E31A1C","black","black"))+ylab("Proportion ASV overlap")+theme_classic()+
    theme(axis.title.x = element_blank(),axis.title.y = element_text(size = 36),axis.text = element_text(size = 30),legend.position = "none"))
#900*700

#Differences in overlap with seeds
subset(Mar_leaf_rain.fung_decon_rar_Mar_field_END_overlap_Rich_M_sum, variable=="overlap_seed_rich_prop")[c("plant_type","mean")]
(0.455-0.585)/0.585


#Proportion ASV overlap

Field_taxa_overlap1000_mod=lmer((value)^-1~plant_type*variable+(1|sampleID_bact)+(1|gh_block),data=Mar_leaf_rain.fung_decon_rar_Mar_field_END_overlap_Rich_M)
#boundary (singular) fit: see ?isSingular
plot(Field_taxa_overlap1000_mod)
hist(resid(Field_taxa_overlap1000_mod))
qqPlot(resid(Field_taxa_overlap1000_mod))
shapiro.test(resid(Field_taxa_overlap1000_mod))
#p-value = 0.1163

#####Table S9 ASV Overlap####
anova(Field_taxa_overlap1000_mod)
#plant_type          0.6681  0.6681     1    17  15.008 0.0012183 ** 
#variable            7.4251  3.7125     2    34  83.404 7.725e-14 ***
#plant_type:variable 1.0835  0.5418     2    34  12.171 0.0001031 ***

emmeans(Field_taxa_overlap1000_mod,pairwise~plant_type|variable)


emmeans(Field_taxa_overlap1000_mod,pairwise~variable)


emmeans(Field_taxa_overlap1000_mod,pairwise~variable|plant_type)


#####POSTHOC TEST: Figure S10b####
emmeans(Field_taxa_overlap1000_mod,pairwise~plant_type*variable, adjust="fdr")

#Abundance overlap

Mar_leaf_rain.fung_decon_rar_Mar_field_END_overlap_Abun_M=melt(Mar_leaf_rain.fung_decon_rar_Mar_field_END_overlap[c("plant_type","sampleID_bact","overlap_rain_sum_prop",
                                                                                                                    "overlap_start_sum_prop","overlap_seed_sum_prop",
                                                                                                                    "gh_block")])




Mar_leaf_rain.fung_decon_rar_Mar_field_END_overlap_Abun_M1=data.frame(Mar_leaf_rain.fung_decon_rar_Mar_field_END_overlap_Abun_M
                                                                      [,c("plant_type","variable","value")])


Mar_leaf_rain.fung_decon_rar_Mar_field_END_overlap_Abun_M_sum=Mar_leaf_rain.fung_decon_rar_Mar_field_END_overlap_Abun_M1%>%
  group_by(plant_type,variable)%>%
  summarise_all(list(~min(.),~mean(.),se=~sd(.)/sqrt(n()),~n(),~sd(.)))
life_order=c("Seedling","Adult")
abun_over_order=c("overlap_start_sum_prop","overlap_seed_sum_prop","overlap_rain_sum_prop")
with(Mar_leaf_rain.fung_decon_rar_Mar_field_END_overlap_Rich_M1,interaction(variable,plant_type))
#brewer.pal(n = 8, name = "Paired")

(abun_overlap_1000_field_p=ggplot(Mar_leaf_rain.fung_decon_rar_Mar_field_END_overlap_Abun_M,aes(x=factor(plant_type,
                                                                                                         levels = life_order),y=value))+
    geom_point(aes(shape=interaction(factor(plant_type,levels = life_order),factor(variable,levels = abun_over_order)),
                   color=interaction(factor(plant_type,levels = life_order),factor(variable,levels = abun_over_order))), 
               size=4,position = position_dodge(0.75),alpha=.15)+
    geom_errorbar(data=Mar_leaf_rain.fung_decon_rar_Mar_field_END_overlap_Abun_M_sum, 
                  aes(group=interaction(factor(plant_type,levels = life_order),factor(variable,levels = abun_over_order)),
                      x=factor(plant_type,levels = life_order),y=mean,ymin=mean-se,ymax=mean+se),width=.5,
                  position = position_dodge(0.75),color="black")+
    geom_point(data=Mar_leaf_rain.fung_decon_rar_Mar_field_END_overlap_Abun_M_sum,
               aes(shape=interaction(factor(plant_type,levels = life_order),factor(variable,levels = abun_over_order)),
                   color=interaction(factor(plant_type,levels = life_order),factor(variable,levels = abun_over_order)),
                   x=factor(plant_type,levels = life_order),y=mean), size=10,
               position = position_dodge(0.75))+scale_shape_manual(values = c(17,15,13,13,19,19),name=NULL)+
    scale_color_manual(values = c("#B2DF8A","#A6CEE3","#E31A1C","#E31A1C","black","black"))+ylab("Proportion read overlap")+theme_classic()+
    theme(axis.title.x = element_blank(),axis.title.y = element_text(size = 36),axis.text = element_text(size = 30),legend.position = "none"))
#900*700




subset(Mar_leaf_rain.fung_decon_rar_Mar_field_END_overlap_Abun_M_sum, variable=="overlap_seed_sum_prop")[c("plant_type","mean")]
(0.521+0.648)/2

Field_abun_overlap1000_mod=lmer((value)~plant_type*variable+(1|sampleID_bact)+(1|gh_block),data=Mar_leaf_rain.fung_decon_rar_Mar_field_END_overlap_Abun_M)
#boundary (singular) fit: see ?isSingular
plot(Field_abun_overlap1000_mod)
hist(resid(Field_abun_overlap1000_mod))
qqPlot(resid(Field_abun_overlap1000_mod))
shapiro.test(resid(Field_abun_overlap1000_mod))
#p-value = 0.1975


#####Table S9 Read Overlap####
anova(Field_abun_overlap1000_mod)
#plant_type          0.08463 0.08463     1    17  4.1841 0.0566016 .  
#variable            0.73778 0.36889     2    34 18.2387 4.152e-06 ***
#plant_type:variable 0.45768 0.22884     2    34 11.3142 0.0001712 ***

emmeans(Field_abun_overlap1000_mod,pairwise~plant_type|variable)


emmeans(Field_abun_overlap1000_mod,pairwise~variable)

#####POSTHOC TEST: Figure S10d####
emmeans(Field_abun_overlap1000_mod,pairwise~plant_type*variable,adjust="fdr")



#






#####Field Indicator analyses for importance of seed versus rain####



ntaxa(Mar_leaf_rain.fung_decon_rar_Mar_field)
#2243
nsamples(Mar_leaf_rain.fung_decon_rar_Mar_field)
#85

#We only need plants
Mar_leaf_rain.fung_decon_rar_Mar_field_Rain_seed=subset_samples(Mar_leaf_rain.fung_decon_rar_Mar_field, plant_type=="Seed"|plant_type=="Rain")
nsamples(Mar_leaf_rain.fung_decon_rar_Mar_field_Rain_seed)
#60

Mar_leaf_rain.fung_decon_rar_Mar_field_Rain_seed=prune_taxa(taxa_sums(Mar_leaf_rain.fung_decon_rar_Mar_field_Rain_seed)>0,Mar_leaf_rain.fung_decon_rar_Mar_field_Rain_seed)


ntaxa(Mar_leaf_rain.fung_decon_rar_Mar_field_Rain_seed)
#1953
Mar_leaf_rain.fung_decon_rar_Mar_field_Rain_seed_otu=data.frame(as(t(otu_table(Mar_leaf_rain.fung_decon_rar_Mar_field_Rain_seed)),"matrix"))
rain_seed_order=sample_data(Mar_leaf_rain.fung_decon_rar_Mar_field_Rain_seed)$plant_type

rain_seed.ind = multipatt(Mar_leaf_rain.fung_decon_rar_Mar_field_Rain_seed_otu, 
                          rain_seed_order,duleg=TRUE, 
                          control = how(nperm=9999), func="IndVal.g")
summary(rain_seed.ind, indvalcomp=TRUE)

head(rain_seed.ind$str)
head(rain_seed.ind$sign)


#Component `A' is the probability that the surveyed
#site belongs to the target site group given the fact that the species has
#been found. This conditional probability is called the specificity or positive
#predictive value of the species as indicator of the site group.

head(rain_seed.ind$A)
#Component `B' is the probability of finding the species in sites belonging to the site group.
#This second conditional probability is called the fidelity or sensitivity of the
#species as indicator of the target site group.
head(rain_seed.ind$B)
head(rain_seed.ind$sign)


#Let's create our own indicator values based off of A and B since multipatt square root transforms and groups automatically 

rain_seed_Specif=rain_seed.ind$A
colnames(rain_seed_Specif)=c("Rain_spec","Seed_spec")
rain_seed_Fidel=rain_seed.ind$B
colnames(rain_seed_Fidel)=c("Rain_fid","Seed_fid")
rain_seed_ind_comb=merge(rain_seed_Specif,rain_seed_Fidel, by="row.names")
colnames(rain_seed_ind_comb)[colnames(rain_seed_ind_comb)=="Row.names"]="OTUs"
rain_seed_ind_comb$Rain_IndV=rain_seed_ind_comb$Rain_spec*rain_seed_ind_comb$Rain_fid
rain_seed_ind_comb$Seed_IndV=rain_seed_ind_comb$Seed_spec*rain_seed_ind_comb$Seed_fid
head(rain_seed_ind_comb)

#Let's add in the significance values for this as well 

rain_seed.ind$sign
rain_seed_ind_comb=merge(rain_seed_ind_comb,rain_seed.ind$sign, by.x = "OTUs",by.y="row.names")
head(rain_seed_ind_comb)

write.csv(rain_seed_ind_comb, here::here("R_file","IndVal_Field_rain_seed_comb_taxa_class.csv"),row.names = F)
#We need OTUs from plants only
Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS=subset_samples(Mar_leaf_rain.fung_decon_rar_Mar_field, plant_type!="Seed"&plant_type!="Rain")
nsamples(Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS)
#25

Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS=prune_taxa(taxa_sums(Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS)>0,Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS)


ntaxa(Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS)
#601

Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_otu=data.frame(otu_table(Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS))
head(Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_otu)
Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_otu$OTUs=row.names(Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_otu)
Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_otu_M=melt(Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_otu)

head(Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_otu_M)
colnames(Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_otu_M)=c("OTUs","sampleID","abund")
nrow(Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_otu_M)
#15025

#let's remove zeros to reduce the computational load
Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_otu_M_NZ=subset(Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_otu_M,abund>0)
nrow(Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_otu_M_NZ)
#2231
rain_seed_ind_comb= read.csv(here::here("R_file","IndVal_Field_rain_seed_comb_taxa_class.csv"))
Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_otu_M_IndV=merge(Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_otu_M_NZ,
                                                                  rain_seed_ind_comb,by="OTUs",all.x = T)


head(Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_otu_M_IndV)

#Now let's calculate the community weighted trait mean 
Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_otu_M_IndV$seed_IndV_abun=Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_otu_M_IndV$Seed_IndV*Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_otu_M_IndV$abund
Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_otu_M_IndV$rain_IndV_abun=Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_otu_M_IndV$Rain_IndV*Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_otu_M_IndV$abund
summary(Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_otu_M_IndV)



Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_otu_M_IndV_sum=Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_otu_M_IndV%>%group_by(sampleID)%>%summarise_at(vars(rain_IndV_abun,seed_IndV_abun),~sum(.,na.rm=TRUE))
Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_IndVal=merge(Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_otu_M_IndV_sum,
                                                              sample_data(Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS),by.x = "sampleID",
                                                              by.y = "row.names")

head(Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_IndVal)



#Make the data long so we can compare rain versus seed importance 

head(Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_IndVal)
Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_IndVal_M=melt(Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_IndVal[,c("sampleID","rain_IndV_abun",
                                                                                                                          "seed_IndV_abun","plant_type",
                                                                                                                          "gh_block",
                                                                                                                          "gh_block","collect_date")])
head(Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_IndVal_M)
Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_IndVal_M_sum=Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_IndVal_M%>%
  group_by(collect_date,plant_type,variable)%>%
  summarise_at(vars("value"),list(~mean(.),~n(),se=~sd(.)/sqrt(n())))



#Make the data long so we can compare rain versus seed importance and average
Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_IndVal$rain_IndV_abun_prop=Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_IndVal$rain_IndV_abun/1000
Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_IndVal$seed_IndV_abun_prop=Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_IndVal$seed_IndV_abun/1000
head(Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_IndVal)
Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_IndVal_M_prop=melt(Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_IndVal[,c("sampleID","rain_IndV_abun_prop",
                                                                                                                               "seed_IndV_abun_prop","plant_type",
                                                                                                                               "gh_block","collect_date")])
head(Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_IndVal_M_prop)
Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_IndVal_M_prop_sum=Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_IndVal_M_prop%>%
  group_by(collect_date,plant_type,variable)%>%
  summarise_at(vars("value"),list(~mean(.),~n(),se=~sd(.)/sqrt(n())))

with(Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_IndVal_M_sum,interaction(collect_date,plant_type))

plant_date_order_1=c("7/16/2018.Seedling","9/5/2018.Seedling","7/16/2018.Adult","9/5/2018.Adult")
brewer.pal(n = 8, name = "Paired")





(comb_IndVal_trait_nb2=ggplot(Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_IndVal_M_prop)+
    geom_point(aes(x=factor(variable,levels = c("seed_IndV_abun_prop","rain_IndV_abun_prop"), labels = c("Seed","Rain")),y=value, shape=plant_type,
                   color=factor(interaction(collect_date,plant_type),levels = plant_date_order_1)),size=4, alpha = 0.15, stroke=1,
               position = position_dodge(0.5))+
    geom_errorbar(data=Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_IndVal_M_prop_sum,
                  aes(x=factor(variable,levels = c("seed_IndV_abun_prop","rain_IndV_abun_prop"), labels = c("Seed","Rain")),
                      y=mean, ymax=mean+se,
                      ymin=mean-se,group=factor(interaction(collect_date,plant_type),levels = plant_date_order_1)),
                  width=0.5,position = position_dodge(0.5))+
    geom_point(data=Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_IndVal_M_prop_sum,
               aes(x=factor(variable,levels = c("seed_IndV_abun_prop","rain_IndV_abun_prop"), labels = c("Seed","Rain")),
                   y=mean, shape=plant_type,
                   color=factor(interaction(collect_date,plant_type),levels = plant_date_order_1)),size=10, stroke=1,position = position_dodge(0.5))+
    scale_shape_manual(values = c(15,17),name=NULL)+ylab("Average indicator value")+
    facet_grid(cols = vars(factor(plant_type,levels = c("Seedling","Adult"))),
               scales = "free_x",space = "free_x")+scale_y_continuous(limits = c(0,0.8))+
    scale_color_manual(values = c("#B2DF8A","#33A02C","#A6CEE3","#1F78B4"),name=NULL)+
    theme_cowplot()+theme(axis.title = element_blank(),axis.text.y =  element_blank(),axis.text.x = element_text(size = 30),
                          legend.position = "none",strip.text = element_text(size = 36),
                          strip.background = element_rect(fill="white",linetype = "solid",colour = "black"),
                          panel.spacing = unit(1, "lines"), panel.grid.major=element_blank(),panel.grid.minor = element_blank()))





#scale_color_manual(values = c("#E31A1C", "black"),name=NULL)+
Field_IndVal_SR_mod=lmer((value)~plant_type*collect_date*variable+
                           (1|sampleID)+(1|gh_block),data=Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_IndVal_M_prop)
plot(Field_IndVal_SR_mod)
hist(resid(Field_IndVal_SR_mod))
qqPlot(resid(Field_IndVal_SR_mod))
shapiro.test(resid(Field_IndVal_SR_mod))
#p-value =  0.6065

####Table S8 Average Indicator Value####
anova(Field_IndVal_SR_mod)

#collect_date:variable            0.284508 0.284508     1    42 25.3808 9.411e-06 ***

#####POSTHOC TEST: Figure 3b####
emmeans(Field_IndVal_SR_mod,pairwise~collect_date*variable*plant_type, adjust="fdr")
#(9/5/2018 rain_IndV_abun_prop Adult) - (9/5/2018 seed_IndV_abun_prop Seedling)       0.224376 0.0596 38.9  3.764  0.0077 
#(7/16/2018 seed_IndV_abun_prop Adult) - (9/5/2018 seed_IndV_abun_prop Seedling)      0.220768 0.0679 40.3  3.249  0.0164
#(9/5/2018 seed_IndV_abun_prop Adult) - (9/5/2018 rain_IndV_abun_prop Seedling)      -0.175210 0.0596 38.9 -2.939  0.0309 
#(9/5/2018 rain_IndV_abun_prop Seedling) - (9/5/2018 seed_IndV_abun_prop Seedling)    0.200792 0.0387 21.0  5.194  0.0011 
#(7/16/2018 seed_IndV_abun_prop Seedling) - (9/5/2018 seed_IndV_abun_prop Seedling)   0.221042 0.0679 40.3  3.253  0.0164 

emmeans(Field_IndVal_SR_mod,pairwise~collect_date*plant_type|variable, adjust="fdr")

emmeans(Field_IndVal_SR_mod,pairwise~collect_date*variable)


#####PLOT: Figure 3####
#####Plot Field and Petri IndVal combine#####

(IndVal_petri_nan_p2=ggplot(Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_IndVal_M_prop)+
   geom_point(aes(x=factor(variable,levels=c("seed_IndV_abun_prop","rain_IndV_abun_prop")),y=value,
                  shape=factor(rain_trt,levels = Rain_order),fill=factor(rain_trt,levels = Rain_order),
                  group=factor(interaction(rain_trt,collect_date),levels = rain_trt_date_order)), 
              color="black",size=4,position = position_dodge(0.65),alpha=.15)+
   geom_errorbar(data=Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_IndVal_M_prop_sum, 
                 aes(x=factor(variable,levels = c("seed_IndV_abun_prop","rain_IndV_abun_prop")),
                     group=factor(interaction(rain_trt,collect_date),levels = rain_trt_date_order),
                     y=mean,ymin=mean-se,ymax=mean+se),width=.5,
                 color="black",position = position_dodge(0.65))+scale_y_continuous(limits = c(0,0.8))+
   geom_point(data=Mar_leaf.fung_decon_pr_petr_S.rar_plants_NS_IndVal_M_prop_sum,
              aes(shape=factor(rain_trt,levels = Rain_order),group=factor(interaction(rain_trt,collect_date),levels = rain_trt_date_order),
                  fill=factor(rain_trt,levels = Rain_order),x=factor(variable,levels = c("seed_IndV_abun_prop","rain_IndV_abun_prop")),y=mean), size=10,
              color="black",position = position_dodge(0.65))+facet_wrap(vars(factor(collect_date,labels = c("Round1","Round2"))),ncol = 2)+
   scale_shape_manual(values = c(25,24),name=NULL)+scale_x_discrete(labels=c("Seed","Rain"))+
   scale_fill_manual(values = c("#6A3D9A","#1F78B4"),name=NULL)+
   ylab("Average indicator value")+theme_cowplot()+
   theme(axis.title.x = element_blank(),axis.title.y = element_text(size = 36),axis.text = element_text(size = 30),
         panel.border =  element_rect(fill = NA,linetype="solid",colour = "black"),
         legend.position = "none",strip.text = element_text(size = 36),strip.placement = "outside",strip.background = element_rect(fill = NA,linetype="blank"),
         panel.spacing = unit(0, "lines")))

(comb_IndVal_trait_nb2=ggplot(Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_IndVal_M_prop)+
    geom_point(aes(x=factor(variable,levels = c("seed_IndV_abun_prop","rain_IndV_abun_prop"), labels = c("Seed","Rain")),y=value, shape=plant_type,
                   color=factor(interaction(collect_date,plant_type),levels = plant_date_order_1)),size=4, alpha = 0.15, stroke=1,
               position = position_dodge(0.5))+
    geom_errorbar(data=Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_IndVal_M_prop_sum,
                  aes(x=factor(variable,levels = c("seed_IndV_abun_prop","rain_IndV_abun_prop"), labels = c("Seed","Rain")),
                      y=mean, ymax=mean+se,
                      ymin=mean-se,group=factor(interaction(collect_date,plant_type),levels = plant_date_order_1)),
                  width=0.5,position = position_dodge(0.5))+
    geom_point(data=Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_IndVal_M_prop_sum,
               aes(x=factor(variable,levels = c("seed_IndV_abun_prop","rain_IndV_abun_prop"), labels = c("Seed","Rain")),
                   y=mean, shape=plant_type,
                   color=factor(interaction(collect_date,plant_type),levels = plant_date_order_1)),size=10, stroke=1,position = position_dodge(0.5))+
    scale_shape_manual(values = c(15,17),name=NULL)+ylab("Average indicator value")+
    facet_grid(cols = vars(factor(plant_type,levels = c("Seedling","Adult"))),
               scales = "free_x",space = "free_x")+scale_y_continuous(limits = c(0,0.8))+
    scale_color_manual(values = c("#B2DF8A","#33A02C","#A6CEE3","#1F78B4"),name=NULL)+
    theme_cowplot()+theme(axis.title = element_blank(),axis.text.y =  element_blank(),axis.text.x = element_text(size = 30),
                          legend.position = "none",strip.text = element_text(size = 36),
                          strip.background = element_rect(fill="white",linetype = "solid",colour = "black"),
                          panel.spacing = unit(1, "lines"), panel.grid.major=element_blank(),panel.grid.minor = element_blank()))




plot_grid(IndVal_petri_nan_p2,comb_IndVal_trait_nb2,rel_widths =c(1,.8),ncol = 2, align="h",axis = "b")
#Petri_field_IndVal_RS_LFE
#1800*600

#####Significant Indicator Field analyses#####
#Extract the otu table with rain and seeds

ntaxa(Mar_leaf_rain.fung_decon_rar_Mar_field)
#2243
nsamples(Mar_leaf_rain.fung_decon_rar_Mar_field)
#85

#We need OTUs from plants only
Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS=subset_samples(Mar_leaf_rain.fung_decon_rar_Mar_field, plant_type!="Seed"&plant_type!="Rain")
nsamples(Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS)
#25

Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS=prune_taxa(taxa_sums(Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS)>0,Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS)


ntaxa(Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS)
#601

Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_otu=data.frame(otu_table(Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS))
head(Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_otu)
Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_otu$OTUs=row.names(Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_otu)

Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_otu_M=melt(Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_otu)

head(Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_otu_M)
colnames(Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_otu_M)=c("OTUs","sampleID","abund")
nrow(Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_otu_M)
#15025

#let's remove zeros to reduce the computational load
Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_otu_M_NZ=subset(Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_otu_M,abund>0)
nrow(Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_otu_M_NZ)
#2231

#Let's load in the indicator matrix and combine that with the below matrices

rain_seed_ind_comb=read.csv(here::here("R_file","IndVal_Field_rain_seed_comb_taxa_class.csv"))

Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_otu_M_IndV=merge(Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_otu_M_NZ,
                                                                  rain_seed_ind_comb,by="OTUs",all.x = T)


head(Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_otu_M_IndV)
summary(Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_otu_M_IndV)
nrow(Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_otu_M_IndV)


#Filter to only significant OTUs

Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_otu_M_IndV_sig=subset(Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_otu_M_IndV, p.value<0.05)

Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_otu_M_IndV_sig$ind_grp=
  ifelse(Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_otu_M_IndV_sig$index==1,"Rain","Seed")
summary(Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_otu_M_IndV_sig)
Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_otu_M_IndV_sig_sum=
  Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_otu_M_IndV_sig%>%group_by(sampleID,ind_grp)%>%summarise_at("abund",~sum(.))

Mar_leaf_fung_decon_rar_Mar_field_IndV_sig_map=merge(Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS_otu_M_IndV_sig_sum,
                                                     sample_data(Mar_leaf_rain.fung_decon_rar_Mar_field_plants_NS),
                                                     by.x = "sampleID",by.y = "row.names")
head(Mar_leaf_fung_decon_rar_Mar_field_IndV_sig_map)

Mar_leaf_fung_decon_rar_Mar_field_IndV_sig_map$prop_abund=Mar_leaf_fung_decon_rar_Mar_field_IndV_sig_map$abund/1000


#Graph the frequency

Mar_leaf_fung_decon_rar_Mar_field_IndV_sig_map_sum=Mar_leaf_fung_decon_rar_Mar_field_IndV_sig_map%>%
  group_by(plant_type,collect_date,ind_grp)%>%summarise_at("prop_abund",list(~mean(.),~n(),se=~sd(.)/sqrt(n())))

plant_date_order_1=c("7/16/2018.Seedling","9/5/2018.Seedling","7/16/2018.Adult","9/5/2018.Adult")




(IndVal_Field_sig_nb=ggplot(Mar_leaf_fung_decon_rar_Mar_field_IndV_sig_map)+
    geom_point(aes(x=factor(ind_grp,levels = c("Seed","Rain")),y=prop_abund, shape=plant_type,
                   color=factor(interaction(collect_date,plant_type),levels = plant_date_order_1)),size=4, alpha = 0.15, stroke=1,
               position = position_dodge(0.5))+
    geom_errorbar(data=Mar_leaf_fung_decon_rar_Mar_field_IndV_sig_map_sum,
                  aes(x=factor(ind_grp,levels = c("Seed","Rain")),
                      y=mean, ymax=mean+se,
                      ymin=mean-se,group=factor(interaction(collect_date,plant_type),levels = plant_date_order_1)),
                  width=0.5,position = position_dodge(0.5))+
    geom_point(data=Mar_leaf_fung_decon_rar_Mar_field_IndV_sig_map_sum,
               aes(x=factor(ind_grp,levels = c("Seed","Rain")),
                   y=mean, shape=plant_type,
                   color=factor(interaction(collect_date,plant_type),levels = plant_date_order_1)),size=10, stroke=1,position = position_dodge(0.5))+
    scale_shape_manual(values = c(15,17),name=NULL)+ylab("Proportion of reads indicator")+
    facet_grid(cols = vars(factor(plant_type,levels = c("Seedling","Adult"),labels = c("Seedling","Adult"))),
               scales = "free_x",space = "free_x")+scale_y_continuous(limits = c(0,0.8))+
    scale_color_manual(values = c("#B2DF8A","#33A02C","#A6CEE3","#1F78B4"),name=NULL)+
    theme_cowplot()+theme(axis.title.y = element_text(size = 34),axis.title.x = element_blank(),axis.text =  element_text(size = 30),
                          legend.position = "none",strip.text = element_text(size = 36),
                          strip.background = element_rect(fill="white",linetype = "solid",colour = "black"),
                          panel.spacing = unit(1, "lines"), panel.grid.major=element_blank(),panel.grid.minor = element_blank()))


#Stats 

Field_Ind_sig_mod=lmer(sqrt(prop_abund)~plant_type*collect_date*ind_grp+
                         (1|sampleID)+(1|gh_block),data=Mar_leaf_fung_decon_rar_Mar_field_IndV_sig_map)
plot(Field_Ind_sig_mod)
hist(resid(Field_Ind_sig_mod))
qqPlot(resid(Field_Ind_sig_mod))
shapiro.test(resid(Field_Ind_sig_mod))
#W = 0.97983, p-value = 0.5447

####Table S8 Abundance of sig Indicator Taxa####
anova(Field_Ind_sig_mod)
#collect_date:ind_grp            0.46912 0.46912     1    42 14.8688 0.0003894 ***

#####POSTHOC TEST: Figure S3b####
emmeans(Field_Ind_sig_mod,pairwise~collect_date*ind_grp*plant_type, adjust="fdr")



#####PLOT: Figure S3####
#####Plot sig IndVal Field and Petri combine#####

(IndVal_sig_petri_nan_p2=ggplot(Mar_leaf.fung_decon_pr_petr_S.rar_IndV_sig_sum_map)+
   geom_point(aes(x=factor(ind_grp,levels=c("Seed","Rain")),y=prop_abund,
                  shape=factor(rain_trt,levels = Rain_order),fill=factor(rain_trt,levels = Rain_order),
                  group=factor(interaction(rain_trt,collect_date),levels = rain_trt_date_order)), 
              color="black",size=4,position = position_dodge(0.65),alpha=.15)+
   geom_errorbar(data=Mar_leaf.fung_decon_pr_petr_S.rar_IndV_sig_sum_map_sum, 
                 aes(x=factor(ind_grp,levels = c("Seed","Rain")),
                     group=factor(interaction(rain_trt,collect_date),levels = rain_trt_date_order),
                     y=mean,ymin=mean-se,ymax=mean+se),width=.5,
                 color="black",position = position_dodge(0.65))+
   geom_point(data=Mar_leaf.fung_decon_pr_petr_S.rar_IndV_sig_sum_map_sum,
              aes(shape=factor(rain_trt,levels = Rain_order),group=factor(interaction(rain_trt,collect_date),levels = rain_trt_date_order),
                  fill=factor(rain_trt,levels = Rain_order),x=factor(ind_grp,levels = c("Seed","Rain")),y=mean), size=10,
              color="black",position = position_dodge(0.65))+facet_wrap(vars(factor(collect_date,labels = c("Round1","Round2"))),ncol = 2)+
   scale_shape_manual(values = c(25,24),name=NULL)+scale_x_discrete(labels=c("Seed","Rain"))+
   scale_fill_manual(values = c("#6A3D9A","#1F78B4"),name=NULL)+
   scale_y_continuous(name = "Proportion of\nreads indicator",limits = c(0,1))+theme_cowplot()+
   theme(axis.title.x = element_blank(),axis.title.y = element_text(size = 36),axis.text = element_text(size = 30),
         panel.border =  element_rect(fill = NA,linetype="solid",colour = "black"),
         legend.position = "none",strip.text = element_text(size = 36),strip.placement = "outside",
         strip.background = element_rect(fill="white"),
         panel.spacing = unit(0, "lines")))

(IndVal_Field_sig_nb2=ggplot(Mar_leaf_fung_decon_rar_Mar_field_IndV_sig_map)+
    geom_point(aes(x=factor(ind_grp,levels = c("Seed","Rain")),y=prop_abund, shape=plant_type,
                   color=factor(interaction(collect_date,plant_type),levels = plant_date_order_1)),size=4, alpha = 0.15, stroke=1,
               position = position_dodge(0.65))+
    geom_errorbar(data=Mar_leaf_fung_decon_rar_Mar_field_IndV_sig_map_sum,
                  aes(x=factor(ind_grp,levels = c("Seed","Rain")),
                      y=mean, ymax=mean+se,
                      ymin=mean-se,group=factor(interaction(collect_date,plant_type),levels = plant_date_order_1)),
                  width=0.5,position = position_dodge(0.65))+
    geom_point(data=Mar_leaf_fung_decon_rar_Mar_field_IndV_sig_map_sum,
               aes(x=factor(ind_grp,levels = c("Seed","Rain")),
                   y=mean, shape=plant_type,
                   color=factor(interaction(collect_date,plant_type),levels = plant_date_order_1)),size=10, stroke=1,position = position_dodge(0.65))+
    scale_shape_manual(values = c(15,17),name=NULL)+
    facet_grid(cols = vars(factor(plant_type,levels = c("Seedling","Adult"),labels = c("Seedling","Adult"))),
               scales = "free_x",space = "free_x")+scale_y_continuous(name = "Proportion of reads indicator",limits = c(0,1))+
    scale_color_manual(values = c("#B2DF8A","#33A02C","#A6CEE3","#1F78B4"),name=NULL)+
    theme_cowplot()+theme(axis.title.y = element_blank(),axis.title.x = element_blank(),axis.text.x =  element_text(size = 30),
                          axis.text.y =  element_blank(),legend.position = "none",strip.text = element_text(size = 36),
                          strip.background = element_rect(fill="white",linetype = "solid",colour = "black"),
                          panel.spacing = unit(1, "lines"), panel.grid.major=element_blank(),panel.grid.minor = element_blank()))


plot_grid(IndVal_sig_petri_nan_p2,IndVal_Field_sig_nb2,rel_widths =c(1,.75),ncol = 2, align="h",axis = "b")
#Petri_field_sig_IndVal_RS_LFE
#1800*600



#####Jaccard Field Pariwise to Sept####
Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt=read.csv(here::here("R_file","Pairwise_turnover_Jaccard_fung_v4.2.2020_decon_pruned_ZOTU_full_rar_Marshall.csv"))
head(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt)
nrow(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt)
#6786

#Subset to only seedling seedling and adult

Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_life_stage=subset(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt,
                                                                   s1_s2_plant_type=="Seedling.Seedling"|s1_s2_plant_type=="Adult.Adult")
nrow(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_life_stage)
#174
unique(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_life_stage$s1_s2_collect_date)

Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_life_stage_sept=subset(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_life_stage,
                                                                        s1_s2_collect_date=="9/5/2018.7/16/2018"|
                                                                          s1_s2_collect_date=="7/16/2018.9/5/2018")
nrow(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_life_stage_sept)
#57

Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_life_stage_sept$treatment=rep("Start")


#Subset Seed to Sept

unique(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt$s1_s2_plant_type)

Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_seed1=subset(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt, 
                                                              s1_s2_plant_type=="Seed.Seedling"|s1_s2_plant_type=="Seedling.Seed"|
                                                                s1_s2_plant_type=="Seed.Adult"|s1_s2_plant_type=="Adult.Seed")
nrow(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_seed1)
#100

unique(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_seed1$s1_s2_collect_date)

Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_seed_sept=subset(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_seed1, 
                                                                  s1_s2_collect_date=="9/5/2018.7/10/2018"|
                                                                    s1_s2_collect_date=="7/10/2018.9/5/2018")


nrow(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_seed_sept)
#76

Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_seed_sept$treatment=rep("Seed")


#Subset Rain to Sept

unique(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt$s1_s2_plant_type)

Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_rain1=subset(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt, 
                                                              s1_s2_plant_type=="Rain.Seedling"|s1_s2_plant_type=="Seedling.Rain"|
                                                                s1_s2_plant_type=="Rain.Adult"|s1_s2_plant_type=="Adult.Rain")
nrow(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_rain1)
#1400

unique(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_rain1$s1_s2_collect_date)

Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_rain_sept=subset(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_rain1, 
                                                                  s1_collect_date=="9/5/2018"|s2_collect_date=="9/5/2018")


nrow(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_rain_sept)
#1064
Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_rain_sept$treatment=rep("Rain")


#combine the matrix 
Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd=rbind(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_life_stage_sept,
                                                               Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_seed_sept,
                                                               Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_rain_sept)


Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd$leaf_sample=ifelse(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd$
                                                                              s1_collect_date=="9/5/2018",as.character(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd$
                                                                                                                         sample1),as.character(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd$
                                                                                                                                                 sample2))
unique(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd$leaf_sample)

Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd$gh_block_comb=ifelse(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd$
                                                                                s1_collect_date=="9/5/2018",as.character(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd$
                                                                                                                           s1_gh_block),as.character(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd$
                                                                                                                                                       s2_gh_block))


Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd$life_stage=ifelse(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd$
                                                                             s1_plant_type=="Seedling"|Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd$
                                                                             s2_plant_type=="Seedling","Seedling","Adult")



head(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd)


#Jaccard

Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum=Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd%>%
  group_by(leaf_sample,treatment,life_stage,gh_block_comb)%>%
  summarise_at("jaccard",list(~mean(.),se=~sd(.)/sqrt(n()),~n()))

Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum1=data.frame(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum
                                                                         [,c("treatment","life_stage","mean")])

colnames(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum1)=c("treatment","life_stage","jaccard")
Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum_sum=Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum1%>%
  group_by(treatment,life_stage)%>%
  summarise_all(list(~mean(.),se=~sd(.)/sqrt(n()),~n(),~sd(.)))
life_order=c("Seedling","Adult")
trt_order=c("Start","Seed","Rain")


(jacc1000_field_p=ggplot(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum,aes(x=factor(life_stage,
                                                                                                   levels = life_order),y=mean))+
    geom_point(aes(shape=interaction(factor(life_stage,levels = life_order),factor(treatment,levels = trt_order)),
                   color=interaction(factor(life_stage,levels = life_order),factor(treatment,levels = trt_order))), 
               size=4,position = position_dodge(0.75),alpha=.15)+
    geom_errorbar(data=Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum_sum, 
                  aes(group=interaction(factor(life_stage,levels = life_order),factor(treatment,levels = trt_order)),
                      x=factor(life_stage,levels = life_order),y=mean,ymin=mean-se,ymax=mean+se),width=.5,
                  position = position_dodge(0.75),color="black")+
    geom_point(data=Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum_sum,
               aes(shape=interaction(factor(life_stage,levels = life_order),factor(treatment,levels = trt_order)),
                   color=interaction(factor(life_stage,levels = life_order),factor(treatment,levels = trt_order)),
                   x=factor(life_stage,levels = life_order),y=mean), size=10,
               position = position_dodge(0.75))+scale_shape_manual(values = c(17,15,13,13,19,19),name=NULL)+
    scale_color_manual(values = c("#B2DF8A","#A6CEE3","#E31A1C","#E31A1C","black","black"))+ylab("Jaccard distance")+theme_classic()+
    theme(axis.title.x = element_blank(),axis.title.y = element_text(size = 36),axis.text = element_text(size = 30),legend.position = "none"))
#900*700

#Distance to seed community

subset(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum_sum, treatment=="Seed")[,c("life_stage","mean")]
(0.816-0.780)/0.780


Field_pair_dist_J1000_mod=lmer((mean)^(-1)~life_stage*treatment+(1|leaf_sample)+(1|gh_block_comb),data=Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum)
plot(Field_pair_dist_J1000_mod)
hist(resid(Field_pair_dist_J1000_mod))
qqPlot(resid(Field_pair_dist_J1000_mod))
shapiro.test(resid(Field_pair_dist_J1000_mod))
#p-value =  0.003167

#####Table S9 Ratio Jaccard####
anova(Field_pair_dist_J1000_mod)
#life_stage           0.010241 0.0102410     1 14.095  9.8526  0.007200 ** 
#treatment            0.028371 0.0141857     2 34.000 13.6476 4.456e-05 ***
#life_stage:treatment 0.011715 0.0058576     2 34.000  5.6354  0.007695 ** 

#####POSTHOC TEST: Figure 4b####
emmeans(Field_pair_dist_J1000_mod,pairwise~life_stage*treatment,adjust="fdr")
emmeans(Field_pair_dist_J1000_mod,pairwise~life_stage|treatment)
emmeans(Field_pair_dist_J1000_mod,pairwise~treatment|life_stage)




emmeans(Field_pair_dist_J1000_mod,pairwise~treatment)


#Nestedness
head(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd)
Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum_nest=Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd%>%
  group_by(leaf_sample,treatment,life_stage,gh_block_comb)%>%
  summarise_at("nestedness",list(~mean(.),se=~sd(.)/sqrt(n()),~n()))
Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum_nest1=data.frame(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum_nest
                                                                              [,c("treatment","life_stage","mean")])
colnames(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum_nest1)=c("treatment","life_stage","nestedness")
Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum_nest_sum=Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum_nest1%>%
  group_by(treatment,life_stage)%>%
  summarise_all(list(~mean(.),se=~sd(.)/sqrt(n()),~n(),~sd(.)))
life_order=c("Seedling","Adult")
trt_order=c("Start","Seed","Rain")

#brewer.pal(n = 8, name = "Paired")

(nested1000_field_p=ggplot(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum_nest,aes(x=factor(life_stage,
                                                                                                          levels = life_order),y=mean))+
    geom_point(aes(shape=interaction(factor(life_stage,levels = life_order),factor(treatment,levels = trt_order)),
                   color=interaction(factor(life_stage,levels = life_order),factor(treatment,levels = trt_order))), 
               size=4,position = position_dodge(0.75),alpha=.15)+
    geom_errorbar(data=Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum_nest_sum, 
                  aes(group=interaction(factor(life_stage,levels = life_order),factor(treatment,levels = trt_order)),
                      x=factor(life_stage,levels = life_order),y=mean,ymin=mean-se,ymax=mean+se),width=.5,
                  position = position_dodge(0.75),color="black")+
    geom_point(data=Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum_nest_sum,
               aes(shape=interaction(factor(life_stage,levels = life_order),factor(treatment,levels = trt_order)),
                   color=interaction(factor(life_stage,levels = life_order),factor(treatment,levels = trt_order)),
                   x=factor(life_stage,levels = life_order),y=mean), size=10,
               position = position_dodge(0.75))+scale_shape_manual(values = c(17,15,13,13,19,19),name=NULL)+
    scale_color_manual(values = c("#B2DF8A","#A6CEE3","#E31A1C","#E31A1C","black","black"))+ylab("Nestedness distance")+theme_classic()+
    theme(axis.title.x = element_blank(),axis.title.y = element_text(size = 36),axis.text = element_text(size = 30),legend.position = "none"))
#900*700

#Distance to Rain community

subset(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum_nest_sum, treatment=="Rain")[,c("life_stage","mean")]
(0.108-0.181)/0.181


Field_pair_dist_nest1000_mod=lmer(sqrt(mean)~life_stage*treatment+(1|leaf_sample)+(1|gh_block_comb),data=Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum_nest)
rePCA(Field_pair_dist_nest1000_mod)
plot(Field_pair_dist_nest1000_mod)
hist(resid(Field_pair_dist_nest1000_mod))
qqPlot(resid(Field_pair_dist_nest1000_mod))
shapiro.test(resid(Field_pair_dist_nest1000_mod))
#p-value = 0.2216


#####Table S9 Nestedness####
anova(Field_pair_dist_nest1000_mod)
#life_stage           0.004956 0.004956     1 48.227  1.0046 0.3211911    
#treatment            0.029414 0.014707     2 48.102  2.9811 0.0601866 .  
#life_stage:treatment 0.085589 0.042795     2 48.102  8.6744 0.0006069 ***


#####POSTHOC TEST: Figure S10a####
emmeans(Field_pair_dist_nest1000_mod,pairwise~life_stage*treatment,adjust="fdr")

emmeans(Field_pair_dist_nest1000_mod,pairwise~life_stage|treatment)

emmeans(Field_pair_dist_nest1000_mod,pairwise~treatment|life_stage)
emmeans(Field_pair_dist_nest1000_mod,pairwise~treatment)

#Turnover

Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum_turn=Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd%>%
  group_by(leaf_sample,treatment,life_stage,gh_block_comb)%>%
  summarise_at("turnover",list(~mean(.),se=~sd(.)/sqrt(n()),~n()))
Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum_turn1=data.frame(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum_turn
                                                                              [,c("treatment","life_stage","mean")])

colnames(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum_turn1)=c("treatment","life_stage","turnover")
Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum_turn_sum=Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum_turn1%>%
  group_by(treatment,life_stage)%>%
  summarise_all(list(~mean(.),se=~sd(.)/sqrt(n()),~n(),~sd(.)))
life_order=c("Seedling","Adult")
trt_order=c("Start","Seed","Rain")

#brewer.pal(n = 8, name = "Paired")

(turnover1000_field_p=ggplot(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum_turn,aes(x=factor(life_stage,
                                                                                                            levels = life_order),y=mean))+
    geom_point(aes(shape=interaction(factor(life_stage,levels = life_order),factor(treatment,levels = trt_order)),
                   color=interaction(factor(life_stage,levels = life_order),factor(treatment,levels = trt_order))), 
               size=4,position = position_dodge(0.75),alpha=.15)+
    geom_errorbar(data=Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum_turn_sum, 
                  aes(group=interaction(factor(life_stage,levels = life_order),factor(treatment,levels = trt_order)),
                      x=factor(life_stage,levels = life_order),y=mean,ymin=mean-se,ymax=mean+se),width=.5,
                  position = position_dodge(0.75),color="black")+
    geom_point(data=Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum_turn_sum,
               aes(shape=interaction(factor(life_stage,levels = life_order),factor(treatment,levels = trt_order)),
                   color=interaction(factor(life_stage,levels = life_order),factor(treatment,levels = trt_order)),
                   x=factor(life_stage,levels = life_order),y=mean), size=10,
               position = position_dodge(0.75))+scale_shape_manual(values = c(17,15,13,13,19,19),name=NULL)+
    scale_color_manual(values = c("#B2DF8A","#A6CEE3","#E31A1C","#E31A1C","black","black"))+ylab("Turnover distance")+theme_classic()+
    theme(axis.title.x = element_blank(),axis.title.y = element_text(size = 36),axis.text = element_text(size = 30),legend.position = "none"))
#900*700

#Distance to Rain community

subset(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum_turn_sum, treatment=="Rain")[,c("life_stage","mean")]
(0.730-0.634)/0.634

Field_pair_dist_turn1000_mod=lmer((mean)~life_stage*treatment+(1|leaf_sample)+(1|gh_block_comb),data=Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum_turn)
plot(Field_pair_dist_turn1000_mod)
hist(resid(Field_pair_dist_turn1000_mod))
qqPlot(resid(Field_pair_dist_turn1000_mod))
shapiro.test(resid(Field_pair_dist_turn1000_mod))
#p-value = 0.2782


#####Table S9 Turnover####
anova(Field_pair_dist_turn1000_mod)
#life_stage           0.0049799 0.0049799     1 14.086  1.7042 0.21266  
#treatment            0.0001026 0.0000513     2 34.000  0.0176 0.98260  
#life_stage:treatment 0.0197438 0.0098719     2 34.000  3.3784 0.04589 *

emmeans(Field_pair_dist_turn1000_mod,pairwise~life_stage|treatment)


#####POSTHOC TEST: Figure S10c####
emmeans(Field_pair_dist_turn1000_mod,pairwise~life_stage*treatment,adjust="fdr")

emmeans(Field_pair_dist_turn1000_mod,pairwise~treatment|life_stage)
#plot_grid(turnover1000_field_p,nested1000_field_p,jacc1000_field_p,ncol = 1)

#800x1800


#plot_grid(nested1000_field_p,turnover1000_field_p,nrow = 2, align="v")
#900*1000

#####PLOT: Figure S10####
#####Plot Nestedness Turnover and overlap combine####

(nested1000_field_p=ggplot(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum_nest,aes(x=factor(life_stage,
                                                                                                          levels = life_order),y=mean))+
   geom_point(aes(shape=interaction(factor(life_stage,levels = life_order),factor(treatment,levels = trt_order)),
                  color=interaction(factor(life_stage,levels = life_order),factor(treatment,levels = trt_order))), 
              size=4,position = position_dodge(0.75),alpha=.15)+
   geom_segment(y=subset(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum_nest_sum,treatment=="Seed"&life_stage=="Adult")$mean,
                yend=subset(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum_nest_sum,treatment=="Rain"&life_stage=="Adult")$mean,
                x=2,xend=2.25, color="#1F78B4",size=1)+
   geom_segment(y=subset(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum_nest_sum,treatment=="Seed"&life_stage=="Seedling")$mean,
                yend=subset(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum_nest_sum,treatment=="Rain"&life_stage=="Seedling")$mean,
                x=1,xend=1.25, color="#33A02C" , size=1)+
   geom_errorbar(data=Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum_nest_sum, 
                 aes(group=interaction(factor(life_stage,levels = life_order),factor(treatment,levels = trt_order)),
                     x=factor(life_stage,levels = life_order),y=mean,ymin=mean-se,ymax=mean+se),width=.5,
                 position = position_dodge(0.75),color="black")+
   geom_point(data=Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum_nest_sum,
              aes(shape=interaction(factor(life_stage,levels = life_order),factor(treatment,levels = trt_order)),
                  color=interaction(factor(life_stage,levels = life_order),factor(treatment,levels = trt_order)),
                  x=factor(life_stage,levels = life_order),y=mean), size=10,
              position = position_dodge(0.75))+scale_shape_manual(values = c(17,15,13,13,19,19),name=NULL)+
   scale_color_manual(values = c("#B2DF8A","#A6CEE3","#E31A1C","#E31A1C","black","black"))+
   scale_y_continuous(name = "Nestedness distance",limits = c(0,0.30))+theme_cowplot()+
   theme(axis.title.x = element_blank(),axis.title.y = element_text(size = 36),axis.text = element_text(size = 30),
         axis.text.x = element_blank(),legend.position = "none"))


(overlap_1000_field_p=ggplot(Mar_leaf_rain.fung_decon_rar_Mar_field_END_overlap_Rich_M,aes(x=factor(plant_type,
                                                                                                    levels = life_order),y=value))+
    geom_point(aes(shape=interaction(factor(plant_type,levels = life_order),factor(variable,levels = rich_over_order)),
                   color=interaction(factor(plant_type,levels = life_order),factor(variable,levels = rich_over_order))), 
               size=4,position = position_dodge(0.75),alpha=.15)+
    geom_segment(y=subset(Mar_leaf_rain.fung_decon_rar_Mar_field_END_overlap_Rich_M_sum,variable=="overlap_seed_rich_prop"&plant_type=="Adult")$mean,
                 yend=subset(Mar_leaf_rain.fung_decon_rar_Mar_field_END_overlap_Rich_M_sum,variable=="overlap_rain_rich_prop"&plant_type=="Adult")$mean,
                 x=2,xend=2.25, color="#1F78B4",size=1)+
    geom_segment(y=subset(Mar_leaf_rain.fung_decon_rar_Mar_field_END_overlap_Rich_M_sum,variable=="overlap_seed_rich_prop"&plant_type=="Seedling")$mean,
                 yend=subset(Mar_leaf_rain.fung_decon_rar_Mar_field_END_overlap_Rich_M_sum,variable=="overlap_rain_rich_prop"&plant_type=="Seedling")$mean,
                 x=1,xend=1.25, color="#33A02C" , size=1)+
    geom_errorbar(data=Mar_leaf_rain.fung_decon_rar_Mar_field_END_overlap_Rich_M_sum, 
                  aes(group=interaction(factor(plant_type,levels = life_order),factor(variable,levels = rich_over_order)),
                      x=factor(plant_type,levels = life_order),y=mean,ymin=mean-se,ymax=mean+se),width=.5,
                  position = position_dodge(0.75),color="black")+
    geom_point(data=Mar_leaf_rain.fung_decon_rar_Mar_field_END_overlap_Rich_M_sum,
               aes(shape=interaction(factor(plant_type,levels = life_order),factor(variable,levels = rich_over_order)),
                   color=interaction(factor(plant_type,levels = life_order),factor(variable,levels = rich_over_order)),
                   x=factor(plant_type,levels = life_order),y=mean), size=10,
               position = position_dodge(0.75))+scale_shape_manual(values = c(17,15,13,13,19,19),name=NULL)+
    scale_color_manual(values = c("#B2DF8A","#A6CEE3","#E31A1C","#E31A1C","black","black"))+
    scale_y_continuous(name = "Proportion ASV overlap",breaks = c(0.25,0.50,0.75,1),limits = c(0.15,1))+theme_cowplot()+
    theme(axis.title.x = element_blank(),axis.title.y = element_text(size = 36),axis.text = element_text(size = 30),
          axis.text.x = element_blank(),legend.position = "none"))


(turnover1000_field_p=ggplot(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum_turn,aes(x=factor(life_stage,
                                                                                                            levels = life_order),y=mean))+
    geom_point(aes(shape=interaction(factor(life_stage,levels = life_order),factor(treatment,levels = trt_order)),
                   color=interaction(factor(life_stage,levels = life_order),factor(treatment,levels = trt_order))), 
               size=4,position = position_dodge(0.75),alpha=.15)+
    geom_segment(y=subset(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum_turn_sum,treatment=="Seed"&life_stage=="Adult")$mean,
                 yend=subset(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum_turn_sum,treatment=="Rain"&life_stage=="Adult")$mean,
                 x=2,xend=2.25, color="#1F78B4",size=1)+
    geom_segment(y=subset(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum_turn_sum,treatment=="Seed"&life_stage=="Seedling")$mean,
                 yend=subset(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum_turn_sum,treatment=="Rain"&life_stage=="Seedling")$mean,
                 x=1,xend=1.25, color="#33A02C" , size=1)+
    geom_errorbar(data=Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum_turn_sum, 
                  aes(group=interaction(factor(life_stage,levels = life_order),factor(treatment,levels = trt_order)),
                      x=factor(life_stage,levels = life_order),y=mean,ymin=mean-se,ymax=mean+se),width=.5,
                  position = position_dodge(0.75),color="black")+
    geom_point(data=Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum_turn_sum,
               aes(shape=interaction(factor(life_stage,levels = life_order),factor(treatment,levels = trt_order)),
                   color=interaction(factor(life_stage,levels = life_order),factor(treatment,levels = trt_order)),
                   x=factor(life_stage,levels = life_order),y=mean), size=10,
               position = position_dodge(0.75))+scale_shape_manual(values = c(17,15,13,13,19,19),name=NULL)+
    scale_color_manual(values = c("#B2DF8A","#A6CEE3","#E31A1C","#E31A1C","black","black"))+
    scale_y_continuous(name = "Turnover distance")+theme_cowplot()+
    theme(axis.title.x = element_blank(),axis.title.y = element_text(size = 36),axis.text = element_text(size = 30),
          axis.text.x = element_text(size = 36),legend.position = "none"))


(abun_overlap_1000_field_p=ggplot(Mar_leaf_rain.fung_decon_rar_Mar_field_END_overlap_Abun_M,aes(x=factor(plant_type,
                                                                                                         levels = life_order),y=value))+
    geom_point(aes(shape=interaction(factor(plant_type,levels = life_order),factor(variable,levels = abun_over_order)),
                   color=interaction(factor(plant_type,levels = life_order),factor(variable,levels = abun_over_order))), 
               size=4,position = position_dodge(0.75),alpha=.15)+
    geom_segment(y=subset(Mar_leaf_rain.fung_decon_rar_Mar_field_END_overlap_Abun_M_sum,variable=="overlap_seed_sum_prop"&plant_type=="Adult")$mean,
                 yend=subset(Mar_leaf_rain.fung_decon_rar_Mar_field_END_overlap_Abun_M_sum,variable=="overlap_rain_sum_prop"&plant_type=="Adult")$mean,
                 x=2,xend=2.25, color="#1F78B4",size=1)+
    geom_segment(y=subset(Mar_leaf_rain.fung_decon_rar_Mar_field_END_overlap_Abun_M_sum,variable=="overlap_seed_sum_prop"&plant_type=="Seedling")$mean,
                 yend=subset(Mar_leaf_rain.fung_decon_rar_Mar_field_END_overlap_Abun_M_sum,variable=="overlap_rain_sum_prop"&plant_type=="Seedling")$mean,
                 x=1,xend=1.25, color="#33A02C" , size=1)+
    geom_errorbar(data=Mar_leaf_rain.fung_decon_rar_Mar_field_END_overlap_Abun_M_sum, 
                  aes(group=interaction(factor(plant_type,levels = life_order),factor(variable,levels = abun_over_order)),
                      x=factor(plant_type,levels = life_order),y=mean,ymin=mean-se,ymax=mean+se),width=.5,
                  position = position_dodge(0.75),color="black")+
    geom_point(data=Mar_leaf_rain.fung_decon_rar_Mar_field_END_overlap_Abun_M_sum,
               aes(shape=interaction(factor(plant_type,levels = life_order),factor(variable,levels = abun_over_order)),
                   color=interaction(factor(plant_type,levels = life_order),factor(variable,levels = abun_over_order)),
                   x=factor(plant_type,levels = life_order),y=mean), size=10,
               position = position_dodge(0.75))+scale_shape_manual(values = c(17,15,13,13,19,19),name=NULL)+
    scale_color_manual(values = c("#B2DF8A","#A6CEE3","#E31A1C","#E31A1C","black","black"))+
    scale_y_continuous(name = "Proportion read overlap",breaks = c(0.25,0.50,0.75,1),limits = c(0.15,1))+theme_cowplot()+
    theme(axis.title.x = element_blank(),axis.title.y = element_text(size = 36),axis.text.y = element_text(size = 30),
          axis.text.x = element_text(size = 36),legend.position = "none"))



plot_grid(nested1000_field_p,overlap_1000_field_p,turnover1000_field_p,abun_overlap_1000_field_p, 
          rel_heights = c(0.9,1))
#1500*1200
#Nest_Turn_OverLap_Field_dot_error_comb_raw



#####Field Nestedness to Turnover Ratio ####
Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd$ratio_nest_turn=Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd$nestedness/Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd$turnover



head(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd)
Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum_rat_NT=Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd%>%
  group_by(leaf_sample,treatment,life_stage,gh_block_comb)%>%
  summarise_at("ratio_nest_turn",list(~mean(.),se=~sd(.)/sqrt(n()),~n()))
Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum_rat_NT1=data.frame(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum_rat_NT
                                                                                [,c("treatment","life_stage","mean")])
colnames(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum_rat_NT1)=c("treatment","life_stage","ratio_nest_turn")
Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum_rat_NT_sum=Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum_rat_NT1%>%
  group_by(treatment,life_stage)%>%
  summarise_all(list(~mean(.),se=~sd(.)/sqrt(n()),~n(),~sd(.)))
life_order=c("Seedling","Adult")
trt_order=c("Start","Seed","Rain")

#brewer.pal(n = 8, name = "Paired")

(rat_NT_1000_field_p=ggplot(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum_rat_NT,aes(x=factor(life_stage,
                                                                                                             levels = life_order),y=mean))+
    geom_point(aes(shape=interaction(factor(life_stage,levels = life_order),factor(treatment,levels = trt_order)),
                   color=interaction(factor(life_stage,levels = life_order),factor(treatment,levels = trt_order))), 
               size=4,position = position_dodge(0.75),alpha=.15)+
    geom_errorbar(data=Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum_rat_NT_sum, 
                  aes(group=interaction(factor(life_stage,levels = life_order),factor(treatment,levels = trt_order)),
                      x=factor(life_stage,levels = life_order),y=mean,ymin=mean-se,ymax=mean+se),width=.5,
                  position = position_dodge(0.75),color="black")+
    geom_point(data=Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum_rat_NT_sum,
               aes(shape=interaction(factor(life_stage,levels = life_order),factor(treatment,levels = trt_order)),
                   color=interaction(factor(life_stage,levels = life_order),factor(treatment,levels = trt_order)),
                   x=factor(life_stage,levels = life_order),y=mean), size=10,
               position = position_dodge(0.75))+scale_shape_manual(values = c(17,15,13,13,19,19),name=NULL)+scale_x_discrete(labels=c("End Seedlings","End Adults"))+
    scale_color_manual(values = c("#B2DF8A","#A6CEE3","#E31A1C","#E31A1C","black","black"))+ylab("Ratio of\nNestedness to Turnover")+theme_classic()+
    theme(axis.title.x = element_blank(),axis.title.y = element_text(size = 36),axis.text = element_text(size = 30),legend.position = "none"))
#900*700
#dot_Ratio_Nest_Turn_Field_dot_error_comb_raw

#####PLOT: Figure 4####
#####Plot field Pair Distance and Ratio####



(bray1000_field_p2=ggplot(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_sept_fd_sum,aes(x=factor(life_stage,
                                                                                                 levels = life_order),y=mean))+
    geom_point(aes(shape=interaction(factor(life_stage,levels = life_order),factor(treatment,levels = trt_order)),
                   color=interaction(factor(life_stage,levels = life_order),factor(treatment,levels = trt_order))), 
               size=4,position = position_dodge(0.75),alpha=.15)+
   geom_segment(y=subset(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_sept_fd_sum_sum,treatment=="Seed"&life_stage=="Adult")$mean,
                yend=subset(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_sept_fd_sum_sum,treatment=="Rain"&life_stage=="Adult")$mean,
                x=2,xend=2.25, color="#1F78B4",size=1)+
   geom_segment(y=subset(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_sept_fd_sum_sum,treatment=="Seed"&life_stage=="Seedling")$mean,
                yend=subset(Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_sept_fd_sum_sum,treatment=="Rain"&life_stage=="Seedling")$mean,
                x=1,xend=1.25, color="#33A02C" , size=1)+
    geom_errorbar(data=Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_sept_fd_sum_sum, 
                  aes(group=interaction(factor(life_stage,levels = life_order),factor(treatment,levels = trt_order)),
                      x=factor(life_stage,levels = life_order),y=mean,ymin=mean-se,ymax=mean+se),width=.5,
                  position = position_dodge(0.75),color="black")+
    geom_point(data=Mar_leaf_rain.fung_decon_rar_Mar.betapair_trt_sept_fd_sum_sum,
               aes(shape=interaction(factor(life_stage,levels = life_order),factor(treatment,levels = trt_order)),
                   color=interaction(factor(life_stage,levels = life_order),factor(treatment,levels = trt_order)),
                   x=factor(life_stage,levels = life_order),y=mean), size=10,
               position = position_dodge(0.75))+scale_shape_manual(values = c(17,15,13,13,19,19),name=NULL)+
    scale_color_manual(values = c("#B2DF8A","#A6CEE3","#E31A1C","#E31A1C","black","black"))+
   scale_y_continuous(name = "Bray-Curtis distance", limits = c(0.42,0.95),breaks = c(seq(0.3,1,by=0.15)))+theme_cowplot()+
    theme(axis.title.y = element_text(size = 36),axis.text.y =  element_text(size = 30),
          axis.title.x = element_blank(),axis.text.x =  element_blank(),legend.position = "none"))




(jacc1000_field_p2=ggplot(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum,aes(x=factor(life_stage,
                                                                                                    levels = life_order),y=mean))+
    geom_point(aes(shape=interaction(factor(life_stage,levels = life_order),factor(treatment,levels = trt_order)),
                   color=interaction(factor(life_stage,levels = life_order),factor(treatment,levels = trt_order))), 
               size=4,position = position_dodge(0.75),alpha=.15)+
    geom_segment(y=subset(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum_sum,treatment=="Seed"&life_stage=="Adult")$mean,
                 yend=subset(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum_sum,treatment=="Rain"&life_stage=="Adult")$mean,
                 x=2,xend=2.25, color="#1F78B4",size=1)+
    geom_segment(y=subset(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum_sum,treatment=="Seed"&life_stage=="Seedling")$mean,
                 yend=subset(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum_sum,treatment=="Rain"&life_stage=="Seedling")$mean,
                 x=1,xend=1.25, color="#33A02C" , size=1)+
    geom_errorbar(data=Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum_sum, 
                  aes(group=interaction(factor(life_stage,levels = life_order),factor(treatment,levels = trt_order)),
                      x=factor(life_stage,levels = life_order),y=mean,ymin=mean-se,ymax=mean+se),width=.5,
                  position = position_dodge(0.75),color="black")+
    geom_point(data=Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum_sum,
               aes(shape=interaction(factor(life_stage,levels = life_order),factor(treatment,levels = trt_order)),
                   color=interaction(factor(life_stage,levels = life_order),factor(treatment,levels = trt_order)),
                   x=factor(life_stage,levels = life_order),y=mean), size=10,
               position = position_dodge(0.75))+scale_shape_manual(values = c(17,15,13,13,19,19),name=NULL)+
    scale_color_manual(values = c("#B2DF8A","#A6CEE3","#E31A1C","#E31A1C","black","black"))+
    scale_y_continuous(name = "Jaccard distance", breaks=c(seq(0.6,1,by=0.05)))+theme_cowplot()+
    theme(axis.title.y = element_text(size = 36),axis.text.y =  element_text(size = 30),
          axis.title.x = element_blank(),axis.text.x =  element_blank(),legend.position = "none"))


(rat_NT_1000_field_p=ggplot(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum_rat_NT)+
    geom_point(aes(x=factor(life_stage,
                                levels = life_order),y=mean,shape=interaction(factor(life_stage,levels = life_order),
                                                                              factor(treatment,levels = trt_order)),
                   color=interaction(factor(life_stage,levels = life_order),factor(treatment,levels = trt_order))), 
               size=4,position = position_dodge(0.75),alpha=.15)+
    geom_segment(y=subset(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum_rat_NT_sum,treatment=="Seed"&life_stage=="Adult")$mean,
                 yend=subset(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum_rat_NT_sum,treatment=="Rain"&life_stage=="Adult")$mean,
                 x=2,xend=2.25, color="#1F78B4",size=1)+
    geom_segment(y=subset(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum_rat_NT_sum,treatment=="Seed"&life_stage=="Seedling")$mean,
                 yend=subset(Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum_rat_NT_sum,treatment=="Rain"&life_stage=="Seedling")$mean,
                 x=1,xend=1.25, color="#33A02C" , size=1)+
    geom_errorbar(data=Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum_rat_NT_sum, 
                  aes(group=interaction(factor(life_stage,levels = life_order),factor(treatment,levels = trt_order)),
                      x=factor(life_stage,levels = life_order),y=mean,ymin=mean-se,ymax=mean+se),width=.5,
                  position = position_dodge(0.75),color="black")+
    geom_point(data=Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum_rat_NT_sum,
               aes(shape=interaction(factor(life_stage,levels = life_order),factor(treatment,levels = trt_order)),
                   color=interaction(factor(life_stage,levels = life_order),factor(treatment,levels = trt_order)),
                   x=factor(life_stage,levels = life_order),y=mean), size=10,
               position = position_dodge(0.75))+scale_shape_manual(values = c(17,15,13,13,19,19),name=NULL)+scale_x_discrete(labels=c("End Seedlings","End Adults"))+
    scale_color_manual(values = c("#B2DF8A","#A6CEE3","#E31A1C","#E31A1C","black","black"))+
    scale_y_continuous(name = "Ratio of\nNestedness to Turnover",breaks=c(seq(0,1,by=0.15)))+theme_cowplot()+
    theme(axis.title.x = element_blank(),axis.title.y = element_text(size = 36),axis.text.x = element_text(size = 36),
          axis.text.y = element_text(size = 30),legend.position = "none"))


plot_grid(bray1000_field_p2,jacc1000_field_p2,rat_NT_1000_field_p,nrow = 3, align="v",rel_heights = c(0.8,0.8,1))
#Field_Bray_Jac_Ratio_dot_raw
#1000*1600








#Ratio Nestedness to Turnover Stats

Field_pair_dist_rat_NT_1000_mod=lmer(sqrt(mean)~life_stage*treatment+(1|leaf_sample)+(1|gh_block_comb),data=Mar_leaf_rain.fung_decon_rar_Mar_PA.betapair_trt_sept_fd_sum_rat_NT)
rePCA(Field_pair_dist_rat_NT_1000_mod)
plot(Field_pair_dist_rat_NT_1000_mod)
hist(resid(Field_pair_dist_rat_NT_1000_mod))
qqPlot(resid(Field_pair_dist_rat_NT_1000_mod))
shapiro.test(resid(Field_pair_dist_rat_NT_1000_mod))
#p-value = 0.2216


#####Table S9 Ratio Nestedness:Turnover####
anova(Field_pair_dist_rat_NT_1000_mod)
#life_stage           0.003811 0.003811     1 48.173  0.3432 0.5607328    
#treatment            0.053334 0.026667     2 48.082  2.4013 0.1013829    
#life_stage:treatment 0.183267 0.091634     2 48.082  8.2516 0.0008301 ***

#####POSTHOC TEST: Figure 4c####
emmeans(Field_pair_dist_rat_NT_1000_mod,pairwise~life_stage*treatment,adjust="fdr")

emmeans(Field_pair_dist_rat_NT_1000_mod,pairwise~life_stage|treatment)

emmeans(Field_pair_dist_rat_NT_1000_mod,pairwise~treatment|life_stage)
emmeans(Field_pair_dist_rat_NT_1000_mod,pairwise~treatment)

