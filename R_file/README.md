# Here are the files generated for the analysis and characterization of rain, Switchgrass leaf and seed endophytic fungi. 

#### Matched_pub_seq_names_TBL.csv
##### File describing the functional role of each published leaf fungal endophyte.

+ target - the NCBI classification of each published OTU
+ Study - study in which OTU was charaterized and published (see below for citation) 
+ Taxa - low classification of OTU from study
+ Result - description of results on characterization of OTU
+ Group - simplified groupin of each OTU based of published results
+ Interaction - combined interaction group limited to Mutualist, Pathogen, Context Mutualist, and Unknown

##### Sequences and functional roles are based off of 
+ Giauque, H. & Hawkes, C.V. 2013 Climate affects symbiotic fungal endophyte diversity and performance. American Journal of Botany 100, 1435-1444. (doi:10.3732/ajb.1200568).
+ Kleczewski, N.M., Bauer, J.T., Bever, J.D., Clay, K. & Reynolds, H.L. 2012 A survey of endophytic fungi of switchgrass (Panicum virgatum) in the Midwest, and their putative roles in plant growth. Fungal Ecology 5, 521-529. (doi:https://doi.org/10.1016/j.funeco.2011.12.006).
+ Giauque, H. 2016 Hierarchical controls of endophyte-mediated drought tolerance : ecological, physiological, and molecular. 
+ Whitaker, B.K., Reynolds, H.L. & Clay, K. 2018 Foliar fungal endophyte communities are structured by environment but not host ecotype in Panicum virgatum (switchgrass). Ecology 99, 2703-2711. (doi:10.1002/ecy.2543).

#### rep_set.FunRainLeaf2019_ZOTU_full.fung_cons_UNKNOWN_UKNOWN_v4.2.2020.fna
##### Representative sequence of OTUs that were not well classified by CONSTAX 

#### rep_set.FunRainLeaf2019_v4.2.2020_ZOTU_rar1000.fna
##### Representative sequence of OTUs used in analyses for the manuscript

#### Mar_leaf_rain_ZOTU_full_v4.2.2020.fung_decon_rar_pruned_phyloseq_obj.RData
##### Phyloseq object with the rarified OTU table, mapping data, and Taxonmy used for the analyses. 

#### Mar_leaf.fung_v4.2.2020_decon_pr_petr_S_nan.rar_phyloseq_obj.RData
##### Phyloseq object with the rarified OTU table, mapping data, and Taxonmy subsetted to only include data for the petri experiment
##### Generated in the "Petri community analyses No Nano" portion of the code


#### Mar_leaf_rain.fung_v4.2.2020_decon_rar_Mar_field_phyloseq_obj.RData
##### Phyloseq object with the rarified OTU table, mapping data, and Taxonmy subsetted to only include data for the field experiment
##### Generated in the "Field experiment analyses" portion of the code

#### Treatments_Mar_leaf.fung_decon_pr_petr_S_nan.rar_map.csv
##### Treatments file created for the PERMANOVA analyses the petri experiment in Primer v6

#### Bray_Mar_leaf.fung_decon_pr_petr_S_nan.rar_dis.csv
##### Bray-Curtis distance matrix created for the PERMANOVA analyses the petri experiment in Primer v6

#### Jaccard_Mar_leaf.fung_decon_pr_petr_S_nan.rar_J_dis.csv
##### Jaccard distance matrix created for the PERMANOVA analyses the petri experiment in Primer v6

#### Mar_leaf_rain.fung_decon_rar_Mar_field_map.csv
##### Treatments file created for the PERMANOVA analyses the field experiment in Primer v6


#### Bray_Mar_leaf_rain.fung_decon_rar_Mar_field_dis.csv
##### Bray-Curtis distance matrix created for the PERMANOVA analyses the petri experiment in Primer v6

#### Jaccard_Mar_leaf_rain.fung_decon_rar_Mar_field_J_dis.csv
##### Jaccard distance matrix created for the PERMANOVA analyses the field experiment in Primer v6

Pairwise_turnover_Jaccard_fung_v4.2.2020_decon_pruned_ZOTU_full_rar_Marshall.csv
Pairwise_turnover_grad_distance_fung_v4.2.2020_decon_pruned_ZOTU_full_rar_Marshall.csv


## Outputs from the comparison of OTUs to previously published OTUs
#### Code run for the matching of OTUs

```

#The program needs for the codons to be uppercase
awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' rep_set.FunRainLeaf2019_v4.2.2020_ZOTU_rar1000.fna > CAP_rep_set.FunRainLeaf2019_v4.2.2020_ZOTU_rar1000.fna
~/HardDrive/Sciencey_Program/usearch11 -usearch_global CAP_rep_set.FunRainLeaf2019_v4.2.2020_ZOTU_rar1000.fna -db switchgrass_compiled_ITS_database_v2.fasta -id 0.97 -strand both -maxaccepts 0 -maxhits 10 -matched MATCHED_pub_seq_rep_set.FunRainLeaf2019_ZOTU_v4.2.2020_rar1000.fa -notmatched NOT_matched_pub_seq_rep_set.FunRainLeaf2019_ZOTU_v4.2.2020_rar1000.fa -userout TBL_pub_seq_rep_set.FunRainLeaf2019_ZOTU_v4.2.2020_rar1000.txt -userfields query+target+id+mid+bits+evalue+ql+ts+qlor+qhir+tlor+thir
#00:02 22Mb    100.0% Searching, 5.2% matched
```

#### CAP_rep_set.FunRainLeaf2019_v4.2.2020_ZOTU_rar1000.fna
##### Representative sequence of OTUs used in analyses for the manuscript with codons capitalized for usearch_global

#### TBL_pub_seq_rep_set.FunRainLeaf2019_ZOTU_v4.2.2020_rar1000.txt
##### Table of matching OTUs see [user field](https://www.drive5.com/usearch/manual/userfields.html) for the description of columns
#### NOT_matched_pub_seq_rep_set.FunRainLeaf2019_ZOTU_v4.2.2020_rar1000.fa
##### Representative sequence of OTUs that did not match the published sequences
#### MATCHED_pub_seq_rep_set.FunRainLeaf2019_ZOTU_v4.2.2020_rar1000.fa
##### Representative sequence of OTUs that matched the published sequences


## Outputs from the classification of OTUs to [FunGuild](https://github.com/UMNFuN/FUNGuild) accessed March 2021
#### FunRainLeaf2019_ZOTU_CONSTAXv4.2.2020_ZOTU_for_FunGuild.csv
##### OTU Taxonomy formated for the FunGuilf classificaiton. See "FUNGuild CONSTAX Dataset creation" section of R code

```
python Guilds_v1.1.py -otu FunRainLeaf2019_ZOTU_CONSTAXv4.2.2020_ZOTU_for_FunGuild.csv -db fungi -m -u

#Found 12894 matching taxonomy records in the database.
#Dereplicating and sorting the result...
#FunGuild tried to assign function to 15504 OTUs in 'FunRainLeaf2019_ZOTU_CONSTAXv4.2.2020_ZOTU_for_FunGuild.csv'.
#FUNGuild made assignments on 8307 OTUs.
#Result saved to 'FunRainLeaf2019_ZOTU_CONSTAXv4.2.2020_ZOTU_for_FunGuild.guilds.txt'

#Additional output:
#  FUNGuild made assignments on 8307 OTUs, these have been saved to FunRainLeaf2019_ZOTU_CONSTAXv4.2.2020_ZOTU_for_FunGuild.guilds_matched.txt.
#7197 OTUs were unassigned, these are saved to FunRainLeaf2019_ZOTU_CONSTAXv4.2.2020_ZOTU_for_FunGuild.guilds_unmatched.txt.

#Total calculating time: 185.96 seconds.
```

#### Files generated by FunGuild

+FunRainLeaf2019_ZOTU_CONSTAXv4.2.2020_ZOTU_for_FunGuild.guilds.txt
+FunRainLeaf2019_ZOTU_CONSTAXv4.2.2020_ZOTU_for_FunGuild.guilds_matched.txt
F+unRainLeaf2019_ZOTU_CONSTAXv4.2.2020_ZOTU_for_FunGuild.guilds_unmatched.txt








