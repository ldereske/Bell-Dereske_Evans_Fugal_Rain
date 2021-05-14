# Raw files from the USEARCHv11 pipeline used to process the raw sequences

## Files made outside of USEARCH 

### Marshall_petri_exp_plant_data.csv
Metadata for mapping the samples sequenced and processed for analyses 

+ sampleID_bact	- unique identfier for samples
+ sampleID_fung - unique identfier for samples used for labeling the raw sequences
+ sample_type - extraction type. "Leaf" were extracted with DNeasy Plant kits, "Rain" was extracted using DNeasy PowerWater Kit, "Blank" were PCR blanks of molecular water, and "Mock" was AMPtk mock community grown and extracted by Bonito lab.
+ Sample_or_Control	- identifer of Sample (plant, Rain, or mock) or control (kit blanks and PCR blanks)
+ extract_num - order in which DNA extraction was conducted.
+ extract_kit - DNA extraction kit used for each sample (Water_1 = lot 160019846, Water_2 = lot 160013097, Water_3 = lot 160021106, Plant_1 = lot 160038217, and Plant_2 = lot 160038217)
+ fungal_PCR_conc - DNA concentration of PCR product amplified using ITS9-ITS4 submitted for squencing. Quantified using Qubit Fluorometric Quantification on a plate reader.
+ raw_DNA_conc - DNA concentration of raw DNA post extraction. Quantified using NanoDrop
+ bact_dilution	- DNA dilution amounts for 16S sequencing not included in this study
+ dil_DNA_conc - DNA dilution amounts for DNA templates used for PCR amplified using ITS9-ITS4
+ description - long unique sample id
+ start_date - date of the start of the sample collection or experimental manipulation
+ collect_date - date of the collection of samples
+ short_description	- unique sample descriptor with project and extraction number
+ wet_mass - wet mass of sample before surface sterilization and DNA extraction
+ vol_filter - volume of rain vacuum filtered for DNA extraction
+ sub_proj - project that each sample was taken from. (GLBRC = Great Lake Bioenergy Research Center Intensive Site, Marshall = Marshall scale-up experiment, and Blank = not associate with any project)
+ plant_type - type of sample in terms of plant and experimental type.
+ field_plot - field pot from which sample was taken 
+ gh_block - combination of field plot and greenhouse block that each sample was from
+ rain_trt - rain treatment for the petri experiment (Live_Rain = unmanipulated rain, Sterile_rain = autoclaved sterilized rain, and Nano = autoclaved sterilized nanopure water (included as another control but not included in analyses))
+ pot_num - pot number or field block for each sample
+ UTM_lat - estimate latitude of the sampling location in UTM
+ UTM_long - estimate longitude of the sampling location in UTM
+ exp_len - length of sampleing
+ seq_plate_no - the plate that each sample was included on for sequencing
+ seq_plate_row	- row that each sample was on the sequecing plate
+ seq_plate_col - column that each sample was on the sequecing plate


### Marshall_petri_exp_plant_data.csv
Data from the plant response of the petri experiment

+ dish - full sample identifier for each sample
+ field_plot - field block from which rain was collected at Marshall Farms
+ pot_num - replicated with in each field plot X rain treatment combination
+ rain_trt - rain treatment for the petri experiment (Live_Rain = unmanipulated rain, Sterile_rain = autoclaved sterilized rain, and Nano = autoclaved sterilized nanopure water (included as another control but not included in analyses))
+ germ_seed	- number of seeds that germinated in a petri dish
+ fung_seed	- number of seeds with visual fungal colonization in each petri dish
+ plant_description - description of plants that were extracted from each petri dish
+ measure_date - date that each petri dish was quanitified
+ start_date - date that the rain sampling begun
+ collect_date - date that the rain sampling was finished
+ note - notes from observer 


## ITS2 USEARCH pipeline Fungal sequences in Rain and Leaves

### Folder with the extracted forward and reverse reads where I need to extract reads
``` 
cd /mnt/research/EvansLab/FungiRainLeaf2019/20191202_Amplicon_PE250
gunzip *.fastq.gz
```
## 1) Quality checking
### 1a) First look at the quality of  raw unmerged seqs for run1
#https://www.drive5.com/usearch/manual/pipe_readprep_understand.html

```
mkdir fastq_info_FungRainLeaf
```

#### make a forloop to run fastx_info on every file
```
nano fasta_info_fq.sh
!#/bin/bash
for fq in *.fastq
do
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -fastx_info $fq -output fastq_info_FungRainLeaf/$fq
done
```

#### make file executable and run the for loop create in your `fasta_info_fq.sh` file
```
chmod +x fasta_info_fq.sh
./fasta_info_fq.sh
```

#### move to the fastq_info directory
```
cd fastq_info_FungRainLeaf/
```

#### Now run the below code to summarize the fastq info for all of the forward and reverse reads

```
grep "^File" * > FungRainLeaf_fastq_lengths.txt
grep "^EE" * > FungRainLeaf_fastq_EE.txt
```

look for any forward and reverse reads that look especially bad in terms of quality (high E is bad quality). This info will also be really helpful for troubleshooting later on (e.g. why some samples have extremely low read numbers)





##  2) Merge the forward and reverse sequences and trim adapters (for each run individually)
### 2a) Merge Pairs
#### Make sure you are in the folder with the extracted forward and reverse reads from Run 1
https://www.drive5.com/usearch/manual/merge_options.html
#### -alnout gives you a human readable text file of the alignment and misalignments for each pair merged.
#### -tabbedout give you extensive information on the quality of the merge
#### This step takes approximately 1 minute
```
cd /mnt/research/EvansLab/Lukas
mkdir ITS2_FungiRainLeaf2019
cd ITS2_FungiRainLeaf2019
mkdir mergedfastq_FungiRainLeaf2019

cd /mnt/research/EvansLab/FungiRainLeaf2019/20191202_Amplicon_PE250
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -fastq_mergepairs *_R1*.fastq -relabel @ -fastqout /mnt/research/EvansLab/Lukas/ITS2_FungiRainLeaf2019/mergedfastq_FungiRainLeaf2019/ITS2_merged_FungiRainLeaf2019.fastq  -tabbedout /mnt/research/EvansLab/Lukas/ITS2_FungiRainLeaf2019/mergedfastq_FungiRainLeaf2019/combined_merged_FungiRainLeaf2019_pair_report.txt -alnout /mnt/research/EvansLab/Lukas/ITS2_FungiRainLeaf2019/mergedfastq_FungiRainLeaf2019/combined_merged_FungiRainLeaf2019_pair_aln.txt
```
```
05:08 720Mb   100.0% 82.2% merged

Totals:
   8601322  Pairs (8.6M)
   7068346  Merged (7.1M, 82.18%)
   4700068  Alignments with zero diffs (54.64%)
   1361203  Too many diffs (> 5) (15.83%)
    171773  No alignment found (2.00%)
         0  Alignment too short (< 16) (0.00%)
    145883  Staggered pairs (1.70%) merged & trimmed
    130.38  Mean alignment length
    362.09  Mean merged length
      0.69  Mean fwd expected errors
      0.61  Mean rev expected errors
      0.32  Mean merged expected errors
```
### 2b) let's check sequence quality of the merged seqs using USEARCH's [fastq_eestats2](https://www.drive5.com/usearch/manual/cmd_fastq_eestats2.html).

```
cd /mnt/research/EvansLab/Lukas/ITS2_FungiRainLeaf2019/mergedfastq_FungiRainLeaf2019
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -fastq_eestats2 ITS2_merged_FungiRainLeaf2019.fastq -output ITS2_merged_FungiRainLeaf2019_eestats2.txt
```
```
7068346 reads, max len 484, avg 362.1

Length         MaxEE 0.50         MaxEE 1.00         MaxEE 2.00
------   ----------------   ----------------   ----------------
    50    6969684( 98.6%)    6994424( 99.0%)    6994517( 99.0%)
   100    6668029( 94.3%)    6919158( 97.9%)    6937881( 98.2%)
   150    6417245( 90.8%)    6842859( 96.8%)    6933313( 98.1%)
   200    6252217( 88.5%)    6757554( 95.6%)    6915271( 97.8%)
   250    6146677( 87.0%)    6694725( 94.7%)    6895115( 97.5%)
   300    5887479( 83.3%)    6519607( 92.2%)    6793079( 96.1%)
   350    3197092( 45.2%)    3684850( 52.1%)    3961786( 56.0%)
   400     990541( 14.0%)    1182958( 16.7%)    1305219( 18.5%)
   450      62327(  0.9%)      94312(  1.3%)     124148(  1.8%)
```



## 3) Let's run the USEARCH version of phix removal
#https://www.drive5.com/usearch/manual/cmd_filter_phix.html

```
mkdir no_phix_USEARCH

/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -filter_phix ITS2_merged_FungiRainLeaf2019.fastq -output no_phix_USEARCH/TS2_merged_FungiRainLeaf2019_filtered.fq -alnout no_phix_USEARCH/TS2_merged_FungiRainLeaf2019_phix_hits.txt
#01:05 875Mb   100.0% Filtering for phix, 77 hits (0.0%)

/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -fastq_eestats2 no_phix_USEARCH/TS2_merged_FungiRainLeaf2019_filtered.fq -output no_phix_USEARCH/TS2_merged_FungiRainLeaf2019_filtered_eestats2.txt
```
```
7068269 reads, max len 484, avg 362.1

Length         MaxEE 0.50         MaxEE 1.00         MaxEE 2.00
------   ----------------   ----------------   ----------------
    50    6969616( 98.6%)    6994347( 99.0%)    6994440( 99.0%)
   100    6667982( 94.3%)    6919089( 97.9%)    6937804( 98.2%)
   150    6417201( 90.8%)    6842796( 96.8%)    6933237( 98.1%)
   200    6252176( 88.5%)    6757491( 95.6%)    6915195( 97.8%)
   250    6146642( 87.0%)    6694665( 94.7%)    6895039( 97.5%)
   300    5887473( 83.3%)    6519575( 92.2%)    6793010( 96.1%)
   350    3197092( 45.2%)    3684847( 52.1%)    3961749( 56.0%)
   400     990541( 14.0%)    1182958( 16.7%)    1305217( 18.5%)
   450      62327(  0.9%)      94312(  1.3%)     124148(  1.8%)
```
### Let's check to see if there are primers or adapters in the fastq file
need to create a random subset to test the frequency of primer sequences
```
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -fastx_subsample no_phix_USEARCH/TS2_merged_FungiRainLeaf2019_filtered.fq -sample_pct 1 -fastaout no_phix_USEARCH/ITS2_merged_FungiRainLeaf2019_filtered_one_per.fq
```
```
#make a fasta file with adapters and primers
nano JGI_CS1-2_ITS2_primers.fa
>ITS9F  
ACACTGACGACATGGTTCTACAGAACGCAGCRAAIIGYGA
>ITS4Rev
TACGGTAGCAGAGACTTGGTCTTCCTCCGCTTATTGATATGC
```
```
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -search_oligodb no_phix_USEARCH/ITS2_merged_FungiRainLeaf2019_filtered_one_per.fq -db JGI_CS1-2_ITS2_primers.fa -strand both -userout no_phix_USEARCH/TS2_merged_FungiRainLeaf2019_filtered_primer_adapter_out_one_per.txt -userfields query+target+diffs+tlo+thi+qlor+qhir
#00:06 667Mb   100.0% Searching ITS2_merged_FungiRainLeaf2019_filtered_one_per.fq, 0.0% matched
```

### Double checking that primers are removed by just search for the adapters
```
#make a fasta file with adapter sequences
nano CS1-2_adapters.fa
>CS1
ACACTGACGACATGGTTCTACA
>CS2
TACGGTAGCAGAGACTTGGTCT 

/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -search_oligodb no_phix_USEARCH/ITS2_merged_FungiRainLeaf2019_filtered_one_per.fq -db CS1-2_adapters.fa -strand both -userout no_phix_USEARCH/TS2_merged_FungiRainLeaf2019_filtered_primer_out_one_per.txt -userfields query+target+diffs+tlo+thi+qlor+qhir

```
### Triple checking that primers are removed by just search for the primers themselves 
```
#make a fasta file with primers
nano JGI_ITS2_primers.fa
>ITS9F  
GAACGCAGCRAAIIGYGA
>ITS4Rev
TCCTCCGCTTATTGATATGC
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -search_oligodb no_phix_USEARCH/ITS2_merged_FungiRainLeaf2019_filtered_one_per.fq -db JGI_ITS2_primers.fa -strand both -userout no_phix_USEARCH/TS2_merged_FungiRainLeaf2019_filtered_primer_only_out_one_per.txt -userfields query+target+diffs+tlo+thi+qlor+qhir
#00:20 668Mb   100.0% Searching ITS2_merged_FungiRainLeaf2019_filtered_one_per.fq, 99.8% matched
```
### Let's use the USEARCH primer remover
#https://www.drive5.com/usearch/manual/cmd_fastx_trim_primer.html
From the user output created using -search_oligodb it looks like the primer sequences start at 17-20 bp
```
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -fastx_trim_primer no_phix_USEARCH/ITS2_merged_FungiRainLeaf2019_filtered_one_per.fq -db JGI_ITS2_primers.fa -strand both  -maxdiffs 5 -width 21 -fastaout no_phix_USEARCH/ITS2_merged_FungiRainLeaf2019_phix_filtered_one_per_primer_re.fa -tabbedout no_phix_USEARCH/ITS2_merged_FungiRainLeaf2019_phix_filtered_one_per_primer_remove.txt

#How many reads were removed
grep -c "^>" no_phix_USEARCH/ITS2_merged_FungiRainLeaf2019_phix_filtered_one_per_primer_re.fa
#70648 reads

grep -c "^>" no_phix_USEARCH/ITS2_merged_FungiRainLeaf2019_filtered_one_per.fq
#70682
```

Before we continue, you may want to check if the sample names are formatted correctly. USEARCH does some funny cutting during the merging step. Any hyphens or underscores can be problematic and you need to remove these (use sed command and merged_cut files)

Additionally, this is a good opportunity to double check that all of your samples merged and have unique IDs using [fastx_get_sample_names](https://www.drive5.com/usearch/manual/cmd_fastx_get_sample_names.html)

```
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -fastx_get_sample_names no_phix_USEARCH/TS2_merged_FungiRainLeaf2019_filtered.fq -output ITS2_merged_FungiRainLeaf2019_phix_filtered_samples.txt
#00:27 38Mb    100.0% 200 samples found
```



## 4) Filtering and Truncate the merged seqs  to MaxEE and set length using [fastq_filter](https://www.drive5.com/usearch/manual/cmd_fastq_filter.html)

#### I need to first filter by quality
```
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -fastq_filter no_phix_USEARCH/TS2_merged_FungiRainLeaf2019_filtered.fq -fastq_maxee 1 -fastaout ITS2_merged_FungiRainLeaf2019_phix_filtered_full.fa
```
```
01:15 628Mb   100.0% Filtering, 93.4% passed
   7068269  Reads (7.1M)
    467864  Discarded reads with expected errs > 1.00
   6600405  Filtered reads (6.6M, 93.4%)

```
#### Now let's try to remove the primers and all reads that lack the primer squence as a final filtering step
```
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -fastx_trim_primer ITS2_merged_FungiRainLeaf2019_phix_filtered_full.fa -db JGI_ITS2_primers.fa -strand both  -maxdiffs 2 -width 21 -fastaout ITS2_merged_FungiRainLeaf2019_phix_primer_filtered_full.fa -tabbedout ITS2_merged_FungiRainLeaf2019_phix_filtered_primer_full_remove.txt
#00:00 44Mb    100.0% Reading JGI_ITS2_primers.fa
#01:01 44Mb    100.0% Processing

grep -c "^>" ITS2_merged_FungiRainLeaf2019_phix_primer_filtered_full.fa
#6588968
```
#### Now I want to truncate the sequences to 350 bp

```

/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -fastx_truncate ITS2_merged_FungiRainLeaf2019_phix_primer_filtered_full.fa -trunclen 350 -padlen 350 -fastaout ITS2_merged_FungiRainLeaf2019_phix_primer_filtered_full350.fa
```
## 5) Filter so we only have unique sequences with [fastx_uniques](https://www.drive5.com/usearch/manual/cmd_fastx_uniques.html)
```
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -fastx_uniques ITS2_merged_FungiRainLeaf2019_phix_primer_filtered_full350.fa  -fastaout uniques_ITS2_merged_FungiRainLeaf2019_phix_primer_filtered_full.fa -sizeout
```
```
00:29 4.7Gb  6588968 seqs, 1412499 uniques, 982842 singletons (69.6%)
00:29 4.7Gb  Min size 1, median 1, max 95096, avg 4.66
00:47 4.1Gb   100.0% Writing uniques_ITS2_merged_FungiRainLeaf2019_phix_primer_filtered_full.fa

```

## 6) Cluster into OTUS and filter out singletons
There are two options here. **(A)** uses the traditional approach and clusters sequences into 0.97 identity cutoff OTUs. **(B)** uses unoise3 to identify ZOTUs.

### 6A) Cluster into 0.97 OTUs using UPARSE and [cluster_otus](https://www.drive5.com/usearch/manual/cmd_cluster_otus.html)
This step will also denovo chimera check and filter out singletons. You can remove single sequences prior to clustering but singletons are also removed at the OTU clustering step (defaults cluster_otus filters out OTUs <2 and unoise3 filters ZOTUs <8)
### 6B) Identify ZOTUs using [unoise3](https://www.drive5.com/usearch/manual/cmd_unoise3.html)
This step will also denovo chimera check and filter out low abundance ZOTUs. IMPORTANT: ZOTUs with less than 3 reads will be filtered out (i.e. -minsize 3) default setting filters out ZOTUs less than 8 reads
Let's use the SLURM job submission

```
nano cluster_FungiRainLeaf2019_Z.OTU_full.sbatch

#!/bin/bash --login
########## SBATCH Lines for Resource Request ##########
 
#SBATCH --time=3:30:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1-2                # number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=2                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=5           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem-per-cpu=20G            # memory required per allocated CPU (or core) - amount of memory (in bytes)
#SBATCH --job-name   cluster_FungiRainLeaf2019_Z.OTU  # you can give your job a name for easier identification (same as -J)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=belldere@msu.edu

cd /mnt/research/EvansLab/Lukas/ITS2_FungiRainLeaf2019/mergedfastq_FungiRainLeaf2019

/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -cluster_otus uniques_ITS2_merged_FungiRainLeaf2019_phix_primer_filtered_full.fa -otus rep_set_ITS2_full_FungiRainLeaf2019_otus.fa -uparseout ITS2_full_FungiRainLeaf2019_otus_uparse.txt -relabel OTU

/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -unoise3 uniques_ITS2_merged_FungiRainLeaf2019_phix_primer_filtered_full.fa -zotus rep_set_ITS2_full_FungiRainLeaf2019_zotus.fa  -tabbedout ITS2_full_FungiRainLeaf2019_otus_zotus_report.txt -minsize 3 

###end of .sbatch

sbatch cluster_FungiRainLeaf2019_Z.OTU_full.sbatch
#Submitted batch job 60224204
```
```
#OTU clustering results
04:20 77Mb    100.0% 7095 OTUs, 3204 chimeras

#ZOTU denoise and filtering results
56:10 745Mb   100.0% 18466 good, 26 chimeras
56:11 745Mb   100.0% Writing zotus
```
## 7) Map reads back to OTUs at a 97% similarity score using [otutab](https://www.drive5.com/usearch/manual/cmd_otutab.html)
**-id 0.97 -strand plus are defaults**
### 7A) Mapping reads to traditional 0.97 OTUS
### 7B) Mapping reads to ZOTUs

```
nano mapping_FungiRainLeaf2019_Z.OTU_full.sbatch

#!/bin/bash --login
########## SBATCH Lines for Resource Request ##########
 
#SBATCH --time=20:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1-2                # number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=2                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=5           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem-per-cpu=20G            # memory required per allocated CPU (or core) - amount of memory (in bytes)
#SBATCH --job-name   mapping_FungiRainLeaf2019_Z.OTU  # you can give your job a name for easier identification (same as -J)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=belldere@msu.edu

cd /mnt/research/EvansLab/Lukas/ITS2_FungiRainLeaf2019/mergedfastq_FungiRainLeaf2019

#OTU
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -otutab ITS2_merged_FungiRainLeaf2019.fastq -otus rep_set_ITS2_full_FungiRainLeaf2019_otus.fa -uc ITS2_full_FungiRainLeaf2019_OTU_map.uc -otutabout OTU_table_ITS2_full_FungiRainLeaf2019.txt -biomout OTU_table_ITS2_full_FungiRainLeaf2019_jsn.biom -notmatchedfq ITS2_full_merged_FungiRainLeaf2019_otu_unmapped.fq

#ZOTU
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -otutab ITS2_merged_FungiRainLeaf2019.fastq -zotus rep_set_ITS2_full_FungiRainLeaf2019_zotus.fa -uc ITS2_FungiRainLeaf2019_ZOTU_map.uc -otutabout ZOTU_table_ITS2_full_FungiRainLeaf2019.txt -biomout ZOTU_table_ITS2_full_FungiRainLeaf2019_jsn.biom -notmatchedfq ITS2_full_merged_FungiRainLeaf2019_ZOTU_unmapped.fq

gzip ITS2_full_FungiRainLeaf2019_ZOTU_map.uc
gzip ITS2_full_FungiRainLeaf2019_OTU_map.uc

gzip ITS2_full_merged_FungiRainLeaf2019_otu_unmapped.fq
gzip ITS2_full_merged_FungiRainLeaf2019_ZOTU_unmapped.fq

###end of .sbatch

sbatch mapping_FungiRainLeaf2019_Z.OTU_full.sbatch
#Submitted batch job 60227843
```
```
#OTU mapping
6934765 / 7068346 mapped to OTUs (98.1%)        
04:39:49 339Mb  Writing OTU_table_ITS2_full_FungiRainLeaf2019.txt
04:39:49 339Mb  Writing OTU_table_ITS2_full_FungiRainLeaf2019.txt ...done.
04:39:50 339Mb  Writing OTU_table_ITS2_full_FungiRainLeaf2019_jsn.biom
04:39:50 339Mb  Writing OTU_table_ITS2_full_FungiRainLeaf2019_jsn.biom ...done.

#ZOTU mapping
6988341 / 7068346 mapped to OTUs (98.9%)        
02:54:31 363Mb  Writing ZOTU_table_ITS2_full_FungiRainLeaf2019.txt
02:54:31 363Mb  Writing ZOTU_table_ITS2_full_FungiRainLeaf2019.txt ...done.
02:54:31 363Mb  Writing ZOTU_table_ITS2_full_FungiRainLeaf2019_jsn.biom
02:54:31 363Mb  Writing ZOTU_table_ITS2_full_FungiRainLeaf2019_jsn.biom ...done.

```

## 8) Classifying taxa against the reference database using [sintax](https://www.drive5.com/usearch/manual/cmd_sintax.html)


### 8A) Classifying the traditional 0.97 OTUs against UNITE8.2 v04.02.2020
### 8B) Classifying the ZOTUs against UNITE8.2 v04.02.2020


```
cd /mnt/research/EvansLab/Lukas/ITS2_FungiRainLeaf2019/mergedfastq_FungiRainLeaf2019

/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -sintax rep_set_ITS2_full_FungiRainLeaf2019_zotus.fa -db /mnt/ufs18/home-087/belldere/Databases/UNITe8.2/sh_general_release_all_04.02.2020/CONSTAX_v2_Nav_training_files/sintax.db -tabbedout taxonomy_ITS2_full_FungiRainLeaf2019_zotus_v04.02.2020.sintax -strand both

/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -sintax rep_set_ITS2_full_FungiRainLeaf2019_zotus.fa -db  /mnt/ufs18/home-087/belldere/Databases/UNITe8.2/sh_general_release_all_04.02.2020/utax_reference_dataset_all_04.02.2020_fix.udb -tabbedout taxonomy_ITS2_full_FungiRainLeaf2019_zotus_v04.02.2020_UTAX.sintax -strand both

```



### 9) Classifying the taxa with the native install and no Bonito lab paths and the utax_reference_dataset_all_04.02.2020_fix.udb [CONSTAX](https://github.com/Gian77/CONSTAXv2) against the eukaryote UNITe 8.2
```

cd /mnt/home/belldere/ITS2_FungiRainLeaf2019

nano CONSTAX_v2_sintax_fix_all_class_FungiRainLeaf2019_ZOTU_full.sbatch

#!/bin/bash --login

#SBATCH --time=12:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --ntasks=1                 # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=1          # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem=32G                  # memory required per node - amount of memory (in bytes)
#SBATCH --job-name costax_v2_sintax_fix_All_v4.2.2020_T         # you can give your job a name for easier identification (same as -J)
#SBATCH --output=%x-%j.SLURMout
#SBATCH --mail-type=ALL
#SBATCH --mail-user=belldere@msu.edu

cd /mnt/home/belldere/ITS2_FungiRainLeaf2019
export PATH=$PATH:$HOME/anaconda3/bin
source activate CONSTAXv2
/mnt/home/belldere/Programs/CONSTAX/CONSTAXv2/constax.sh --input /mnt/research/EvansLab/Lukas/ITS2_FungiRainLeaf2019/mergedfastq_FungiRainLeaf2019/rep_set_ITS2_full_FungiRainLeaf2019_zotus.fasta --db /mnt/ufs18/home-087/belldere/Databases/UNITe8.2/sh_general_release_all_04.02.2020/sh_general_release_dynamic_all_04.02.2020_fix.fasta --trainfile /mnt/home/belldere/Databases/UNITe8.2/sh_general_release_all_04.02.2020/CONSTAX_v2_Nav_training_files/ --tax ZOTU_constax_V2_sintax_fix_classification_all_v4.2.2020_taxa/ -o ZOTU_constax_V2_sintax_fix_classification_all_v4.2.2020_taxa/ --conf 0.8 -b --pathfile /mnt/home/belldere/Programs/CONSTAX/CONSTAXv2/pathfile.txt 

scontrol show job $SLURM_JOB_ID     ### write job information to output file
####

sbatch CONSTAX_v2_sintax_fix_all_class_FungiRainLeaf2019_ZOTU_full.sbatch

#Submitted batch job 6657328
```