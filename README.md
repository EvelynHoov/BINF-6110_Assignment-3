# BINF-6110_Assignment-3
Comparison of gut microbiota between individuals maintaining vegan and omnivorous diets.
# Introduction  
  Gut microbiota diversity is not only influenced by how someone lives, but also incredibly influential to their wellbeing. For instance, the health of an individual's microbiome is an incredibly influential factor to the development of chronic diseases (Singh et al., 2017). For this reason, characterizing human gut microbiota is incredibly important. One step in characterizing gut microbiota is evaluating the influences of different diets on the diversity of microorganisms and species abundance in the human gut. De Filippis et al. examined the influence of western omnivorous diets, vegetarian, and vegan diets on bacterial strain abundance using shotgun metagenomics (De Filippis et al., 2019). In this project the metagenome data from the vegan and omnivorous diets of this paper will be examined.
  I will be using fastqc for quality control of the reads as it is an efficient method to evaluate fastq files. I will be using kraken2 and braken for taxonomic classification of the reads as this process is efficient, in spite of the possibilitiy of a high false positive rate. Using RStudio I will also compute richness and dominance using Chao1, Berger-Parker, and the Shannon index to thoroughly evaluate the diversity within each of the samples. To measure beta diversity I will be performing a Bray-Curtis analysis to determine how similar the samples from the vegan and omnivorous groups are. Aldex3 will be used to determine differential abundance of species between the groups as it is only relatively conservative, which should work better for the small sample size used in this project.
- n=3, vegan and omnivore factors
- shotgun genomics: needed for genome reconstruction, using full reference genomes, need sufficient database & sample size
- gut biome samples
- short-read illumina samples
- use core_nt database because it also evaluates eukaryotic genomes
# Methods  
 To do this, fastqc will be used to determine the quality of the raw reads provided from the paper. Then, kraken2 and bracken will be used to provide and refine the taxonomic abundance of the reads using k-mer-based matching (https://github.com/nmoragas/metagenomics_pipeline). From there the files will be brought into an RStudio environment where taxonomic abundance within the groups will be compared alongside evaluation of alpha and beta metrics. Finally, aldex3 will be used to run a differential abundance analysis on the data to determine differences in  
- human reads already removed, qc done using fastqc
- kraken2 for taxonomy indexing
- bracken for measuring taxa abundance
- Rstudio packages for alpha and beta diversity metrics
- aldex3 for differential abundance analysis between vegan and omnivore
</ins>Figure 1.<\ins> <img width="3000" height="1500" alt="taxonomic_abundance_by_group" src="https://github.com/user-attachments/assets/ddefd97f-b37b-49c2-9a5b-0e683da04c82" />
<img width="2400" height="1500" alt="shannon_barplot" src="https://github.com/user-attachments/assets/e046241c-07ab-4a92-a7bb-19eda740b1b0" />
<img width="2400" height="1500" alt="chao1_barplot" src="https://github.com/user-attachments/assets/9bcdf5b6-b602-4dd0-b0aa-596bb9eceb29" />
<img width="2100" height="1500" alt="bray_curtis_pcoa" src="https://github.com/user-attachments/assets/f306db10-2b8f-43ba-b99b-feb740080eec" />
<img width="2400" height="1500" alt="berger_parker_barplot" src="https://github.com/user-attachments/assets/929f6611-2667-44a7-a770-2c0b59220ceb" />
<img width="2700" height="1650" alt="aldex3_top_taxa" src="https://github.com/user-attachments/assets/f947e7ae-c7cc-4ff5-a808-c7a98ffd6280" />

# Results  

# Discussion
# References
