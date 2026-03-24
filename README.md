# BINF-6110_Assignment-3
Comparison of gut microbiota between individuals maintaining vegan and omnivorous diets.
# Introduction  
  Gut microbiota diversity is not only influenced by how someone lives, but also incredibly influential to their wellbeing. For instance, the health of an individual's microbiome is an incredibly influential factor to the development of chronic diseases (Singh et al., 2017). For this reason, characterizing human gut microbiota is incredibly important. One step in characterizing gut microbiota is evaluating the influences of different diets on the diversity of microorganisms and species abundance in the human gut. De Filippis et al. examined the influence of western omnivorous diets, vegetarian, and vegan diets on bacterial strain abundance using shotgun metagenomics (De Filippis et al., 2019). In this project the metagenome data from the vegan and omnivorous diets of this paper will be examined.
  I will be using fastqc for quality control of the reads as it is an efficient method to evaluate FASTQ files. I will be using kraken2 and braken for taxonomic classification of the reads as this process is efficient, in spite of the possibilitiy of a high false positive rate (https://www.nature.com/articles/s41596-022-00738-y). Using RStudio I will also compute richness and dominance using Chao1, Berger-Parker, and the Shannon index to thoroughly evaluate the diversity within each of the samples (https://link.springer.com/protocol/10.1007/978-1-0716-5009-7_11). To measure beta diversity I will be performing a Bray-Curtis analysis to determine how similar the samples from the vegan and omnivorous groups are. Aldex3 will be used to determine differential abundance of species between the groups as it is only relatively conservative, which should work better for the small sample size used in this project (https://github.com/jsilve24/ALDEx3).
- n=3, vegan and omnivore factors
- shotgun genomics: needed for genome reconstruction, using full reference genomes, need sufficient database & sample size
- gut biome samples
- short-read illumina samples
- use core_nt database because it also evaluates eukaryotic genomes
# Methods  
 ## Raw Reads
   Raw reads were obtained from De Filippis et al. using SRA toolkit to access the FASTQ files (De Filippis et al., 2019). The reads were sequenced by an Illumina sequencer and the human reads were removed before the reads were published. Six sequences total were used (n = 3) for both vegan and omnivorous groups. 
  ## Quality Control  
  Quality of the reads was determined by FASTQC (FastQC v0.12.1). Once the reads were all determined to be of good quality, the reads were then prepared for taxonomic classification.
  ## Taxonomic Classification  
  Taxonomic classification was done in a Compute Canada Nibi cluster where all samples were run using an array. Kraken2 (kraken2 v2.1.6) and bracken (bracken v3.0) were used to perform the taxonomic classification, the abundance tables of which were then fed into RStudio for statistical analysis. The kraken2 standard database was used to perform the taxonomic classification of the reads (https://github.com/DerrickWood/kraken2/wiki/Manual#kraken-2-databases).
  ## Statistical Analysis  
  RStudio (RStudio 2026.01.1+403) was the environment used to perform statistical analysis on the samples alongside the packages dplyr (dplyr v1.2.0), tidyr (tidyr v1.3.2), tibble (tibble v3.3.1), and grid (grid v4.5.0). Alpha and beta diversity metrics were measured alongside differential abundance of species between the samples using the package vegan (vegan v2.7.3). The chosen alpha diversity tests were Chao1 for species richness, Berger-Parker for evenness, and Shannon Index for a combined metric of both.
  Beta diversity was measured using a Bray-Curtis test. All tests were run on default perameters.
## Differential abundance of species
ALDEx3 (ALDEx3 v1.0.2) will be used to run a differential abundance analysis on the data to determine differences in species abundance between the groups. 
## Plots
GGplot2 (ggplot2 v4.0.2) will be used to form the plots displayed in this project.  \
# Results  
  Based on the analyses performed in this project, there is a great amount of taxonomic abundance and diversity between the groups, but significantly more so within the groups. 
<ins>Figure 1.</ins> Taxonomic abundance of the top 20 samples observed in both groups of samples. Based on this graph it can be observed that generally species abundance has a larger difference within the samples than across them.  \
<img width="3000" height="1500" alt="taxonomic_abundance_by_group" src="https://github.com/user-attachments/assets/ddefd97f-b37b-49c2-9a5b-0e683da04c82" />  \
<ins>Figure 2.</ins> Barplot of the Chao1 richness alpha diversity metric of each of the samples. In this graph it can be observed that there is a lot of variation in richness within the groups, however the richness of species within the vegan group could potentially have more consistent species richness across the samples.  \
<img width="2400" height="1500" alt="chao1_barplot" src="https://github.com/user-attachments/assets/9bcdf5b6-b602-4dd0-b0aa-596bb9eceb29" />  \
<ins>Figure 3.</ins> Barplot of the Berger-Parker evenness alpha diversity metric for each of the samples. Based on this plot there is not very much consistency for evenness within the groups as one sample from each group is significantly less even than the rest.  \
<img width="2400" height="1500" alt="berger_parker_barplot" src="https://github.com/user-attachments/assets/929f6611-2667-44a7-a770-2c0b59220ceb" />  \
<ins>Figure 3.</ins> A barplot visualizing the Shannon diversity index for each of the samples. As demonstrated, Shannon diversity varies heavily within both treatment groups, however their indeces are similar across the groups.  \
<img width="2400" height="1500" alt="shannon_barplot" src="https://github.com/user-attachments/assets/e046241c-07ab-4a92-a7bb-19eda740b1b0" />  \
<ins>Figure 4</ins> PCoA plot for the measured Bray-Curtis beta diversity analysis. This metric demonstrates that there are few similarities in relative abundance within the treatment groups and that there is generally more similarity between samples in different treatment groups.  \
<img width="2100" height="1500" alt="bray_curtis_pcoa" src="https://github.com/user-attachments/assets/f306db10-2b8f-43ba-b99b-feb740080eec" />  \
<ins>Figure 5.</ins> Visualisation of the top 25 differentially abundant species determined using Aldex3. Considering even the most differentially abundant species are still sometimes observed in both groups determines that there is likely a high level of similarity in terms of the taxa observed in each of the groups, especially of those found in the omnivorous group compared to the vegan group.
<img width="2700" height="1650" alt="aldex3_top_taxa" src="https://github.com/user-attachments/assets/f947e7ae-c7cc-4ff5-a808-c7a98ffd6280" />

# Discussion
# References
