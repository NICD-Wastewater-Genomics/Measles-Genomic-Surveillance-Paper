## Integrating measles wastewater and clinical whole-genome sequencing enables high-resolution tracking of virus evolution and transmission
Abstract
---
Measles outbreaks have surged globally in recent years, but current surveillance systems have limited capacity to monitor measles virus (MeV) transmission and evolution at population scale. Although MeV can be detected in wastewater, the public health potential of wastewater genomic surveillance for MeV remains largely unexplored. Here, we deploy sensitive, low-cost MeV wastewater genomic surveillance combining virus concentration, whole-genome amplicon sequencing, and bioinformatic analysis alongside routine clinical genomic surveillance during the 2024-25 outbreak in South Africa. Integrated phylogenetic analyses of wastewater and clinical MeV genomes revealed previously undetected interprovincial spread and transmission links not captured by standard N450 sequencing. Our findings demonstrate that wastewater-integrated whole-genome surveillance expands the coverage and resolution of routine MeV monitoring and provides a scalable tool to advance measles control and elimination efforts.  

---
This repository contains the following folders: 

- `src`: R and Python scripts needed for all basic bioinformatic analyses
- `variants`: single nucleotide variant (SNV) frequencies for each wastewater and clinical sample
- `depths`: coverage depth, by nucleotide position relative to Hu-1 reference
- `tree`: Additional scripts used for phylogenetic analyses 
- `n450_tree`: Additional scripts used for N450-only phylogenetic analyses

The `agg_demixed.tsv` file contains aggregated Freyja outputs for sequenced wastewater samples. 

All wastewater sequencing raw data is available via NCBI SRA (PRJNA1377662), and consensus sequences are available via Pathoplexus (https://pathoplexus.org/seqsets/PP_SS_892). 