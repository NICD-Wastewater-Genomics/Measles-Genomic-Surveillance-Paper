## tree

Scripts used to build and analyze the whole-genome phylogenetic trees (background context sequences + South African wastewater/clinical sequences).

### Building trees
- `background_tree.sh`: Downloads background (GenBank) sequences for the relevant taxon, determines genotypes for unlabeled sequences, builds a VCF for UShER, and plots the resulting tree by genotype
- `fetch_sequences.py`: Pulls all sequences matching a query from GenBank; used by `background_tree.sh`
- `get_clades_genbank.py`: Extracts clade/genotype assignments from metadata where possible, correcting incorrect GenBank entries
- `add_dates_to_fasta.py`: Adds collection date information to FASTA headers
- `label_new_consensus.py`: Adds date information to new South African consensus sequences
- `subtrees_withprivate.sh`: Core script for whole-genome sequence alignment and tree inference, including private clinic-collected sequences

### Mutation-annotated trees
- `mutation_annotated_subtrees.sh`: Converts trees to UShER mutation-annotated tree (MAT) format and extracts sequence-specific paths for South African sequences
- `get_key_mutations.py`: Extracts key mutations from mutation-annotated trees

### Plotting and clade extraction
- `plot_subtrees_withprivate.py`: Plots whole-genome subtrees with private/South African sequences highlighted, including Mantel test comparisons between distance matrices
- `tree_plot_and_clade_extract.py`: Extracts genotypes for background sequences not already present in metadata

### N450 vs. whole-genome comparisons
- `prune_trees.R`: Trims trees to tips present in both N450 and whole-genome datasets, for PS/AU test comparisons
- `clustering_metrics.R`: Calculates parsimony scores comparing whole-genome and N450 sequencing approaches
- `convert_nexus.R`: Converts Nexus files to Newick format, for use with IQ-TREE AU testing

### Supporting data
- `measles_bg_metadata_cleaned.csv`: Cleaned metadata for background sequences
- `measles_tree_bg.pdf`: Rendered background tree figure
