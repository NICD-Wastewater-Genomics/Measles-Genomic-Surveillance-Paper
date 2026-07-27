## n450_tree

Scripts and outputs for N450-only phylogenetic analyses (trimmed to the N450 region used by standard measles genotyping), for comparison against the whole-genome trees in `tree/`.

### Building the trees
- `n450_subtrees.sh`: For each genotype (B3, D8), pulls the whole-genome tree sequences from `tree/`, aligns and trims to N450 coordinates with `minimap2`/`gofasta`, removes the reference sequence, drops sequences without full N450 coverage, builds a tree with IQ-TREE, and infers a time-resolved tree with TreeTime
- `clean_seqs.py`: helper script, removes sequences lacking N450 coverage (i.e. regions of all Ns) from an alignment

### Plotting
- `plot_subtrees.py`: Plots the N450 subtrees (by genotype), matched to clinical metadata (ED Fig 10)

### `B3/` and `D8/`
Per-genotype outputs from `n450_subtrees.sh`:
- `<GT>_aligned_wo_ref_clean.fasta.treefile`: IQ-TREE maximum-likelihood tree
- `<GT>_time_resolved_tree.pdf`: Rendered time-resolved tree figure
- `<GT>_timetree/`: TreeTime output directory (`timetree.nexus`/`.pdf`, `divergence_tree.nexus`, `n450_<GT>_timetree.nwk` and pruned variant, `ancestral_sequences.fasta`, `branch_mutations.txt`, `dates.tsv`, `molecular_clock.txt`, `sequence_evolution_model.txt`, `root_to_tip_regression.pdf`, `auspice_tree.json`)
