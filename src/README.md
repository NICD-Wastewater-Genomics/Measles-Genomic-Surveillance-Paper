## src

R and Python scripts used for bioinformatic processing and epidemiological/figure analyses.

### Bioinformatic bash scripts
- `run_samples.sh`: Runs the Freyja demixing pipeline on raw wastewater FASTQs
- `run_consensus.sh`: Generates consensus sequences from trimmed BAMs using `samtools mpileup` + `ivar consensus`


### Figures and analyses
- `plot_clinical_map_district_level.py`: District-level map of clinical case data (Fig 1a)
- `plot_ww_pos_map_district_level.py`: District-level map of wastewater positivity (Fig 1b, ED Fig 3)
- `wastewater_vs_clin_pos_rate.R`: Compares wastewater and clinical positivity rates over time (Fig 1d, ED Figs 2a-d)
- `concentration_and_coverage.py`: Wastewater concentration vs. sequencing coverage analyses (ED Fig 1b, ED Fig 5c)
- `plot_variability.sh`: Aligns whole-genome sequences to the reference and calls variable sites (mafft + snp-sites), then calls `highlight_sa_snps.py`

- `calc_diversity.py`: Within-sample diversity calculations (ED Fig 7b,c)
- `highlight_sa_snps.py`: Highlights South African-specific SNPs from VCF output (ED Figs 7a, 7d)


- `clincal_genotype_counts.R`: Sequencing outcomes and genotype counts (D8/B3) from clinical data (Fig 2a)
- `dynamics_over_time.py`: Case/lineage dynamics over time (Fig 2b,c)
- `plot_ww_genotype_on_map.py`: Maps wastewater genotype detections (Fig 2d,e; Fig 3b,c)
- `plot_clinical_genotype_district.py`: District-level clinical genotype summaries (Fig 3a)

- `plot_site_variation.py`: Visualization of non-dominant SNVs present in wastewater samples for both genotype B3 and D8 (ED Fig. 9)

### `simulations/`
Scripts and inputs for read-simulation benchmarking of the Freyja mixture/demixing workflow (ED Fig. 6):
#### Mixture analysis
- `Snakefile-simulate`: Simulates amplicon reads per isolate (via `bygul`) for pure samples
- `snakefile-freyja`: Aligns simulated reads and runs Freyja demixing
- `summarise_freyja_mixtures.R`: Aggregates and summarizes Freyja demixing results across simulated mixtures

#### Amplicon dropout analysis
- `check_primers.sh`: Amplicon dropout analysis full pipeline, calls below helper scripts
- `extract_seqs.py`: Extract amplicons from sequences in MSA
- `pull_fasta.py`: pulls the fasta sequences from NCBI
- `create_mixtures.py`: Builds simulated mixtures of reference isolates
- `amplicon_dropouts.R`: Create dropout heatmaps
- `assemblies/`: Reference isolate assemblies used to build simulated mixtures
