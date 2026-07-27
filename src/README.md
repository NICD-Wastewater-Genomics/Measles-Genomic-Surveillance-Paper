## src

R and Python scripts used for bioinformatic processing and epidemiological/figure analyses.

### Bioinformatic bash scripts
- `run_samples.sh`: Runs the Freyja demixing pipeline on raw wastewater FASTQs
- `run_consensus.sh`: Generates consensus sequences from trimmed BAMs using `samtools mpileup` + `ivar consensus`


### Figures and analyses
- `plot_variability.sh`: Aligns whole-genome sequences to the reference and calls variable sites (mafft + snp-sites), then calls `highlight_sa_snps.py`
- `highlight_sa_snps.py`: Highlights South African-specific SNPs from VCF output (ED Figs 7a, 7d)

- `concentration_and_coverage.py`: Wastewater concentration vs. sequencing coverage analyses (SFig 1b, SFig 3b)
- `calc_diversity.py`: Within-sample diversity calculations (SFig 5a,b)
- `wastewater_vs_clin_pos_rate.R`: Compares wastewater and clinical positivity rates over time (Fig 1d, SFigs 2a-d)
- `clincal_genotype_counts.R`: Sequencing outcomes and genotype counts (D8/B3) from clinical data (Fig 2a)
- `dynamics_over_time.py`: Case/lineage dynamics over time
- `plot_ww_genotype_on_map.py`: Maps wastewater genotype detections (Fig 2d,e; Fig 3b,c)
- `plot_ww_pos_map_district_level.py`: District-level map of wastewater positivity
- `plot_clinical_genotype_district.py`: District-level clinical genotype summaries
- `plot_clinical_map_district_level.py`: District-level map of clinical case data

### `simulations/`
Scripts and inputs for read-simulation benchmarking of the Freyja mixture/demixing workflow:
- `create_mixtures.py`: Builds simulated mixtures of reference isolates
- `Snakefile-simulate`: Simulates amplicon reads per isolate (via `bygul`)
- `snakefile-freyja`: Aligns simulated reads and runs Freyja demixing
- `summarise_freyja_mixtures.R`: Aggregates and summarizes Freyja demixing results across simulated mixtures
- `aggregated_result.tsv`, `lineages.tsv`, `rep_isolates.txt`, `reference.fasta`, `primer_v3_400.bed`: supporting inputs
- `assemblies/`: Reference isolate assemblies used to build simulated mixtures
