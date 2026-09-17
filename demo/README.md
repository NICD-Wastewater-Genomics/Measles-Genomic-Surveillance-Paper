## Demo

This is small, self-contained demo indended to demonstrate the bioinformatic pipeline used to generate the wastewater sequencing results in this paper. The demo downloads a handful of raw FASTQs directly from NCBI SRA and runs them through the same mapping, primer trimming, consensus-calling, and Freyja variant-calling/demixing steps as `../src/run_samples.sh` and `../src/run_consensus.sh`, using the reference, primer scheme, and lineage definitions in `../assets`.

This demo only processes 3 samples, and thus is not meant to reproduce paper figures. It is intended as a quick way to test out the approach on your system. 

### Pipeline steps
For each sample:
1. Download paired-end reads from SRA (`prefetch` + `fasterq-dump`)
2. Map reads to the measles whole-genome reference with `minimap2`
3. Trim amplicon primers with `ivar trim`
4. Generate a consensus sequence with `samtools mpileup` + `ivar consensus`
5. Call variants and sequencing depth with `freyja variants`
6. Estimate lineage/genotype relative abundance with `freyja demix`

Outputs from all samples are then combined using `freyja aggregate` into `agg_demixed.tsv`.

### Samples
Three small paired-end amplicon WGS runs from wastewater samples in BioProject [PRJNA1377662](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1377662) (see `accessions.txt`): `SRR37574359`, `SRR37574380`, `SRR37574398`. Each is only a few MB, so the whole demo should run in a few minutes.

### Requirements
- `sra-tools` (`prefetch`, `fasterq-dump`)
- `minimap2`
- `samtools`
- `ivar`
- `freyja` >= 2.0 — earlier releases don't support the `--pathogen` and `--autoadapt` options used here. See the [Freyja repo](https://github.com/andersen-lab/Freyja) for installation instructions, compatible Python versions, and  details on the demixing method. Freyja has been widely validated on Linux and Mac operating systems. 

The easiest way to get all of these is via the included conda environment file:

```bash
conda env create -f environment.yml
conda activate measles-demo
```
The environment solving process usually completes within 1-2 minutes, and the download takes ~30 seconds provided a reliable internet connection. 

### Usage
From this directory:

```bash
bash run_demo.sh
```

The script will skip any sample whose reads are already downloaded or whose demixed output already exists, and can be re-run if the run is interrupted. 

### Outputs
```
demo/
├── data/       # downloaded FASTQs
├── bams/       # raw mapped BAMs
├── trimmed/    # primer-trimmed, sorted BAMs
├── sequences/  # per-sample consensus sequences (ivar consensus)
├── variants/   # per-sample SNV calls (freyja variants)
├── depths/     # per-sample depth-of-coverage
├── outputs/    # per-sample demixed lineage abundances
└── agg_demixed.tsv   # aggregated demixing results across all samples
```

### Next steps (although not run by this demo)
1. Pull background sequences for phylogenetic context. `../tree/fetch_sequences.py` retrieves genotype-diverse measles genomes from GenBank by taxid, filtering by length/base-composition and (optionally) human host, as used in `../tree/background_tree.sh`:

   ```bash
   python ../tree/fetch_sequences.py --taxid 11234 --min_length_fraction 0.8 \
       --min_gatc_fraction 0.8 --output_prefix measles_bg --human_only
   ```

2. Combine the demo consensus sequences with those background sequences and align them with `mafft`, the same way whole-genome consensus sequences are aligned ahead of phylogenetic inference in `../tree/subtrees_withprivate.sh`:

   ```bash
   cat sequences/*.fa measles_bg.fasta > demo_consensus_plus_bg.fasta
   mafft --maxiterate 1000 --thread 4 demo_consensus_plus_bg.fasta > demo_consensus_aligned.fasta
   ```

This aligned FASTA is the starting point for phylogenetic tree inference (e.g. with `iqtree`/`treetime`), as done for the full dataset in `../tree`.
