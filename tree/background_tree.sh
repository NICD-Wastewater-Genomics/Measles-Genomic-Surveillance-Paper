## script for downloading all background sequences, and determing genotypes of sequences if not labeled

# download seqs from taxon, human only
python fetch_sequences.py --taxid 11234 --min_length_fraction 0.8 --min_gatc_fraction 0.8 --output_prefix measles_bg --human_only

minimap2 -a ../assets/wg_reference.fasta measles_bg.fasta | gofasta sam toMultiAlign --pad > background/bg_aligned.fasta

# create vcf for usher. 
faToVcf background/bg_aligned.fasta background/bg_aligned.vcf -ref=NC_001498.1

# ##remove reference seq. 
awk 'BEGIN{while((getline<"seqs_to_filter.txt")>0)l[">"$1]=1}/^>/{f=!l[$1]}f' background/bg_aligned.fasta > background/bg_aligned_wo_ref.fasta
iqtree -s background/bg_aligned_wo_ref.fasta -T 10 -m GTR -redo 

# plots tree and separates out sequences by genotype
python tree_plot_and_clade_extract.py
