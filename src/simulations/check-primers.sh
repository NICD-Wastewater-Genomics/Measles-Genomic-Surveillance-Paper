# pull sequences
python pull_fasta.py nextclade_sequence_lineages.tsv all_seqs.fasta
# run bygul check primers
bygul check-primers all_seqs.fasta primer_v3_400.bed reference.fasta
# align requences
minimap2 -a -x asm20 --sam-hit-only --secondary=no --score-N=0 reference.fasta  all_seqs.fasta -o aligned.sam
gofasta sam toMultiAlign -s aligned.sam -o aligned.fasta
# extract primer sequences in the alignment
python extract_seqs.py
# creat heatmap of dropout ratios for each amplicon per clade
Rscript amplicon_dropouts.R