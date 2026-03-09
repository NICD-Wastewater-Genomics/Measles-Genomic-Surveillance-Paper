# script to convert trees to UShER format and
# extract sequence-specific paths for SA sequences

for gt in B3 D8;
do
echo "working on genotype ${gt}"
all_seqs="${gt}"/"${gt}"_bg_plus_nicd_withprivate.fasta
mafft_aligned="${gt}"/"${gt}"_wg_all_dated_aligned_withprivate.fasta
mafft_tree="${gt}"/"${gt}"_wg_all_dated_aligned_withprivate.fasta.treefile
all_meta="${gt}"/"${gt}"_all_meta.tsv
cat $all_seqs ../assets/wg_reference.fasta > "${gt}"/"${gt}"_all_wref.fasta
minimap2 -a ../assets/wg_reference.fasta "${gt}"/"${gt}"_all_wref.fasta | gofasta sam toMultiAlign --pad > "${gt}"/"${gt}"_all_ref_aligned.fasta
faToVcf "${gt}"/"${gt}"_all_ref_aligned.fasta "${gt}"/"${gt}"_all_ref_aligned.vcf -ref=NC_001498.1

usher -t $mafft_tree -v "${gt}"/"${gt}"_all_ref_aligned.vcf -o "${gt}"/"${gt}"_all_ref_aligned.pb --threads 4
# for viz, can extract to taxonium format
# usher_to_taxonium --input "${gt}"/"${gt}"_all_ref_aligned.pb --genbank ../assets/measles_ref.gb --output "${gt}"/"${gt}"_all_ref_aligned.jsonl.gz

# use matutils to extract path to sequence tip
matUtils extract -i "${gt}"/"${gt}"_all_ref_aligned.pb  -S "${gt}"/"${gt}"_sequence_paths.txt
python get_key_mutations.py --input "${gt}"/"${gt}"_sequence_paths.txt --output "${gt}"/"${gt}"_mutation_counts.csv
done
