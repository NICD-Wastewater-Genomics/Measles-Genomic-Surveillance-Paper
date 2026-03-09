## core script for WG sequence alignment and tree inference
for gt in B3 D8;
do
echo "working on genotype ${gt}"
# NOTE: bg metadata has some incorrect genotypes assignments that were manually fixed. 
python add_dates_to_fasta.py measles_bg_metadata_cleaned.csv "${gt}"/measles_"${gt}"_bg.fasta "${gt}"/measles_"${gt}"_bg_dated.fasta collection_date_iso seqs_to_filter.txt

cat  ../sequences/seqs_ww_MEASLESgenome-"${gt}"/* > "${gt}"_ww_seqs.fasta
if [ -f "${gt}/${gt}_timetree/outliers.tsv" ]; then
    python label_new_consensus.py ../metadata/Measles_seqdata_02022026.csv "${gt}"_ww_seqs.fasta "${gt}"_ww_seqs_dated.fasta SampleID SampleCollectionDate "${gt}/${gt}_timetree_withprivate/outliers.tsv"
else
    python label_new_consensus.py ../metadata/Measles_seqdata_02022026.csv "${gt}"_ww_seqs.fasta "${gt}"_ww_seqs_dated.fasta SampleID SampleCollectionDate
fi
cat  ../sequences/seqs_clinical_MEASLESgenome-"${gt}"/* > "${gt}"_clinical_seqs.fasta
python label_new_consensus.py ../metadata/mev\ clinical\ metadata_updated\ 040226\ ks_JL_090226.csv "${gt}"_clinical_seqs.fasta "${gt}"_clinical_seqs_dated_withprivate.fasta Seq\ ID Onset\ date

# run aligment and fit ML tree with iqtree
cat "${gt}"/measles_"${gt}"_bg_dated.fasta "${gt}"_ww_seqs_dated.fasta "${gt}"_clinical_seqs_dated_withprivate.fasta > "${gt}"/"${gt}"_bg_plus_nicd_withprivate.fasta
mafft --maxiterate 1000 --thread 20 "${gt}"/"${gt}"_bg_plus_nicd_withprivate.fasta > "${gt}"/"${gt}"_wg_all_dated_aligned_withprivate.fasta
iqtree -s "${gt}"/"${gt}"_wg_all_dated_aligned_withprivate.fasta -T 20 -m GTR -redo 

head -n1 "${gt}"/measles_"${gt}"_bg_dated_meta.tsv > "${gt}"/"${gt}"_all_meta.tsv  # Get the header from the first file
#concatenating all metadata, dropping header after the first one
for f in "${gt}"/measles_"${gt}"_bg_dated_meta.tsv "${gt}"_ww_seqs_dated_meta.tsv "${gt}"_clinical_seqs_dated_withprivate_meta.tsv; do tail -n+2 "$f" >> "${gt}"/"${gt}"_all_meta.tsv; done 
# drop sequences exceeding the clock filter
# infer time-calibrated tree

treetime --aln "${gt}"/"${gt}"_wg_all_dated_aligned_withprivate.fasta --dates "${gt}"/"${gt}"_all_meta.tsv --tree "${gt}"/"${gt}"_wg_all_dated_aligned_withprivate.fasta.treefile --outdir "${gt}"/"${gt}"_timetree_withprivate/
done
