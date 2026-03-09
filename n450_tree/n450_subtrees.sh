for gt in B3 D8;
do
    echo "working on genotype ${gt}"
    # pull from full tree, including sequences from private clinics
    fn="../tree/${gt}"/"${gt}"_bg_plus_nicd_withprivate.fasta
    mkdir -p "${gt}"
    # align and trim to N450 coordinates. 
    minimap2 -a ../assets/wg_reference.fasta $fn | gofasta sam toMultiAlign --trim --trimstart 1233 --trimend 1682 --pad > "${gt}"/"${gt}"_aligned.fasta
    # remove ref seq, since it's not B3 or D8
    awk '/^>/ {f=($1!=">NC_001498.1")} f' "${gt}"/"${gt}"_aligned.fasta > "${gt}"/"${gt}"_aligned_wo_ref.fasta

    # remove sequences without the N450 present
    python clean_seqs.py "${gt}"/"${gt}"_aligned_wo_ref.fasta "${gt}"/"${gt}"_aligned_wo_ref_clean.fasta

    iqtree -s "${gt}"/"${gt}"_aligned_wo_ref_clean.fasta -T 10 -m GTR -redo

    treetime --aln "${gt}"/"${gt}"_aligned_wo_ref_clean.fasta --dates "../tree/${gt}"/"${gt}"_all_meta.tsv --tree "${gt}"/"${gt}"_aligned_wo_ref_clean.fasta.treefile --outdir "${gt}"/"${gt}"_timetree/

done
