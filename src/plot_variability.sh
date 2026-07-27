for gt in B3 D8;
do
    # mafft --auto --thread 5 --keeplength --addfull \
    #     "../tree/${gt}/${gt}_wg_all_dated_aligned_withprivate.fasta" \
    #     ../assets/wg_reference.fasta \
    #     > "../ref_alignments/${gt}_wg_aligned_to_ref.fasta"

    # snp-sites -v -o "../ref_alignments/${gt}_variable_sites.vcf" "../ref_alignments/${gt}_wg_aligned_to_ref.fasta"

    python highlight_sa_snps.py "../ref_alignments/${gt}_variable_sites.vcf" "${gt}"

done