dir0="../data/"
for fn in ../data/*R1_001.fastq.gz;
do
    fn_out=${fn##*/}
    baseName=${fn%%_R1_001.fastq.gz*}
    base0=${baseName##*/}
    echo $base0 $fn
    r1="${dir0}${base0}_R1_001.fastq.gz"
    r2="${dir0}${base0}_R2_001.fastq.gz"
    ref="../assets/wg_reference.fasta"
    if [ -f "../sequences/${base0}.fa" ]; then
        echo "skipping ${base0}"
    else
        echo $base0
        samtools mpileup -d 1000 -A -Q 0 "../trimmed/${base0}.sorted.bam"  | ivar consensus -p "../sequences/${base0}" -q 20 -t 0.9
    fi
done