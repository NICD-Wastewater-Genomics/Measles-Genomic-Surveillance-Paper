dir0="../data/"
for fn in ../data/*R1_001.fastq.gz;
do
    fn_out=${fn##*/}
    baseName=${fn%%_R1_001.fastq.gz*}
    base0=${baseName##*/}

    if [ -f "../outputs/${base0}.demixed.tsv" ]; then
    echo "Already processed ${fn}"
    else
    echo "Running ${fn}"
    r1="${dir0}${base0}_R1_001.fastq.gz"
    r2="${dir0}${base0}_R2_001.fastq.gz"
    ref="../assets/wg_reference.fasta"
    bed="../assets/measles_primers.bed"
    echo $r1 $r2
    minimap2 -ax sr -O16,36 $ref $r1 $r2 | samtools view -bS - >"../bams/${base0}.bam" #increase gap opening param, since we're seeing some artifactual 1nt indels
    samtools sort -o "../bams/${base0}.sorted.bam" "../bams/${base0}.bam" && samtools index "../bams/${base0}.sorted.bam"
    ivar trim -x 3 -e -m 80 -i "../bams/${base0}.sorted.bam" -b $bed -p "../trimmed/${base0}.bam"
    samtools sort -o "../trimmed/${base0}.sorted.bam" "../trimmed/${base0}.bam"
    samtools index "../trimmed/${base0}.sorted.bam"

    freyja variants "../trimmed/${base0}.sorted.bam" --variants "../variants/${base0}.tsv" --depths "../depths/${base0}.tsv" --ref $ref
    freyja demix "../variants/${base0}.tsv" "../depths/${base0}.tsv" --output "../outputs/${base0}.demixed.tsv" --depthcutoff 5 --autoadapt --pathogen MEASLES --lineageyml "../assets/lineagesMEASLES.yml"
    fi
done

freyja aggregate ../outputs/ --output ../agg_demixed.tsv --ext tsv
