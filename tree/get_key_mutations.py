import argparse
from datetime import datetime
from freyja import convert_paths2barcodes
import pandas as pd
from functools import reduce
## script for mutation extraction
# to be used with mutation annotated trees (mutation_annotated_subtrees.sh)

def clean_mut_list(x):
    mutgroups = x.strip().split(' ')
    mutgroups = [m.split(':')[-1].strip() for m in mutgroups]
    return mutgroups

def reversion_checking(mut_list):
    # check if a reversion is present.
    flipPairs = [(d, d[-1] + d[1:len(d)-1]+d[0]) for d in mut_list
                 if (d[-1] + d[1:len(d)-1]+d[0]) in mut_list]
    flipPairs = [list(fp) for fp in list(set(flipPairs))]
    # subtract lower of two pair counts to get the lineage defining mutations
    if len(flipPairs)>0:
        for f in flipPairs:
            p0 = mut_list.count(f[0])
            p1 = mut_list.count(f[1])
            if p0==p1:
                mut_list = [ml for ml in mut_list if (ml !=f[0]) or (ml!= f[1])]
            elif p1>p0:
                mut_list = [ml for ml in mut_list if (ml !=f[0])]
            elif p0>p1:
                mut_list = [ml for ml in mut_list if (ml!= f[1])]
        return mut_list
        
    else:
        return mut_list

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", help="sequence paths",required=True)
    parser.add_argument("--specific_seqs", help="seq list",default=None)
    parser.add_argument("--output", help="mutations and counts",required=True)
    args = parser.parse_args()
    df = pd.read_csv(args.input, sep='\t',header=None)
    df.iloc[:,1] = df.iloc[:,1].apply(lambda x:clean_mut_list(x))
    df['all_muts'] = df.iloc[:,1].apply(lambda x:[x0 for x_ in x for x0 in x_.split(',')])
    df['all_muts'] = df['all_muts'].apply(lambda x: reversion_checking(x))
    df['Country'] = df.iloc[:,0].apply(lambda x:x.split('|')[1])

    df = df[df['Country']=='SouthAfrica']


    common = reduce(lambda a, b: set(a) & set(b), df["all_muts"])
    df["muts_filtered"] = df["all_muts"].apply(lambda x: list(set([m for m in x if m not in common])))
    from collections import Counter
    from itertools import chain
    if args.specific_seqs is not None:
        print('not yet implemented')
    else:
        counts = Counter(chain.from_iterable(df["muts_filtered"]))
        mutation_counts = (
            pd.DataFrame(counts.items(), columns=["mutation", "count"])
            .sort_values("count", ascending=False)
            .reset_index(drop=True))

        # Write to CSV
        mutation_counts["pos"] = mutation_counts["mutation"].str.extract(r"(\d+)").astype(int)
        mutation_counts = mutation_counts.sort_values("pos").reset_index(drop=True)
        mutation_counts.to_csv(args.output, index=False)

