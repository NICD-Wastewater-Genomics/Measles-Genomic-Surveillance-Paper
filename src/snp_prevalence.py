### This script reports how commonly specific SNPs (e.g. A746G, T3348C) are observed
### across the sequences in a snp-sites VCF, overall and broken down by group.
import argparse
import os
import re
import sys

import numpy as np
import pandas as pd

from highlight_sa_snps import GENE_REGIONS, _get_gene_idx, gt_to_float, parse_date, parse_vcf

VALID_BASES = set('ACGT')
WW_PREFIXES = ('ENV', 'CST', 'NAT', 'ART_MEV')

# e.g. A746G (ref+alt), 746G (alt only), A746 (ref only), 746 (any alt at that site)
SNP_RE = re.compile(r'^([ACGT]?)(\d+)([ACGT]?)$', re.IGNORECASE)


def parse_snp(spec):
    """Parse a mutation string into (ref_or_None, pos, alt_or_None)."""
    m = SNP_RE.match(spec.strip())
    if not m:
        raise ValueError(f"Cannot parse SNP '{spec}' (expected e.g. A746G, 746G or 746)")
    ref, pos, alt = m.group(1).upper(), int(m.group(2)), m.group(3).upper()
    return (ref or None), pos, (alt or None)


def gene_of(pos):
    j = _get_gene_idx(pos)
    return GENE_REGIONS[j][0] if j < len(GENE_REGIONS) else 'intergenic'


def sample_metadata(sample_cols):
    """Split 'name|country|date' headers into a per-sample metadata frame."""
    recs = []
    for c in sample_cols:
        parts = c.split('|')
        name    = parts[0]
        country = parts[1] if len(parts) > 1 else 'unknown'
        date    = parts[2] if len(parts) > 2 else 'unknown'
        recs.append({
            'sample':  c,
            'name':    name,
            'country': country,
            'date':    date,
            'year':    date.split('-')[0] if date != 'unknown' else 'unknown',
            'region':  'SouthAfrica' if country == 'SouthAfrica' else 'Other',
            'type':    'wastewater' if name.startswith(WW_PREFIXES) else 'clinical',
        })
    return pd.DataFrame(recs)


def expand_sites(df, sample_cols):
    """One row per (POS, ALT) allele with a 1/0/NaN call per sample.

    1 = carries this alt, 0 = carries some other clean allele, NaN = uncalled or
    resolved to an ambiguity/gap allele (treated as missing, as in highlight_sa_snps).
    """
    gt_mat = np.column_stack([df[c].apply(gt_to_float).values for c in sample_cols])

    positions, refs, alts, rows = [], [], [], []
    for i, (pos, ref, alt_str) in enumerate(zip(df['POS'].astype(int), df['REF'], df['ALT'])):
        gt_row   = gt_mat[i]
        site_alts = alt_str.split(',')
        ambig_gt = {j + 1 for j, a in enumerate(site_alts) if a.upper() not in VALID_BASES}
        is_nan   = np.isnan(gt_row)
        is_ambig = np.isin(gt_row, list(ambig_gt)) & ~is_nan
        clean    = ~is_nan & ~is_ambig
        for alt_idx, alt in enumerate(site_alts, start=1):
            if alt.upper() not in VALID_BASES:
                continue
            row = np.full(len(gt_row), np.nan)
            row[clean] = (gt_row[clean] == alt_idx).astype(float)
            positions.append(pos)
            refs.append(ref.upper())
            alts.append(alt.upper())
            rows.append(row)

    return (np.array(positions), np.array(refs), np.array(alts),
            np.vstack(rows) if rows else np.empty((0, len(sample_cols))))


def match_snp(spec, positions, refs, alts):
    """Indices of expanded sites matching a parsed SNP spec."""
    ref, pos, alt = spec
    hit = positions == pos
    if ref is not None:
        hit &= refs == ref
    if alt is not None:
        hit &= alts == alt
    return np.where(hit)[0]


def summarise(calls, meta):
    """Overall counts for one allele across the samples in `meta`."""
    called  = ~np.isnan(calls)
    carrier = called & (calls == 1)
    n_called, n_alt = int(called.sum()), int(carrier.sum())
    dates = [d for d in meta.loc[carrier, 'date'] if 'X' not in d.upper() and d != 'unknown']
    return {
        'n_alt':       n_alt,
        'n_called':    n_called,
        'n_missing':   int((~called).sum()),
        'AF':          round(n_alt / n_called, 4) if n_called else np.nan,
        'n_countries': int(meta.loc[carrier, 'country'].nunique()),
        'first_date':  min(dates) if dates else '',
        'last_date':   max(dates) if dates else '',
    }


def breakdown(calls, meta, key):
    """Per-group counts for one allele, one row per level of `meta[key]`."""
    called  = ~np.isnan(calls)
    carrier = called & (calls == 1)
    out = []
    for level, idx in meta.groupby(key, sort=True).groups.items():
        sel   = meta.index.isin(idx)
        n_c   = int((called & sel).sum())
        n_a   = int((carrier & sel).sum())
        out.append({
            'group':    key,
            'level':    level,
            'n_alt':    n_a,
            'n_called': n_c,
            'AF':       round(n_a / n_c, 4) if n_c else np.nan,
        })
    return out


def main():
    ap = argparse.ArgumentParser(
        description='Report how commonly specific SNPs are observed in a snp-sites VCF.')
    ap.add_argument('vcf', help='VCF from snp-sites, e.g. ../ref_alignments/D8_variable_sites.vcf')
    ap.add_argument('snps', nargs='*',
                    help='mutations to query, e.g. A746G T3348C (also 746G or 746)')
    ap.add_argument('--top', type=int, default=0, metavar='N',
                    help='additionally report the N most commonly observed alt alleles')
    ap.add_argument('--group-by', nargs='*', default=['region', 'type'],
                    choices=['region', 'type', 'year', 'country'],
                    help='sample groupings to break each SNP down by (default: region type)')
    ap.add_argument('--min-called', type=int, default=1, metavar='N',
                    help='drop breakdown rows with fewer than N called samples (default: 1)')
    ap.add_argument('--keep-reference', action='store_true',
                    help='keep the reference column (header without |country|date) as a sample')
    ap.add_argument('--out-prefix', default=None,
                    help='write <prefix>_snp_prevalence_by_group.tsv '
                         '(default: alongside the input VCF, e.g. ../ref_alignments/D8)')
    ap.add_argument('--no-save', action='store_true', help='print results without writing the TSV')
    args = ap.parse_args()

    if not args.snps and not args.top:
        ap.error('give at least one SNP, or use --top N')

    try:
        specs = [parse_snp(s) for s in args.snps]
    except ValueError as e:
        sys.exit(str(e))

    # Default to writing beside the input VCF, stripping the snp-sites suffix
    if args.out_prefix is None:
        stem = os.path.basename(args.vcf)
        for suffix in ('_variable_sites.vcf', '.vcf'):
            if stem.endswith(suffix):
                stem = stem[:-len(suffix)]
                break
        args.out_prefix = os.path.join(os.path.dirname(args.vcf) or '.', stem)

    df = parse_vcf(args.vcf)
    sample_cols = df.columns[9:].tolist()
    if not args.keep_reference:
        sample_cols = [c for c in sample_cols if '|' in c]

    meta = sample_metadata(sample_cols)
    positions, refs, alts, bin_mat = expand_sites(df, sample_cols)
    print(f"{args.vcf}: {len(sample_cols)} samples, {len(df)} variable sites, "
          f"{len(positions)} biallelic-expanded alleles", file=sys.stderr)

    # Resolve requested SNPs to expanded-site indices, keeping the query order
    queries = []
    for spec, raw in zip(specs, args.snps):
        idx = match_snp(spec, positions, refs, alts)
        if len(idx) == 0:
            print(f"  {raw}: not observed in any sample (no matching ALT allele at this site)",
                  file=sys.stderr)
            continue
        queries.extend((raw, i) for i in idx)

    if args.top:
        af = np.array([np.nanmean(bin_mat[i]) if np.any(~np.isnan(bin_mat[i])) else np.nan
                       for i in range(len(positions))])
        order = np.argsort(-np.nan_to_num(af, nan=-1))[:args.top]
        seen  = {i for _, i in queries}
        queries.extend((f"{refs[i]}{positions[i]}{alts[i]}", i)
                       for i in order if i not in seen)

    if not queries:
        sys.exit('No requested SNP was found in the VCF.')

    summary_rows, group_rows = [], []
    for label, i in queries:
        calls = bin_mat[i]
        mut   = f"{refs[i]}{positions[i]}{alts[i]}"
        row   = {'query': label, 'mutation': mut, 'POS': int(positions[i]),
                 'REF': refs[i], 'ALT': alts[i], 'gene': gene_of(int(positions[i]))}
        row.update(summarise(calls, meta))
        summary_rows.append(row)
        for key in args.group_by:
            for g in breakdown(calls, meta, key):
                if g['n_called'] >= args.min_called:
                    group_rows.append({'mutation': mut, **g})

    summary = pd.DataFrame(summary_rows)
    groups  = pd.DataFrame(group_rows)

    pd.set_option('display.width', 200, 'display.max_rows', None)
    print('\n=== Overall prevalence ===')
    print(summary.drop(columns='query').to_string(index=False))
    for mut in summary['mutation']:
        sub = groups[groups['mutation'] == mut]
        if sub.empty:
            continue
        print(f"\n=== {mut} by group ===")
        print(sub.drop(columns='mutation').to_string(index=False))

    if not args.no_save:
        out_tsv = f"{args.out_prefix}_snp_prevalence_by_group.tsv"
        groups.to_csv(out_tsv, sep='\t', index=False)
        print(f"\nWrote {out_tsv}")


if __name__ == '__main__':
    main()
