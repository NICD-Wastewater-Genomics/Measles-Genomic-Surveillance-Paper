### this script generates new ED Figs 7a and 7d. 
import sys
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import numpy as np

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['savefig.bbox'] = 'tight'
matplotlib.rcParams['font.family'] = 'Arial'

GLOBAL_AF_MAX = 0.9  # exclude sites where alt is in >90% of called samples

GENE_REGIONS = [
    ('N',     56,  1744),
    ('P/V/C', 1748,  3402),
    ('M',    3406,  4872),
    ('F',    4876,  7247),
    ('H',    7251,  9208),
    ('L',    9212, 15854),
]
_GENE_COLORS    = ['#888888'] * 6
_INTERGENIC_CLR = '#888888'

def parse_date(s):
    try:
        return pd.to_datetime(s)
    except Exception:
        year = s.split('-')[0]
        return pd.to_datetime(f"{year}-07-01")


def _get_gene_idx(pos):
    for j, (_, start, end) in enumerate(GENE_REGIONS):
        if start <= pos <= end:
            return j
    return len(GENE_REGIONS)  # intergenic


def parse_vcf(path):
    header, rows = [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                header = line.lstrip('#').rstrip('\n').split('\t')
            else:
                rows.append(line.rstrip('\n').split('\t'))
    return pd.DataFrame(rows, columns=header)


def gt_to_float(gt):
    if gt in ('.', './.', '.|.'):
        return np.nan
    return float(int(gt.split('/')[0].split('|')[0]))

if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.exit(f"Usage: {sys.argv[0]} <input.vcf> <genotype_label>")

    vcf_path, gt_label = sys.argv[1], sys.argv[2]
    df = parse_vcf(vcf_path)

    sample_cols = df.columns[9:].tolist()
    _cutoff = pd.Timestamp('2024-01-01')
    sa_idx     = [i for i, c in enumerate(sample_cols)
                  if '|SouthAfrica|' in c
                  and parse_date(c.split('|')[2]) >= _cutoff]
    non_sa_idx = [i for i, c in enumerate(sample_cols) if '|SouthAfrica|' not in c]

    print(f"Total samples: {len(sample_cols)}  SA 2024 or after: {len(sa_idx)}  non-SA: {len(non_sa_idx)}")

    # Raw genotype matrix: shape (n_sites, n_samples)
    gt_mat = np.column_stack([df[c].apply(gt_to_float).values for c in sample_cols])

    # Explode multi-allelic sites: one row per alt allele
    valid_bases = set('ACGT')

    exp_positions, exp_refs, exp_alts, exp_rows = [], [], [], []
    for i, (pos, ref, alt_str) in enumerate(zip(
            df['POS'].astype(int), df['REF'], df['ALT'])):
        gt_row = gt_mat[i]
        alts = alt_str.split(',')
        ambig_gt = {j + 1 for j, a in enumerate(alts) if a.upper() not in valid_bases}
        for alt_idx, alt in enumerate(alts, start=1):
            row     = np.full(len(gt_row), np.nan)
            is_nan  = np.isnan(gt_row)
            is_ambig = np.isin(gt_row, list(ambig_gt)) & ~is_nan
            clean   = ~is_nan & ~is_ambig
            row[clean] = (gt_row[clean] == alt_idx).astype(float)
            exp_positions.append(pos)
            exp_refs.append(ref)
            exp_alts.append(alt)
            exp_rows.append(row)

    positions = np.array(exp_positions)
    refs      = np.array(exp_refs)
    alts_col  = np.array(exp_alts)
    bin_mat   = np.vstack(exp_rows)  # (n_expanded_sites, n_samples)

    sa_mat     = bin_mat[:, sa_idx]
    non_sa_mat = bin_mat[:, non_sa_idx]

    # Classify SA samples as wastewater or clinical before filtering
    WW_PREFIXES = ('ENV', 'CST', 'NAT', 'ART_MEV')
    sa_names = [sample_cols[i].split('|')[0] for i in sa_idx]
    sa_is_ww = [name.startswith(WW_PREFIXES) for name in sa_names]
    ww_sa_idx  = [j for j, ww in enumerate(sa_is_ww) if ww]

    sa_ww_mat = sa_mat[:, ww_sa_idx]

    sa_alt_count   = np.nansum(sa_mat == 1, axis=1).astype(int)
    sa_ww_has_alt  = np.any(sa_ww_mat == 1, axis=1)
    non_sa_has_alt = np.any(non_sa_mat == 1, axis=1)
    global_af      = np.nanmean(bin_mat, axis=1)

    alt_is_clean = np.array([
        all(b.upper() in valid_bases for b in alt.split(','))
        for alt in alts_col
    ])

    # observed in > 1 wastewater SA sample AND not globally common
    keep       = alt_is_clean & (sa_alt_count > 1)
    sa_private = keep & ~non_sa_has_alt

    print(f"SA-present SNPs (global AF ≤ {GLOBAL_AF_MAX}): {keep.sum()}")
    print(f"SA-private SNPs: {sa_private.sum()}")

    # Write SA-private TSV
    private_rows = []
    for i in np.where(sa_private)[0]:
        n_alt    = int(np.nansum(sa_mat[i] == 1))
        n_called = int(np.sum(~np.isnan(sa_mat[i])))
        private_rows.append({
            'POS': positions[i], 'REF': refs[i], 'ALT': alts_col[i],
            'n_SA_alt': n_alt, 'n_SA_called': n_called,
            'SA_AF': round(n_alt / n_called, 3) if n_called else np.nan,
        })

    # Build heatmap matrix: SA sequences only
    keep_idx     = np.where(sa_private)[0]
    hm           = bin_mat[np.ix_(keep_idx, sa_idx)].T   # (n_sa, n_sites)
    site_pos     = positions[keep_idx]
    site_private = sa_private[keep_idx]
    site_labels  = [f"{r}{p}{a}" for r, p, a in
                    zip(refs[keep_idx], site_pos, alts_col[keep_idx])]

    n_samp, n_sites = hm.shape
    print(f"Heatmap NaN (missing): {np.sum(np.isnan(hm))} / {hm.size} cells")

    names  = sa_names
    dates  = [sample_cols[i].split('|')[2] for i in sa_idx]
    is_ww  = sa_is_ww

    # Sort rows by sampling date, using mid-year for partial dates
    order  = sorted(range(n_samp), key=lambda i: parse_date(dates[i]))
    hm     = hm[order]
    is_ww  = [is_ww[i] for i in order]

    # Square cells: fixed cell size drives both axes dimensions
    cell   = 0.25         # inches per cell
    strip  = 0.2          # width of each annotation strip (inches)
    gene_h = 0.25         # height of gene annotation strip (inches)
    margin_l, margin_r = 0.3, 2.5
    margin_t, margin_b = 0.8 + gene_h + 0.05, 2.5   # extra room for gene strip above heatmap
    hm_w  = n_sites * cell
    hm_h  = n_samp  * cell
    fig_w = hm_w + 2 * strip + margin_l + margin_r
    fig_h = hm_h + margin_t + margin_b

    sorted_dates = [dates[i] for i in order]
    sorted_names = [names[i] for i in order]
    date_vals  = np.array([parse_date(d).toordinal() for d in sorted_dates], dtype=float)
    date_vmin  = pd.Timestamp('2024-01-01').toordinal()
    date_vmax  = pd.Timestamp('2025-08-31').toordinal()
    date_norm  = (date_vals - date_vmin) / (date_vmax - date_vmin)

    # Layout: margins expressed as figure fractions so axes are exactly cell-sized
    l = margin_l / fig_w
    r = 1 - margin_r / fig_w
    b = margin_b / fig_h
    t = 1 - margin_t / fig_h

    fig = plt.figure(figsize=(fig_w, fig_h))
    gs  = GridSpec(1, 3, width_ratios=[hm_w, strip, strip], wspace=0.02,
                   left=l, right=r, top=t, bottom=b)
    ax_hm   = fig.add_subplot(gs[0])
    ax_type = fig.add_subplot(gs[1])
    ax_date = fig.add_subplot(gs[2])

    # ── Gene region annotation strip (above heatmap) ─────────────────────
    gene_idx_arr = np.array([_get_gene_idx(p) for p in site_pos], dtype=float)
    gene_cmap    = matplotlib.colors.ListedColormap(_GENE_COLORS + [_INTERGENIC_CLR])
    hm_pos  = ax_hm.get_position()
    ax_gene = fig.add_axes([hm_pos.x0, hm_pos.y1 + 0.005,
                             hm_pos.width, gene_h / fig_h])
    ax_gene.pcolormesh(np.arange(n_sites + 1) - 0.5, np.array([-0.5, 0.5]),
                       gene_idx_arr.reshape(1, -1),
                       cmap=gene_cmap, vmin=0, vmax=len(GENE_REGIONS))
    ax_gene.set_xlim(-0.5, n_sites - 0.5)
    ax_gene.set_ylim(-0.5, 0.5)
    ax_gene.set_xticks([])
    ax_gene.set_yticks([])

    # Label each contiguous gene block; draw white dividers at boundaries
    prev_g, seg_start = gene_idx_arr[0], 0
    for j in range(1, n_sites):
        if gene_idx_arr[j] != prev_g:
            mid  = (seg_start + j - 1) / 2
            name = GENE_REGIONS[int(prev_g)][0] if int(prev_g) < len(GENE_REGIONS) else ''
            ax_gene.text(mid, 0, name, ha='center', va='center',
                         fontsize=7, fontweight='bold', color='white', clip_on=True)
            ax_gene.axvline(j - 0.5, color='white', linewidth=1)
            seg_start, prev_g = j, gene_idx_arr[j]
    if n_sites > 0:
        mid  = (seg_start + n_sites - 1) / 2
        name = GENE_REGIONS[int(prev_g)][0] if int(prev_g) < len(GENE_REGIONS) else ''
        ax_gene.text(mid, 0, name, ha='center', va='center',
                     fontsize=7, fontweight='bold', color='white', clip_on=True)

    # make heatmap
    # color coding 0=ref, 1=alt, 2=missing 
    cmap_hm = matplotlib.colors.ListedColormap(['#f0f0f0', '#c0392b', '#b0b0b0'])
    norm_hm = matplotlib.colors.BoundaryNorm([-0.5, 0.5, 1.5, 2.5], cmap_hm.N)
    hm_plot = np.where(np.isnan(hm), 2, hm)
    ax_hm.pcolormesh(np.arange(n_sites + 1) - 0.5, np.arange(n_samp + 1) - 0.5,
                     hm_plot, cmap=cmap_hm, norm=norm_hm)
    ax_hm.invert_yaxis()
    ax_hm.set_xticks(range(n_sites))
    ax_hm.set_xticklabels(site_labels, rotation=90, fontsize=10)
    ax_hm.set_xlabel('Mutation', fontsize=14, labelpad=4)
    ax_hm.set_yticklabels([])
    ax_hm.set_xticks(np.arange(-0.5, n_sites, 1), minor=True)
    # ax_hm.set_yticks(np.arange(-0.5, n_samp, 1), minor=False)
    ax_hm.set_yticks(np.arange(-0.5, n_samp, 1), minor=True)
    ax_hm.grid(which='minor', color='white', linewidth=0.5)
    ax_hm.tick_params(which='minor', length=0)

    def _strip_axes(ax):
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_yticks(np.arange(-0.5, n_samp, 1), minor=True)
        ax.grid(which='minor', color='white', linewidth=0.5, axis='y')
        ax.tick_params(which='minor', length=0)

    strip_X = np.array([-0.5, 0.5])
    strip_Y = np.arange(n_samp + 1) - 0.5

    # ── Sample type strip (clinical=red, wastewater=blue) ───
    cmap_type = matplotlib.colors.ListedColormap(['red', 'royalblue'])
    ax_type.pcolormesh(strip_X, strip_Y,
                       np.array(is_ww, dtype=float).reshape(-1, 1),
                       cmap=cmap_type, vmin=0, vmax=1)
    ax_type.invert_yaxis()
    _strip_axes(ax_type)

 # ── Collection date strip ─────────────────────────────────────────────
    cmap_date = plt.get_cmap('Oranges')
    ax_date.pcolormesh(strip_X, strip_Y, date_norm.reshape(-1, 1),
                       cmap=cmap_date, vmin=0, vmax=1)
    ax_date.invert_yaxis()
    _strip_axes(ax_date)

    # Draw X over cells where the month is unknown (date string contains 'XX')
    sorted_raw_dates = [dates[i] for i in order]
    for row_i, raw_date in enumerate(sorted_raw_dates):
        parts = raw_date.split('-')
        if len(parts) >= 2 and parts[1].upper() == 'XX':
            ax_date.plot([-0.4, 0.4], [row_i - 0.4, row_i + 0.4],
                         color='grey', linewidth=0.8, solid_capstyle='round')
            ax_date.plot([-0.4, 0.4], [row_i + 0.4, row_i - 0.4],
                         color='grey', linewidth=0.8, solid_capstyle='round')

    # ── Date colorbar (placed in right margin, aligned to date strip) ─────
    cbar_x = r + 0.01
    cbar_w = max(1.0 - r - 0.02, 0.06)
    cbar_ax = fig.add_axes([cbar_x-0.15, 0.97, cbar_w, 0.012])
    sm = matplotlib.cm.ScalarMappable(
        cmap=cmap_date,
        norm=matplotlib.colors.Normalize(vmin=date_vmin, vmax=date_vmax))
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cbar_ax, orientation='horizontal')
    tick_vals = np.linspace(date_vmin, date_vmax, 4)
    cbar.set_ticks(tick_vals)
    cbar.set_ticklabels([pd.Timestamp.fromordinal(int(v)).strftime('%Y-%m') for v in tick_vals],
                        fontsize=7)
    cbar.outline.set_visible(False)
    cbar.ax.set_title('Collection Date', fontsize=10, pad=4)

    # ── Legend (below plot) ───────────────────────────────────────────────
    fig.legend(handles=[
        mpatches.Patch(color='#c0392b', label='Alt allele'),
        mpatches.Patch(color='#f0f0f0', label='Ref allele', ec='#aaaaaa'),
        mpatches.Patch(color='#b0b0b0', label='Missing'),
    ], loc='lower left', bbox_to_anchor=(l, 0.99), ncol=3, fontsize=9, frameon=False)
    fig.legend(handles=[
        mpatches.Patch(color='red', label='Clinical'),
        mpatches.Patch(color='royalblue', label='Wastewater'),
    ], loc='lower right', bbox_to_anchor=(r, 0.99), ncol=2, fontsize=9, frameon=False)

    ax_gene.set_title(f'Genotype {gt_label}', fontsize=15, pad=15)

    out_pdf = f"../figures/{gt_label}_sa_snp_heatmap.pdf"
    fig.savefig(out_pdf, bbox_inches='tight')
    plt.close()
