CELL1

import subprocess, sys
subprocess.check_call([sys.executable, "-m", "pip", "install", "scipy", "matplotlib", "pandas", "numpy", "seaborn", "-q"])

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker
from matplotlib.gridspec import GridSpec
from scipy.optimize import curve_fit
import warnings, random
warnings.filterwarnings('ignore')


CELL2

matplotlib.rcParams.update({
    'font.family':       'sans-serif',
    'font.sans-serif':   ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size':         11,
    'axes.titlesize':    13,
    'axes.labelsize':    12,
    'axes.linewidth':    1.2,
    'axes.spines.top':   False,
    'axes.spines.right': False,
    'xtick.direction':   'out',
    'ytick.direction':   'out',
    'xtick.major.width': 1.2,
    'ytick.major.width': 1.2,
    'legend.frameon':    True,
    'legend.framealpha': 0.9,
    'legend.edgecolor':  '#cccccc',
    'legend.fontsize':   10,
    'figure.dpi':        150,
    'savefig.dpi':       300,
    'savefig.bbox':      'tight',
    'pdf.fonttype':      42,
    'svg.fonttype':      'none',
})

COLORS = {
    'pan':       '#2166ac',
    'core':      '#4dac26',
    'accessory': '#d01c8b',
    'unique':    '#f4a742',
    'new_genes': '#8856a7',
    'fit':       '#d7191c',
    'soft_core': '#80cdc1',
    'shell':     '#a6611a',
}

# ── 3. Upload CSV ────────────────────────────────────────────
from google.colab import files
print("Upload your Panaroo gene_presence_absence.csv:")
uploaded = files.upload()
csv_path = list(uploaded.keys())[0]
print(f"Loaded: {csv_path}")

# ── 4. Parse Panaroo CSV ─────────────────────────────────────
df = pd.read_csv(csv_path, low_memory=False)

meta_cols = [
    'Gene', 'Non-unique Gene name', 'Annotation',
    'No. isolates', 'No. sequences', 'Avg sequences per isolate',
    'Genome fragment', 'Order within fragment', 'Accessory Fragment',
    'Accessory Order with Fragment', 'QC',
    'Min group size nuc', 'Max group size nuc', 'Avg group size nuc'
]
genome_cols = [c for c in df.columns if c not in meta_cols]

# Binary presence/absence matrix (genes x genomes)
pa = df[genome_cols].notna().astype(int)
pa.index = df['Gene'].values

N_GENOMES = len(genome_cols)
N_GENES   = len(pa)
print(f"Genomes: {N_GENOMES}  |  Gene families: {N_GENES}")

# ── 5. USER SETTINGS — edit these ───────────────────────────
N_PERMUTATIONS  = 100     # resampling iterations
CORE_THRESHOLD  = 0.99    # fraction of genomes for core (0.99 = 99%)
SOFTCORE_THRESH = 0.95
SHELL_THRESH    = 0.15
SPECIES_NAME    = ""      # e.g. "Escherichia coli"
OUTPUT_PREFIX   = "pangenome"
FIGURE_FORMAT   = "svg"   # "svg", "png", or "pdf"
RANDOM_SEED     = 42

random.seed(RANDOM_SEED)
np.random.seed(RANDOM_SEED)

# ── 6. Resampling ────────────────────────────────────────────
def accumulate(pa_df, genomes):
    order = list(genomes)
    random.shuffle(order)
    pan_acc, core_acc, acc_acc, uniq_acc, new_acc = [], [], [], [], []
    seen = set()
    for i, g in enumerate(order):
        current = order[:i+1]
        n = i + 1
        sub = pa_df[current]
        freq = sub.sum(axis=1)
        present = set(sub.index[freq > 0])
        pan_acc.append(len(present))
        core_acc.append(int((freq == n).sum()))
        acc_acc.append(int(((freq > 1) & (freq < n)).sum()))
        uniq_acc.append(int((freq == 1).sum()))
        new_acc.append(len(present - seen))
        seen = present
    return {
        'pan':       np.array(pan_acc),
        'core':      np.array(core_acc),
        'accessory': np.array(acc_acc),
        'unique':    np.array(uniq_acc),
        'new_genes': np.array(new_acc),
    }

print(f"Running {N_PERMUTATIONS} permutations...")
results = []
for i in range(N_PERMUTATIONS):
    results.append(accumulate(pa, genome_cols))
    if (i+1) % 25 == 0:
        print(f"  {i+1}/{N_PERMUTATIONS}")

def agg(key):
    mat = np.vstack([r[key] for r in results])
    return dict(mean=mat.mean(0), sd=mat.std(0),
                q5=np.percentile(mat,5,0), q95=np.percentile(mat,95,0))

stats = {k: agg(k) for k in ['pan','core','accessory','unique','new_genes']}
x = np.arange(1, N_GENOMES + 1)

print("\nPangenome summary (all genomes):")
for k, label in [('pan','Pan'), ('core','Core'), ('accessory','Accessory'), ('unique','Unique')]:
    m, s = stats[k]['mean'][-1], stats[k]['sd'][-1]
    print(f"  {label:12s}: {m:.0f} ± {s:.0f} genes")

# ── 7. Heap's Law fitting ────────────────────────────────────
def heaps(n, k, gamma):    return k * np.power(n, gamma)
def power_decay(n, k, d):  return k * np.power(n, -d)

def fit(func, xd, yd, p0, bounds):
    try:
        popt, pcov = curve_fit(func, xd, yd, p0=p0, bounds=bounds, maxfev=20000)
        perr = np.sqrt(np.diag(pcov))
        yp   = func(xd, *popt)
        r2   = 1 - np.sum((yd-yp)**2) / np.sum((yd-yd.mean())**2)
        return popt, perr, r2
    except Exception as e:
        print(f"  Fit failed: {e}"); return None, None, None

pan_p,  pan_e,  pan_r2  = fit(heaps,       x,    stats['pan']['mean'],         [stats['pan']['mean'][0], 0.5],  ([1,0.01],[1e7,1.0]))
core_p, core_e, core_r2 = fit(power_decay, x,    stats['core']['mean'],        [stats['core']['mean'][0], 0.1], ([1,0.0], [1e7,5.0]))
new_p,  new_e,  new_r2  = fit(power_decay, x[1:],stats['new_genes']['mean'][1:],[stats['new_genes']['mean'][1], 0.5], ([0,0.0],[1e6,2.0]))

print("\n═══ Heap's Law Results ═══")
if pan_p is not None:
    k, g = pan_p
    status = "OPEN" if g >= 0.5 else "CLOSED"
    print(f"  Pan:  P(n) = {k:.2f} × n^{g:.4f}  |  γ={g:.4f}±{pan_e[1]:.4f}  R²={pan_r2:.4f}  → {status}")
if core_p is not None:
    print(f"  Core: C(n) = {core_p[0]:.2f} × n^-{core_p[1]:.4f}  |  R²={core_r2:.4f}")
if new_p is not None:
    print(f"  New:  G(n) = {new_p[0]:.2f} × n^-{new_p[1]:.4f}  |  R²={new_r2:.4f}")

# ── 8. Gene classification ───────────────────────────────────
gene_freq = pa.sum(axis=1) / N_GENOMES
core_genes     = gene_freq[gene_freq >= CORE_THRESHOLD].index
softcore_genes = gene_freq[(gene_freq >= SOFTCORE_THRESH) & (gene_freq < CORE_THRESHOLD)].index
shell_genes    = gene_freq[(gene_freq >= SHELL_THRESH) & (gene_freq < SOFTCORE_THRESH)].index
cloud_genes    = gene_freq[(gene_freq > 0) & (gene_freq < SHELL_THRESH)].index

print(f"\nGene classification:")
for name, genes in [('Core',core_genes),('Soft-core',softcore_genes),('Shell',shell_genes),('Cloud',cloud_genes)]:
    print(f"  {name:10s}: {len(genes):6d}  ({len(genes)/N_GENES*100:.1f}%)")

# ── 9. Figure 1: Pangenome accumulation curves ───────────────
title_str = f" — {SPECIES_NAME}" if SPECIES_NAME else ""
x_fit = np.linspace(1, N_GENOMES, 400)

fig = plt.figure(figsize=(14, 10))
gs  = GridSpec(2, 2, figure=fig, hspace=0.38, wspace=0.32)
fig.suptitle(f"Pangenome Analysis{title_str}", fontsize=15, fontweight='bold', y=1.01)

# Panel A: Pan + Core
ax1 = fig.add_subplot(gs[0, 0])
for key, label, col in [('pan','Pan-genome',COLORS['pan']),('core','Core genome',COLORS['core'])]:
    m = stats[key]['mean']
    ax1.fill_between(x, stats[key]['q5'], stats[key]['q95'], color=col, alpha=0.15)
    ax1.plot(x, m, color=col, lw=2.2, label=f"{label} (mean±95% CI)")
if pan_p is not None:
    ax1.plot(x_fit, heaps(x_fit, *pan_p), '--', color=COLORS['fit'], lw=1.6,
             label=f"Heap fit: P(n)={pan_p[0]:.1f}n$^{{{pan_p[1]:.3f}}}$, R²={pan_r2:.4f}")
if core_p is not None:
    ax1.plot(x_fit, power_decay(x_fit, *core_p), ':', color='#555', lw=1.6,
             label=f"Core fit: C(n)={core_p[0]:.1f}n$^{{-{core_p[1]:.3f}}}$, R²={core_r2:.4f}")
ax1.set_xlabel('Number of genomes'); ax1.set_ylabel('Gene families')
ax1.set_title('A   Pangenome accumulation curves')
ax1.legend(fontsize=8); ax1.set_xlim(1, N_GENOMES)
ax1.yaxis.set_minor_locator(ticker.AutoMinorLocator())
ax1.xaxis.set_minor_locator(ticker.AutoMinorLocator())

# Panel B: New genes per genome
ax2 = fig.add_subplot(gs[0, 1])
ax2.fill_between(x, stats['new_genes']['q5'], stats['new_genes']['q95'], color=COLORS['new_genes'], alpha=0.15)
ax2.plot(x, stats['new_genes']['mean'], color=COLORS['new_genes'], lw=2.2, label='New genes per genome')
if new_p is not None:
    ax2.plot(x_fit, power_decay(x_fit, *new_p), '--', color=COLORS['fit'], lw=1.6,
             label=f"Fit: G(n)={new_p[0]:.1f}n$^{{-{new_p[1]:.3f}}}$, R²={new_r2:.4f}")
ax2.set_xlabel('Number of genomes'); ax2.set_ylabel('New gene families added')
ax2.set_title('B   New genes per genome')
ax2.legend(fontsize=8); ax2.set_xlim(1, N_GENOMES)
ax2.yaxis.set_minor_locator(ticker.AutoMinorLocator())
ax2.xaxis.set_minor_locator(ticker.AutoMinorLocator())

# Panel C: Accessory + Unique
ax3 = fig.add_subplot(gs[1, 0])
for key, label, col in [('accessory','Accessory',COLORS['accessory']),('unique','Unique',COLORS['unique'])]:
    ax3.fill_between(x, stats[key]['q5'], stats[key]['q95'], color=col, alpha=0.15)
    ax3.plot(x, stats[key]['mean'], color=col, lw=2.2, label=label)
ax3.set_xlabel('Number of genomes'); ax3.set_ylabel('Gene families')
ax3.set_title('C   Accessory & unique gene content')
ax3.legend(fontsize=9); ax3.set_xlim(1, N_GENOMES)
ax3.yaxis.set_minor_locator(ticker.AutoMinorLocator())
ax3.xaxis.set_minor_locator(ticker.AutoMinorLocator())

# Panel D: Donut chart
ax4 = fig.add_subplot(gs[1, 1])
sizes  = [len(core_genes), len(softcore_genes), len(shell_genes), len(cloud_genes)]
labels = ['Core','Soft-core','Shell','Cloud']
cols   = [COLORS['core'], COLORS['soft_core'], COLORS['shell'], COLORS['accessory']]
wedges, texts, autotexts = ax4.pie(
    sizes, colors=cols,
    autopct=lambda p: f'{p:.1f}%' if p > 2 else '',
    startangle=90, pctdistance=0.75,
    wedgeprops=dict(width=0.55, edgecolor='white', linewidth=1.5)
)
for t in autotexts: t.set_fontsize(9); t.set_fontweight('bold')
legend_labels = [f"{l} ({s:,}, {s/N_GENES*100:.1f}%)" for l,s in zip(labels,sizes)]
ax4.legend(wedges, legend_labels, loc='center left', bbox_to_anchor=(0.85,0.5), fontsize=8)
ax4.set_title('D   Gene classification')

plt.savefig(f'{OUTPUT_PREFIX}_curves.{FIGURE_FORMAT}', bbox_inches='tight')
plt.savefig(f'{OUTPUT_PREFIX}_curves.png', dpi=300, bbox_inches='tight')
plt.show()
print(f"Saved: {OUTPUT_PREFIX}_curves.{FIGURE_FORMAT} + .png")


CELL3

# ── 11. Figure 3: Gene frequency histogram ───────────────────
fig, ax = plt.subplots(figsize=(9, 5))
freq_vals = pa.sum(axis=1).values
bins = np.arange(0.5, N_GENOMES + 1.5, 1)
bar_colors = []
for v in range(1, N_GENOMES + 1):
    f = v / N_GENOMES
    if   f >= CORE_THRESHOLD:  bar_colors.append(COLORS['core'])
    elif f >= SOFTCORE_THRESH: bar_colors.append(COLORS['soft_core'])
    elif f >= SHELL_THRESH:    bar_colors.append(COLORS['shell'])
    else:                      bar_colors.append(COLORS['accessory'])

counts, _, patches = ax.hist(freq_vals, bins=bins, edgecolor='white', linewidth=0.4)
for patch, col in zip(patches, bar_colors): patch.set_facecolor(col)
ax.legend(handles=[
    mpatches.Patch(color=COLORS['core'],      label=f'Core (≥{CORE_THRESHOLD*100:.0f}%)'),
    mpatches.Patch(color=COLORS['soft_core'], label=f'Soft-core ({SOFTCORE_THRESH*100:.0f}–{CORE_THRESHOLD*100:.0f}%)'),
    mpatches.Patch(color=COLORS['shell'],     label=f'Shell ({SHELL_THRESH*100:.0f}–{SOFTCORE_THRESH*100:.0f}%)'),
    mpatches.Patch(color=COLORS['accessory'], label=f'Cloud (<{SHELL_THRESH*100:.0f}%)'),
], fontsize=9)
ax.set_xlabel('Genomes containing gene family')
ax.set_ylabel('Number of gene families')
ax.set_title(f'Gene family frequency distribution{title_str}')
ax.set_xlim(0.5, N_GENOMES + 0.5)
ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
plt.tight_layout()
plt.savefig(f'{OUTPUT_PREFIX}_histogram.{FIGURE_FORMAT}', bbox_inches='tight')
plt.savefig(f'{OUTPUT_PREFIX}_histogram.png', dpi=300, bbox_inches='tight')
plt.show()
print(f"Saved: {OUTPUT_PREFIX}_histogram.{FIGURE_FORMAT} + .png")

# ── 12. Export CSV summary ───────────────────────────────────
results_df = pd.DataFrame({
    'n_genomes':      x,
    'pan_mean':       stats['pan']['mean'],   'pan_sd':  stats['pan']['sd'],
    'pan_q5':         stats['pan']['q5'],     'pan_q95': stats['pan']['q95'],
    'core_mean':      stats['core']['mean'],  'core_sd': stats['core']['sd'],
    'accessory_mean': stats['accessory']['mean'],
    'unique_mean':    stats['unique']['mean'],
    'new_genes_mean': stats['new_genes']['mean'],
})
results_df.to_csv(f'{OUTPUT_PREFIX}_stats.csv', index=False)
print(f"Saved: {OUTPUT_PREFIX}_stats.csv")

# ── 13. Download all files ───────────────────────────────────
from google.colab import files
import glob
for f in sorted(glob.glob(f'{OUTPUT_PREFIX}*')):
    print(f"Downloading {f}...")
    files.download(f)
print("\n✅ Done!")
