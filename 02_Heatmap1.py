CELL1
!pip install -q openpyxl

CELL2
from google.colab import files
uploaded = files.upload()
EXCEL_FILE = list(uploaded.keys())[0]
print(f"Uploaded: {EXCEL_FILE}")

CELL3
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list
from scipy.spatial.distance import pdist
import warnings
warnings.filterwarnings("ignore")

# ── Sheet names ───────────────────────────────────────────────────
GENES_SHEET   = "genes"
GENOMES_SHEET = "genomes"

# ── Column names in the genes sheet ──────────────────────────────
COL_GENE     = "Gene"
COL_FUNC_TAG = "Functional_tag"
COL_CATEGORY = "Category"

# ── Column names in the genomes sheet ────────────────────────────
COL_GENOME = "Genome"
COL_STRAIN = "Strain"
COL_SOURCE = "Source"

# ── Clustering ────────────────────────────────────────────────────
CLUSTER_METHOD = "ward"
CLUSTER_METRIC = "euclidean"   # euclidean works for continuous values

# ── Figure size ───────────────────────────────────────────────────
FIG_WIDTH  = 22
FIG_HEIGHT = 16

# ── Colourmap range (min/max of your data values) ────────────────
VMIN = 0
VMAX = 1    # change to 4 if your matrix has copy numbers, keep 1 for presence/absence

# ── Output filenames ──────────────────────────────────────────────
OUTPUT_PDF = "pangenome_heatmap.pdf"
OUTPUT_PNG = "pangenome_heatmap.png"

# ── Category colour palette ───────────────────────────────────────
CATEGORY_PALETTE = [
    "#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00",
    "#A65628","#F781BF","#999999","#66C2A5","#FC8D62",
    "#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494",
]

np.random.seed(42)
print("Settings ready.")


CELL4
xls        = pd.ExcelFile(EXCEL_FILE)
genes_df   = xls.parse(GENES_SHEET)
genomes_df = xls.parse(GENOMES_SHEET)

print("=== GENES sheet columns ===")
for i, c in enumerate(genes_df.columns.tolist()):
    print(f"  [{i}] '{c}'")
print("\n=== GENOMES sheet columns ===")
for i, c in enumerate(genomes_df.columns.tolist()):
    print(f"  [{i}] '{c}'")

def find_col(df, target):
    clean = {c.strip().lower(): c for c in df.columns}
    match = clean.get(target.strip().lower())
    if match is None:
        raise KeyError(
            f"\n  Column '{target}' NOT found.\n"
            f"  Available: {list(df.columns)}\n"
            f"  --> Update the COL_* variable in CELL 3."
        )
    return match

COL_GENE     = find_col(genes_df, COL_GENE)
COL_FUNC_TAG = find_col(genes_df, COL_FUNC_TAG)
COL_CATEGORY = find_col(genes_df, COL_CATEGORY)
COL_GENOME   = find_col(genomes_df, COL_GENOME)
COL_STRAIN   = find_col(genomes_df, COL_STRAIN)
COL_SOURCE   = find_col(genomes_df, COL_SOURCE)

print(f"\nMatched -> Gene:'{COL_GENE}'  FuncTag:'{COL_FUNC_TAG}'  Category:'{COL_CATEGORY}'")
print(f"           Genome:'{COL_GENOME}'  Strain:'{COL_STRAIN}'  Source:'{COL_SOURCE}'")

# Build matrix (rows=genes, cols=genomes)
meta_cols   = [COL_GENE, COL_FUNC_TAG, COL_CATEGORY]
genome_cols = [c for c in genes_df.columns if c not in meta_cols]

matrix_df = genes_df[genome_cols].copy()
matrix_df = matrix_df.apply(pd.to_numeric, errors="coerce").fillna(0)
matrix_df.index = genes_df[COL_GENE].values

gene_meta   = genes_df[[COL_GENE, COL_FUNC_TAG, COL_CATEGORY]].drop_duplicates(subset=COL_GENE).set_index(COL_GENE)
genome_meta = genomes_df.set_index(COL_GENOME)

print(f"\nGenes: {len(genes_df)}  |  Genomes: {len(genome_cols)}")
print(f"Categories: {genes_df[COL_CATEGORY].nunique()}  |  Func tags: {genes_df[COL_FUNC_TAG].nunique()}")
print(f"Data range: min={matrix_df.values.min():.2f}  max={matrix_df.values.max():.2f}")
print("Data loaded.")


CELL5

def safe_linkage(data, axis="row"):
    arr = data.values if axis == "row" else data.T.values
    try:
        dist = pdist(arr, metric=CLUSTER_METRIC)
    except Exception:
        dist = pdist(arr, metric="euclidean")
    dist = np.nan_to_num(dist, nan=0.0)
    return linkage(dist, method=CLUSTER_METHOD)

print("Clustering genes (rows) and genomes (columns)...")
row_link      = safe_linkage(matrix_df, axis="row")
col_link      = safe_linkage(matrix_df, axis="col")
row_order     = leaves_list(row_link)
col_order     = leaves_list(col_link)
sorted_matrix = matrix_df.iloc[row_order, col_order]
print("Clustering done.")

CELL6
categories_list = gene_meta[COL_CATEGORY].unique().tolist()
cat_colors = {c: CATEGORY_PALETTE[i % len(CATEGORY_PALETTE)] for i, c in enumerate(categories_list)}
print(f"Categories: {len(cat_colors)}")
for name, col in cat_colors.items():
    print(f"  {col}  {name}")
print("Colour maps ready.")


CELL7
print("Rendering figure...")
n_genes, n_genomes = sorted_matrix.shape

ordered_genes   = sorted_matrix.index.tolist()
ordered_genomes = [genome_cols[i] for i in col_order]

fig = plt.figure(figsize=(FIG_WIDTH, FIG_HEIGHT), dpi=150)
fig.patch.set_facecolor("white")

# Layout:
#   Row 0 (thin): top dendrogram
#   Row 1 (main): left_dendro | category_bar | heatmap | legend
gs = gridspec.GridSpec(
    2, 4, figure=fig,
    height_ratios=[0.15, 1],
    width_ratios=[0.10, 0.025, 1, 0.20],
    hspace=0.01, wspace=0.01,
)
ax_top_den  = fig.add_subplot(gs[0, 2])
ax_left_den = fig.add_subplot(gs[1, 0])
ax_cat_bar  = fig.add_subplot(gs[1, 1])
ax_heat     = fig.add_subplot(gs[1, 2])
ax_leg      = fig.add_subplot(gs[:, 3])
ax_leg.axis("off")

# ── Top dendrogram (over genome columns) ─────────────────────────────────────
dendrogram(col_link, ax=ax_top_den, orientation="top", color_threshold=0,
           above_threshold_color="#333333",
           link_color_func=lambda k: "#333333", no_labels=True)
ax_top_den.set_axis_off()

# ── Left dendrogram (over gene rows) ─────────────────────────────────────────
dendrogram(row_link, ax=ax_left_den, orientation="left", color_threshold=0,
           above_threshold_color="#333333",
           link_color_func=lambda k: "#333333", no_labels=True)
ax_left_den.set_axis_off()

# ── Main heatmap — continuous white→dark-blue scale ──────────────────────────
blue_cmap = LinearSegmentedColormap.from_list(
    "white_to_navy",
    ["#FFFFFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"],
    N=256
)
im = ax_heat.imshow(
    sorted_matrix.values,
    aspect="auto",
    cmap=blue_cmap,
    vmin=VMIN, vmax=VMAX,
    interpolation="nearest",
)

# Gene names on RIGHT y-axis (like reference image)
ax_heat.yaxis.set_label_position("right")
ax_heat.yaxis.tick_right()
ax_heat.set_yticks(range(n_genes))
gene_fontsize = max(2.5, min(7, 200 / n_genes))  # auto-scale font to number of genes
ax_heat.set_yticklabels(ordered_genes, fontsize=gene_fontsize, fontfamily="monospace")
ax_heat.tick_params(axis="y", which="both", length=0, pad=2)

# Genome/strain names on BOTTOM x-axis, rotated (like reference image)
ax_heat.set_xticks(range(n_genomes))
genome_labels = []
for g in ordered_genomes:
    label = str(genome_meta.loc[g, COL_STRAIN]) if g in genome_meta.index else g
    genome_labels.append(label)
genome_fontsize = max(2.5, min(7, 300 / n_genomes))
ax_heat.set_xticklabels(genome_labels, fontsize=genome_fontsize, rotation=90, ha="center", va="top")
ax_heat.tick_params(axis="x", which="both", length=0, pad=2)

for sp in ax_heat.spines.values(): sp.set_visible(False)

# ── Category bar (left of heatmap, one colour per gene row) ──────────────────
def safe_cat(gene):
    if gene not in gene_meta.index: return "Unknown"
    val = gene_meta.loc[gene, COL_CATEGORY]
    return str(val.iloc[0]) if hasattr(val, "iloc") else str(val)

cat_vals = [safe_cat(g) for g in ordered_genes]
cat_rgb  = np.array([[mcolors.to_rgb(cat_colors.get(c, "#CCCCCC"))] for c in cat_vals])
ax_cat_bar.imshow(cat_rgb, aspect="auto", interpolation="none")
ax_cat_bar.set_yticks([])
ax_cat_bar.set_xticks([0])
ax_cat_bar.set_xticklabels(["Category"], fontsize=6, rotation=90, va="top")
ax_cat_bar.tick_params(bottom=False)
for sp in ax_cat_bar.spines.values(): sp.set_visible(False)

# ── Colourbar ─────────────────────────────────────────────────────────────────
# Place it inside the legend axes
cbar_ax = fig.add_axes([0.815, 0.62, 0.018, 0.25])
cb = fig.colorbar(im, cax=cbar_ax)
cb.set_ticks([VMIN, (VMIN+VMAX)/2, VMAX])
cb.ax.tick_params(labelsize=6)
cb.outline.set_visible(False)
label_txt = "Presence/Absence (0/1)" if VMAX == 1 else "Gene copy number"
cb.set_label(label_txt, fontsize=6, labelpad=4)

# ── Category legend ───────────────────────────────────────────────────────────
leg_left   = 0.815   # figure x coordinate
leg_top    = 0.60    # figure y coordinate (starts below colourbar)
dy         = 0.042
fig.text(leg_left, leg_top, "Category", fontsize=7, fontweight="bold",
         transform=fig.transFigure, va="top")
ly = leg_top - 0.04
for cat, color in cat_colors.items():
    fig.patches.append(mpatches.Rectangle(
        (leg_left, ly - 0.014), 0.016, 0.022,
        transform=fig.transFigure, color=color, clip_on=False
    ))
    fig.text(leg_left + 0.022, ly - 0.003, cat[:34],
             transform=fig.transFigure, fontsize=5.5, va="center")
    ly -= dy
    if ly < 0.04: break

# ── Title ─────────────────────────────────────────────────────────────────────
fig.suptitle(
    "Pangenome Presence/Absence — Neuroactive Genes",
    fontsize=12, fontweight="bold", y=0.998, color="#111111"
)

plt.show()
print("Figure rendered successfully.")


CELL8
from google.colab import files

fig.savefig(OUTPUT_PDF, dpi=300, bbox_inches="tight", facecolor="white", format= "png")
print(f"Saved: {OUTPUT_PDF}")
print(f"Saved: {OUTPUT_PNG}")

files.download(OUTPUT_PDF)
files.download(OUTPUT_PNG)
print("Download started!")
