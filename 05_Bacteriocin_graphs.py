CELL1
from google.colab import files
uploaded = files.upload()

CELL2


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
# 0. 
# ----------------------------------------------------------------------
EXCEL_PATH = "Bacteriocins.xlsx"   # <-- change path if needed
OUTDIR = "."                       # <-- where PNGs are saved

# Fixed category order + colors used across every figure (matches originals)
CATEGORY_ORDER = ["Fermented food", "Food", "Dairy", "Gut",
                   "Plant", "Fruit", "Other", "Unknown"]
CATEGORY_COLORS = {
    "Fermented food": "#D9541E",
    "Food":           "#B8860B",
    "Dairy":          "#F2A93B",
    "Gut":            "#0E6E5C",
    "Plant":          "#6FA83A",
    "Fruit":          "#2E86C1",
    "Other":          "#7F7F7F",
    "Unknown":        "#BDBDBD",
}

# Raw "Isolation Source" strings (Sheet3) -> the 8 manuscript categories.
# Feaces / Human Source / Environment / Animal / Insect / Probiotic / Fish /
# Vegetable are all pooled into "Other".
RAW_TO_CATEGORY = {
    "Fermented food": "Fermented food",
    "Food":           "Food",
    "Dairy":          "Dairy",
    "Gut":            "Gut",
    "Plant":          "Plant",
    "Fruit":          "Fruit",
    "Unknown":        "Unknown",
    "Feaces":         "Other",
    "Human Source":   "Other",
    "Environment":    "Other",
    "Animal":         "Other",
    "Insect":         "Other",
    "Probiotic":      "Other",
    "Fish":           "Other",
    "Vegetable":      "Other",
}

# Short display names + fixed column order for the class heatmap/stacked bar
CLASS_SHORT_NAMES = {
    "Class IIb (Two-peptide, synergistic)":       "IIb Two-peptide",
    "Class IIa (Pediocin-like, antilisterial)":   "IIa Pediocin-like",
    "Class IId (Non-Pediocin-like/Other Linear)": "IId Linear",
    "Class IIc (Circular)":                       "IIc Circular",
    "Class III  (Bacteriolysins)":                "III Bacteriolysin",
}
CLASS_ORDER = ["IIb Two-peptide", "IId Linear", "IIa Pediocin-like",
               "IIc Circular", "III Bacteriolysin"]

TOP_N_BACTERIOCINS = 10

# ----------------------------------------------------------------------
# 1. LOAD + MERGE DATA
# ----------------------------------------------------------------------
def load_data(excel_path=EXCEL_PATH):
    bgc = pd.read_excel(excel_path, sheet_name="Sheet1")
    bgc.columns = [c.strip() for c in bgc.columns]
    bgc["Bacteriocin_Class"] = bgc["Bacteriocin_Class"].str.strip()
    bgc["Class_Short"] = bgc["Bacteriocin_Class"].map(CLASS_SHORT_NAMES)

    meta = pd.read_excel(excel_path, sheet_name="Sheet3")
    meta.columns = [c.strip() for c in meta.columns]
    meta["Isolation Source"] = meta["Isolation Source"].str.strip()
    meta["Category"] = meta["Isolation Source"].map(RAW_TO_CATEGORY).fillna("Unknown")

    # NCBI accessions carry a version suffix in Sheet1 (GCF_xxx.1) but not
    # in the Sheet3 metadata export (GCF_xxx) -> strip before joining.
    bgc["GenomeBase"] = bgc["Genome"].str.split(".").str[0]
    meta["GenomeBase"] = meta["Genome"].str.split(".").str[0]

    genome_cat = meta[["GenomeBase", "Category"]].drop_duplicates("GenomeBase")
    bgc = bgc.merge(genome_cat, on="GenomeBase", how="left")
    bgc["Category"] = bgc["Category"].fillna("Unknown")

    # per-genome table: n_bgc = number of bacteriocin/BGC hits in that genome
    per_genome = (bgc.groupby(["GenomeBase", "Category"])
                     .size().rename("n_bgc").reset_index())

    return bgc, per_genome


# ----------------------------------------------------------------------
# FIGURE 1 — Mean BGC diversity by isolation source (bar + SD error bars)
# ----------------------------------------------------------------------
def fig1_mean_bgc(per_genome, outpath=f"{OUTDIR}/Fig1_BGC_mean_by_source.png"):
    stats = (per_genome.groupby("Category")["n_bgc"]
             .agg(["mean", "std", "count"])
             .reindex(CATEGORY_ORDER))
    overall_mean = per_genome["n_bgc"].mean()

    fig, ax = plt.subplots(figsize=(11, 7))
    colors = [CATEGORY_COLORS[c] for c in CATEGORY_ORDER]
    bars = ax.bar(stats.index, stats["mean"], yerr=stats["std"],
                   color=colors, capsize=5, error_kw={"linewidth": 1.3})

    for bar, n in zip(bars, stats["count"]):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() * 0.5,
                f"n={n}", ha="center", va="center",
                color="white", fontweight="bold", fontsize=12)

    ax.axhline(overall_mean, color="black", linestyle="--", linewidth=1.3)
    ax.text(len(CATEGORY_ORDER) - 0.4, overall_mean + 0.15,
            f"overall mean ({overall_mean:.2f})", ha="right", fontsize=11)

    ax.set_ylabel("Mean BGC types per genome")
    ax.set_title("Mean BGC diversity by isolation source", fontsize=16)
    ax.set_ylim(0, max(stats["mean"] + stats["std"]) * 1.25)
    plt.xticks(rotation=25, ha="right")
    ax.yaxis.grid(True, alpha=0.4)
    ax.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(outpath, dpi=200)
    plt.show()
    return stats


# ----------------------------------------------------------------------
# FIGURE 2 — Bacteriocin CLASS prevalence heatmap (% of genomes with >=1 hit)
# ----------------------------------------------------------------------
def fig2_class_heatmap(bgc, outpath=f"{OUTDIR}/Fig2_Class_heatmap_by_source.png"):
    genome_class = bgc[["GenomeBase", "Category", "Class_Short"]].drop_duplicates()
    n_genomes = bgc.groupby("Category")["GenomeBase"].nunique().reindex(CATEGORY_ORDER)

    prevalence = (genome_class.groupby(["Category", "Class_Short"])["GenomeBase"]
                  .nunique().unstack(fill_value=0)
                  .reindex(index=CATEGORY_ORDER, columns=CLASS_ORDER, fill_value=0))
    prevalence_pct = prevalence.div(n_genomes, axis=0) * 100

    fig, ax = plt.subplots(figsize=(11, 7))
    im = ax.imshow(prevalence_pct.values, cmap="Greens", vmin=0, vmax=100, aspect="auto")
    ax.set_xticks(range(len(CLASS_ORDER)))
    ax.set_xticklabels(CLASS_ORDER, rotation=25, ha="right")
    ax.set_yticks(range(len(CATEGORY_ORDER)))
    ax.set_yticklabels(CATEGORY_ORDER)

    for i in range(prevalence_pct.shape[0]):
        for j in range(prevalence_pct.shape[1]):
            val = prevalence_pct.values[i, j]
            color = "white" if val > 55 else "black"
            ax.text(j, i, f"{val:.0f}%", ha="center", va="center",
                     color=color, fontsize=12)

    cbar = fig.colorbar(im, ax=ax, shrink=0.85)
    cbar.set_label("% genomes")
    ax.set_title("Bacteriocin class prevalence by isolation source (%)", fontsize=16)
    fig.tight_layout()
    fig.savefig(outpath, dpi=200)
    plt.show()
    return prevalence_pct


# ----------------------------------------------------------------------
# FIGURE 3 — Top-10 individual bacteriocin prevalence heatmap
# ----------------------------------------------------------------------
def fig3_top10_heatmap(bgc, outpath=f"{OUTDIR}/Fig3_Bacteriocin_heatmap_by_source.png"):
    top10 = bgc["Bacteriocin"].value_counts().head(TOP_N_BACTERIOCINS).index.tolist()

    genome_bact = bgc[bgc["Bacteriocin"].isin(top10)][
        ["GenomeBase", "Category", "Bacteriocin"]].drop_duplicates()
    n_genomes = bgc.groupby("Category")["GenomeBase"].nunique().reindex(CATEGORY_ORDER)

    prevalence = (genome_bact.groupby(["Category", "Bacteriocin"])["GenomeBase"]
                  .nunique().unstack(fill_value=0)
                  .reindex(index=CATEGORY_ORDER, columns=top10, fill_value=0))
    prevalence_pct = prevalence.div(n_genomes, axis=0) * 100

    display_names = [b.replace("_", " ") for b in top10]

    fig, ax = plt.subplots(figsize=(14, 7))
    im = ax.imshow(prevalence_pct.values, cmap="Oranges", vmin=0, vmax=100, aspect="auto")
    ax.set_xticks(range(len(top10)))
    ax.set_xticklabels(display_names, rotation=30, ha="right", style="italic")
    ax.set_yticks(range(len(CATEGORY_ORDER)))
    ax.set_yticklabels(CATEGORY_ORDER)

    for i in range(prevalence_pct.shape[0]):
        for j in range(prevalence_pct.shape[1]):
            val = prevalence_pct.values[i, j]
            color = "white" if val > 60 else "black"
            ax.text(j, i, f"{val:.0f}%", ha="center", va="center",
                     color=color, fontsize=11)

    cbar = fig.colorbar(im, ax=ax, shrink=0.85)
    cbar.set_label("% genomes")
    ax.set_title("Top-10 bacteriocin prevalence by isolation source (%)", fontsize=16)
    fig.tight_layout()
    fig.savefig(outpath, dpi=200)
    plt.show()
    return prevalence_pct


# ----------------------------------------------------------------------
# FIGURE 4 — BGC richness distribution (violin + boxplot overlay)
# ----------------------------------------------------------------------
def fig4_richness_violin(per_genome, outpath=f"{OUTDIR}/Fig4_BGC_richness_violin_by_source.png"):
    overall_mean = per_genome["n_bgc"].mean()
    data = [per_genome.loc[per_genome["Category"] == cat, "n_bgc"].values
            for cat in CATEGORY_ORDER]
    counts = [len(d) for d in data]

    fig, ax = plt.subplots(figsize=(13, 7))
    positions = np.arange(1, len(CATEGORY_ORDER) + 1)

    vp = ax.violinplot(data, positions=positions, showextrema=False, widths=0.8)
    for body, cat in zip(vp["bodies"], CATEGORY_ORDER):
        body.set_facecolor(CATEGORY_COLORS[cat])
        body.set_alpha(0.55)
        body.set_edgecolor("black")
        body.set_linewidth(0.8)

    bp = ax.boxplot(data, positions=positions, widths=0.25, patch_artist=True,
                     medianprops={"color": "white", "linewidth": 2},
                     boxprops={"facecolor": "none", "edgecolor": "black"},
                     whiskerprops={"color": "black"}, capprops={"color": "black"},
                     flierprops={"markerfacecolor": "gray", "markeredgecolor": "gray",
                                 "markersize": 5})
    for patch, cat in zip(bp["boxes"], CATEGORY_ORDER):
        patch.set_facecolor(CATEGORY_COLORS[cat])

    ax.axhline(overall_mean, color="black", linestyle="--", linewidth=1.3)
    ax.text(len(CATEGORY_ORDER) + 0.55, overall_mean, f"mean {overall_mean:.2f}",
            va="center", fontsize=11)

    for pos, n in zip(positions, counts):
        ax.text(pos, -0.35, f"n={n}", ha="center", va="top", fontsize=10)

    ax.set_xticks(positions)
    ax.set_xticklabels(CATEGORY_ORDER, rotation=25, ha="right")
    ax.set_ylabel("BGC types per genome")
    ax.set_title("BGC richness distribution by isolation source", fontsize=16)
    ax.set_ylim(bottom=-0.8)
    ax.yaxis.grid(True, alpha=0.4)
    ax.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(outpath, dpi=200)
    plt.show()


# ----------------------------------------------------------------------
# FIGURE 5 — Bacteriocin class COMPOSITION (stacked %, of BGC hits not genomes)
# ----------------------------------------------------------------------
def fig5_class_stacked(bgc, outpath=f"{OUTDIR}/Fig5_Class_stacked_by_source.png"):
    comp = (bgc.groupby(["Category", "Class_Short"]).size()
            .unstack(fill_value=0)
            .reindex(index=CATEGORY_ORDER, columns=CLASS_ORDER, fill_value=0))
    comp_pct = comp.div(comp.sum(axis=1), axis=0) * 100

    fig, ax = plt.subplots(figsize=(13, 8))
    bottom = np.zeros(len(CATEGORY_ORDER))
    class_colors = {
        "IIb Two-peptide":   "#0E7A5F",
        "IId Linear":        "#B8860B",
        "IIa Pediocin-like": "#D9541E",
        "IIc Circular":      "#5B4B9E",
        "III Bacteriolysin": "#8B2E4E",
    }
    for cls in CLASS_ORDER:
        vals = comp_pct[cls].values
        bars = ax.bar(CATEGORY_ORDER, vals, bottom=bottom,
                       color=class_colors[cls], label=cls)
        for bar, v in zip(bars, vals):
            if v >= 5:
                ax.text(bar.get_x() + bar.get_width() / 2,
                        bar.get_y() + bar.get_height() / 2,
                        f"{v:.0f}%", ha="center", va="center",
                        color="white", fontweight="bold", fontsize=11)
        bottom += vals

    ax.set_ylabel("% of BGC hits")
    ax.set_title("Bacteriocin class composition by isolation source", fontsize=16)
    ax.set_ylim(0, 100)
    plt.xticks(rotation=25, ha="right")
    ax.legend(loc="upper right", bbox_to_anchor=(1.0, 1.0))
    ax.yaxis.grid(True, alpha=0.4)
    ax.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(outpath, dpi=200)
    plt.show()
    return comp_pct


# ----------------------------------------------------------------------
# RUN EVERYTHING
# ----------------------------------------------------------------------
if __name__ == "__main__":
    bgc, per_genome = load_data(EXCEL_PATH)

    print("Genomes by category:")
    print(per_genome["Category"].value_counts().reindex(CATEGORY_ORDER))

    stats1 = fig1_mean_bgc(per_genome)
    heat2 = fig2_class_heatmap(bgc)
    heat3 = fig3_top10_heatmap(bgc)
    fig4_richness_violin(per_genome)
    stack5 = fig5_class_stacked(bgc)

    print("\nAll 5 figures saved to:", OUTDIR)