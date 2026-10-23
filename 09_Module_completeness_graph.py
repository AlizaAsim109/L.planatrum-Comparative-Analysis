# ============================================================
# 1. SETUP
# ============================================================
import sys, subprocess

def _pip_install(pkgs):
    subprocess.run([sys.executable, "-m", "pip", "install", "-q", *pkgs], check=True)

try:
    import openpyxl  # noqa
except ImportError:
    _pip_install(["openpyxl"])

import textwrap
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

IN_COLAB = "google.colab" in sys.modules
print("Running in Google Colab:", IN_COLAB)
print("pandas", pd.__version__, "| matplotlib", matplotlib.__version__)


CELL2
# ============================================================
# 2. UPLOAD YOUR WORKBOOK
# ============================================================
if IN_COLAB:
    from google.colab import files
    print("Select your GBM_completeness_from_eggNOG.xlsx file...")
    uploaded = files.upload()
    assert len(uploaded) >= 1, "No file uploaded."
    workbook_path = list(uploaded.keys())[0]
else:
    workbook_path = "GBM_completeness_from_eggNOG.xlsx"

print("Using workbook:", workbook_path)


CELL3
# ============================================================
# 3. LOAD Module_Completeness_Summary AND FILTER TO MODULES WITH DATA
# ============================================================
SHEET_NAME = "Module_Completeness_Summary"
HEADER_ROW_1INDEXED = 3  # row 1: free-text note; row 2: blank; row 3: real column headers

xl = pd.ExcelFile(workbook_path)
assert SHEET_NAME in xl.sheet_names, (
    f"Sheet '{SHEET_NAME}' not found. Sheets present: {xl.sheet_names}. "
    f"Update SHEET_NAME above to match your workbook."
)
raw = xl.parse(SHEET_NAME, header=HEADER_ROW_1INDEXED - 1)
raw = raw.dropna(how="all")

required_cols = {"Module_ID", "Module_name", "N_required_steps", "Mean_completeness_pct"}
missing_cols = required_cols - set(raw.columns)
assert not missing_cols, (
    f"Expected columns {sorted(missing_cols)} not found in '{SHEET_NAME}'. "
    f"Columns present: {list(raw.columns)}."
)

df = raw[raw["Mean_completeness_pct"] > 0].copy().reset_index(drop=True)
assert len(df) > 0, "No modules with completeness > 0% were found - nothing to plot."
print(f"[ok] {len(raw)} total modules in the sheet; {len(df)} have completeness > 0% and will be plotted.")
print(df[["Module_ID", "Module_name", "N_required_steps", "Mean_completeness_pct"]]
      .sort_values("Mean_completeness_pct", ascending=False).to_string(index=False))



CELL4
# ============================================================
# 4. ASSIGN FUNCTIONAL CATEGORIES (by Module_ID, not by guessing from the name)
# ============================================================
# Categories are this notebook's own functional grouping for figure readability -
# NOT an officially published GBM taxonomy. Update/extend this dict if your dataset
# surfaces different modules; the assertion below will tell you exactly which Module_ID
# is missing rather than silently mis-plotting it.
CATEGORY_MAP = {
    "MGB004": "Amino-acid-derived\nneuroactive precursors",
    "MGB005": "Amino-acid-derived\nneuroactive precursors",
    "MGB006": "Amino-acid-derived\nneuroactive precursors",
    "MGB022": "GABAergic metabolism",
    "MGB019": "GABAergic metabolism",
    "MGB032": "NAD / quinolinate\nmetabolism",
    "MGB033": "NAD / quinolinate\nmetabolism",
    "MGB031": "Steroid metabolism",
    "MGB036": "Vitamin / cofactor\nbiosynthesis",
    "MGB038": "Vitamin / cofactor\nbiosynthesis",
    "MGB037": "Vitamin / cofactor\nbiosynthesis",
    "MGB040": "Vitamin / cofactor\nbiosynthesis",
    "MGB041": "Vitamin / cofactor\nbiosynthesis",
    "MGB029": "Stress-response\nchaperone",
    "MGB043": "Short-chain fatty acid\nmetabolism",
    "MGB044": "Short-chain fatty acid\nmetabolism",
    "MGB045": "Short-chain fatty acid\nmetabolism",
    "MGB052": "Short-chain fatty acid\nmetabolism",
    "MGB053": "Short-chain fatty acid\nmetabolism",
    "MGB034": "Branched-chain fatty\nacid metabolism",
    "MGB035": "Branched-chain fatty\nacid metabolism",
}
CATEGORY_ORDER = [
    "Amino-acid-derived\nneuroactive precursors", "GABAergic metabolism",
    "NAD / quinolinate\nmetabolism", "Steroid metabolism",
    "Vitamin / cofactor\nbiosynthesis", "Stress-response\nchaperone",
    "Short-chain fatty acid\nmetabolism", "Branched-chain fatty\nacid metabolism",
]

unmapped = sorted(set(df["Module_ID"]) - set(CATEGORY_MAP.keys()))
assert not unmapped, (
    f"Module(s) {unmapped} appeared in your data but have no entry in CATEGORY_MAP. "
    f"Add them to CATEGORY_MAP in this cell (pick an existing category, or a new one - "
    f"remember to also add any new category name to CATEGORY_ORDER), then re-run from here."
)

df["Category"] = df["Module_ID"].map(CATEGORY_MAP)

def status(pct):
    if pct == 100:
        return "Complete"
    if pct >= 50:
        return "Near-complete"
    return "Incomplete"

df["Status"] = df["Mean_completeness_pct"].apply(status)
df["cat_rank"] = df["Category"].apply(lambda c: CATEGORY_ORDER.index(c))
df = df.sort_values(["cat_rank", "Mean_completeness_pct"], ascending=[True, True]).reset_index(drop=True)

# minor label-space abbreviation only - does not change any underlying value
df["Module_name"] = df["Module_name"].replace(
    "Menaquinone synthesis (vitamin K2) II (alternative pathway: futalosine pathway)",
    "Menaquinone synthesis (vitamin K2) II (futalosine pathway)",
)

print("[ok] Categories assigned. Plot order:")
print(df[["Module_ID", "Module_name", "Category", "Mean_completeness_pct"]].to_string(index=False))


CELL5
# ============================================================
# 5. BUILD THE FIGURE
# ============================================================
plt.rcParams.update({"font.family": "DejaVu Serif", "font.size": 10,
                      "pdf.fonttype": 42, "svg.fonttype": "none"})

STATUS_COLORS = {"Complete": "#1b6b45", "Near-complete": "#c98a1a", "Incomplete": "#a4302a"}

# Pleasant, distinguishable pastel band colors + matching darker text colors, auto-generated
# per category actually present (so this still works if CATEGORY_MAP above is edited/extended).
_BAND_PALETTE = ["#eef3fb", "#eaf6ee", "#fbeeee", "#f7f0e3", "#f1edf7",
                  "#eef7f6", "#fdf3e7", "#f3f3e9", "#eef0f7", "#f6eef2"]
_TEXT_PALETTE = ["#2a4d8f", "#1b6b45", "#9c3b36", "#8a6a1f", "#5b3f82",
                  "#1a6e68", "#9c6a1f", "#6d6b2f", "#3c4a8a", "#8a3c63"]

cats_present = []
for c in df["Category"]:
    if c not in cats_present:
        cats_present.append(c)
cat_band_colors = {c: _BAND_PALETTE[i % len(_BAND_PALETTE)] for i, c in enumerate(cats_present)}
cat_text_colors = {c: _TEXT_PALETTE[i % len(_TEXT_PALETTE)] for i, c in enumerate(cats_present)}

plan = []
for c in cats_present:
    plan.append({"type": "header", "category": c})
    sub = df[df["Category"] == c]
    for _, row in sub.iterrows():
        plan.append({"type": "module", "category": c, **row.to_dict()})

n_rows = len(plan)
y_pos = np.arange(n_rows)[::-1]

fig = plt.figure(figsize=(14.2, 0.40 * n_rows + 2.4))
ax = fig.add_axes([0.41, 0.14, 0.44, 0.74])

i = 0
while i < n_rows:
    if plan[i]["type"] == "header":
        c = plan[i]["category"]
        j = i
        while j + 1 < n_rows and plan[j + 1]["type"] == "module":
            j += 1
        y_top = y_pos[i] + 0.5
        y_bot = y_pos[j] - 0.5
        ax.axhspan(y_bot, y_top, color=cat_band_colors[c], zorder=0, clip_on=False)
        i = j + 1
    else:
        i += 1

bar_height = 0.64
for idx, item in enumerate(plan):
    yi = y_pos[idx]
    if item["type"] == "header":
        ax.text(-1.0, yi, item["category"].replace("\n", " "), va="center", ha="left",
                fontsize=9.8, fontweight="bold", color=cat_text_colors[item["category"]], clip_on=False)
    else:
        pct = item["Mean_completeness_pct"]
        steps = int(item["N_required_steps"])
        color = STATUS_COLORS[item["Status"]]
        ax.barh(yi, 100, height=bar_height, color="#e8e6e0", zorder=1, edgecolor="none")
        ax.barh(yi, pct, height=bar_height, color=color, zorder=2, edgecolor="white", linewidth=0.6)
        label = f"{pct:.0f}%" if pct == int(pct) else f"{pct:.1f}%"
        ax.text(pct + 2, yi, f"{label}  ({steps} step{'s' if steps > 1 else ''})",
                va="center", ha="left", fontsize=8.7, color="#333333", zorder=4, clip_on=False)
        name_label = f"{item['Module_ID']} \u2013 {item['Module_name']}"
        ax.text(-1.0, yi, name_label, va="center", ha="right", fontsize=8.8, color="#1a1a1a", clip_on=False)

ax.set_yticks([])
ax.set_xlim(0, 100)
ax.set_ylim(-0.7, n_rows - 0.3)
ax.set_xticks([0, 25, 50, 75, 100])
ax.set_xticklabels(["0%", "25%", "50%", "75%", "100%"], fontsize=9.3)
ax.set_xlabel("Functional-step completeness (% of official GBM steps reconstructed)", fontsize=10, labelpad=9)
for spine in ["top", "right", "left"]:
    ax.spines[spine].set_visible(False)
ax.spines["bottom"].set_color("#888888")
ax.tick_params(axis="y", length=0)
ax.tick_params(axis="x", length=3, color="#888888")

legend_elems = [mpatches.Patch(color=STATUS_COLORS["Complete"], label="Complete (100%)"),
                mpatches.Patch(color=STATUS_COLORS["Near-complete"], label="Near-complete (\u226550%)"),
                mpatches.Patch(color=STATUS_COLORS["Incomplete"], label="Incomplete (<50%)")]
leg = fig.legend(handles=legend_elems, title="Reconstruction status", loc="lower center",
                  frameon=True, fontsize=9.3, title_fontsize=9.8, ncol=3, bbox_to_anchor=(0.63, 0.012))
leg.get_frame().set_edgecolor("#cccccc")

n_present = len(df)
n_total = len(raw)
fig.suptitle("Pangenome-level reconstruction of Gut-Brain Module (GBM) functional potential",
             fontsize=14, fontweight="bold", x=0.63, y=0.975)
fig.text(0.63, 0.925,
         f"{n_present} of {n_total} official GBM modules with \u22651 required functional step supported by the "
         "combined, non-redundant\npangenome gene catalogue (eggNOG-mapper KO annotation; AND/OR step logic per "
         "Valles-Colomer et al., 2019)",
         ha="center", va="top", fontsize=9.3, color="#444444", style="italic")

print("[ok] Figure built.")


CELL6
# ============================================================
# 6. SAVE AND DOWNLOAD
# ============================================================
png_path = "Figure_GBM_pangenome_reconstruction.png"
pdf_path = "Figure_GBM_pangenome_reconstruction.pdf"
fig.savefig(png_path, dpi=350)
fig.savefig(pdf_path)
print(f"[ok] Saved {png_path} (350 dpi, for drafts/Word) and {pdf_path} (vector, for journal submission).")

if IN_COLAB:
    from google.colab import files
    files.download(png_path)
    files.download(pdf_path)