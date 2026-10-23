CELL 1

import sys, subprocess

def _pip_install(pkgs):
    subprocess.run([sys.executable, "-m", "pip", "install", "-q", *pkgs], check=True)

try:
    import openpyxl  # noqa
except ImportError:
    _pip_install(["openpyxl"])

try:
    import xlsxwriter  # noqa
except ImportError:
    _pip_install(["xlsxwriter"])

import io, re, json, base64, hashlib, textwrap, warnings
from collections import OrderedDict, defaultdict
import numpy as np
import pandas as pd
import requests

warnings.filterwarnings("ignore", category=UserWarning, module="openpyxl")

IN_COLAB = "google.colab" in sys.modules
print("Running in Google Colab:", IN_COLAB)
print("pandas", pd.__version__, "| numpy", np.__version__)


CELL 2
# ============================================================
# 2. CONFIG  (edit only if your file layout is non-standard)
# ============================================================
CONFIG = {
    # --- eggNOG-mapper input ---
    # Column in your eggNOG export that holds KEGG Orthology hits.
    # Standard eggNOG-mapper column name is "KEGG_ko" (values like "ko:K00001,ko:K00002" or "-").
    # If your column is named differently, add it to this list (checked in order, case-insensitive).
    "ko_column_candidates": ["KEGG_ko", "KEGG_KOs", "KO", "KOs", "ko"],

    # How to split multiple genomes out of the workbook:
    #   "auto"        -> if the workbook has >1 sheet, one sheet = one genome (sheet name = genome ID);
    #                     if it has 1 sheet, look for a genome/sample ID column (see genome_id_candidates);
    #                     if none found, the whole sheet is treated as ONE genome.
    #   "per_sheet"   -> force: one sheet = one genome, regardless of sheet count.
    #   "single_genome" -> force: treat the entire input as ONE genome (use if you have exactly one
    #                       eggNOG output file for one genome/MAG).
    "genome_split_mode": "auto",

    # Only used when genome_split_mode == "auto" and there is exactly one sheet:
    "genome_id_candidates": ["genome_id", "genome", "sample", "sample_id", "assembly", "strain", "genome_name"],

    # Rows in raw eggNOG-mapper .emapper.annotations output sometimes carry '##' comment/header lines
    # when exported directly (e.g. '## whatever' at top, '## N queries scanned' at bottom). These are
    # dropped automatically by column-name detection below; no action needed.

    # --- GBM database ---
    # Pinned, reproducible source (exact commit, will not silently change):
    "gbm_url": ("https://raw.githubusercontent.com/omixer/omixer-rpmR/"
                "184f1cc99aefed722e20eb00eda082348a064c4e/inst/extdata/GBMs.v1.0.txt"),
    # Expected SHA-256 of that exact file, used to verify the download is byte-identical to the
    # embedded offline fallback (see cell 3). If this ever fails, the notebook automatically falls
    # back to the embedded copy rather than silently using a possibly-different file.
    "gbm_sha256_expected": None,  # computed and printed after cell 3 runs; filled in automatically

    # --- output ---
    "output_xlsx": "GBM_completeness_from_eggNOG.xlsx",
    "completeness_analytical_threshold_pct": 80,  # NOT an official GBM cutoff; reported only as a label
}
print(json.dumps({k: v for k, v in CONFIG.items() if k != "gbm_sha256_expected"}, indent=2))


CELL 3

# ============================================================
# 3. GET THE OFFICIAL GBM DATABASE  (network fetch with verified offline fallback)
# ============================================================
# Embedded verbatim copy of GBMs.v1.0.txt (base64-encoded to avoid any whitespace/escaping issues),
# retrieved from https://github.com/omixer/omixer-rpmR, commit 184f1cc99aefed722e20eb00eda082348a064c4e.
# This guarantees the notebook can always reproduce the exact same analysis even with no internet
# access, and lets us verify that any freshly downloaded copy is byte-identical before trusting it.
_GBM_EMBEDDED_B64 = (
    "TUdCMDAxCVNlcm90b25pbiBzeW50aGVzaXMgSQpDT0cyNTA5CkswMTU5MwlDT0cwMDc2Ci8vLwpNR0IwMDIJU2Vyb3RvbmluIHN5bnRoZXNpcyBJSQpLMDE1OTMJQ09HMDA3NgpLMDA0OTIKLy8vCk1HQjAwMwlNZWxhdG9uaW4gc3ludGhlc2lzCmJhY3ROT0cyNjMxNglwcm9OT0czMDE5MQpLMDA1NDMJQ09HMDUwMAovLy8KTUdCMDA0CUt5bnVyZW5pbmUgc3ludGhlc2lzCkswMDQ1MwlLMDA0NjMJYmFjdE5PRzEyNDUzCUNPRzM0ODMKSzAxNDMyCUswNzEzMAlDT0cxODc4Ci8vLwpNR0IwMDUJVHJ5cHRvcGhhbiBzeW50aGVzaXMKSzAxNjU3LEswMTY1OAlLMTM0OTcJSzEzNTAxCUsxMzUwMwpLMDA3NjYJSzEzNDk3CkswMTgxNwlLMTM0OTgJSzEzNTAxCkswMTYwOQlLMTM0OTgJSzEzNTAxCkswMTY5NSxLMDE2OTYJSzAxNjk1LEswNjAwMQovLy8KTUdCMDA2CUdsdXRhbWF0ZSBzeW50aGVzaXMgSQpLMDAyNjEJSzAwMjYyCi8vLwpNR0IwMDcJR2x1dGFtYXRlIHN5bnRoZXNpcyBJSQpLMDAyNjUsSzAwMjY2CUswMDI4NAovLy8KTUdCMDA4CUFjZXR5bGNob2xpbmUgc3ludGhlc2lzCkswMDYyNAliYWN0Tk9HODI2MTcJYmFjdE5PRzAyMTU5Ci8vLwpNR0IwMDkJSGlzdGFtaW5lIHN5bnRoZXNpcwpLMDE1OTAJVElHUjAwNTQxCWJhY3ROT0cwNTgyOQliYWN0Tk9HMTY2MjcKLy8vCk1HQjAxMAlIaXN0YW1pbmUgZGVncmFkYXRpb24KSzAwMjc2CWJhY3ROT0cwODM0MAliYWN0Tk9HMDE2NjIKLy8vCk1HQjAxMQlDb3J0aXNvbCBkZWdyYWRhdGlvbgpDT0czOTU5CmJhY3ROT0cwMTg0NApDT0cxMDYzCmJhY3ROT0cwMzI2OQovLy8KTUdCMDEyCURvcGFtaW5lIHN5bnRoZXNpcwpLMDA1MDUJVElHUjAxMjY5CWJhY3ROT0cxNDA5OAliYWN0Tk9HMjAyNDkKSzAxNTkzCUswMTU5MglUSUdSMDM4MTEJVElHUjAzODEyCUNPRzAwNzYKLy8vCk1HQjAxMwlOb3JhZHJlbmFsaW5lIHN5bnRoZXNpcwpOT0cxMzI1NTMKLy8vCk1HQjAxNAlBZHJlbmFsaW5lIHN5bnRoZXNpcwpOT0cxMzM2NjMKLy8vCk1HQjAxNQlwLUNyZXNvbCBzeW50aGVzaXMKSzE4NDI3LEsxODQyOAlLMTU3NjAsSzE1NzYxLEsxNTc2MixLMTU3NjMsSzE1NzY0CUswMzE1MAliYWN0Tk9HMjg4NzYKLy8vCk1HQjAxNglwLUNyZXNvbCBkZWdyYWRhdGlvbgpLMDU3OTcJSzAzMzgwCUswMDQ0OCxLMDA0NDkJVElHUjAyNDIzLFRJR1IwMjQyMglLMDA0NDgsVElHUjAyNDIyCVRJR1IwMjQyMyxLMDA0NDkKLy8vCk1HQjAxNwlJUEEgc3ludGhlc2lzIChDbG9zdHJpZGl1bSBzcG9yb2dlbmVzIHBhdGh3YXkpCmZpcm1OT0cwNDI5MApmaXJtTk9HMDA2MjYsYmFjdE5PRzA3ODgxLGJhY3ROT0cxNTM0MQovLy8KTUdCMDE4CUt5bnVyZW5pbmUgZGVncmFkYXRpb24KSzAwNDg2CWJhY3ROT0cwMjU5MgovLy8KTUdCMDE5CUdBQkEgZGVncmFkYXRpb24KSzA3MjUwCUsxNDI2OAlLMDA4MjMJSzEzNTI0CUswMDYxMwlLMTY4NzEKLy8vCk1HQjAyMAlHQUJBIHN5bnRoZXNpcyBJCkswOTQ3MApLMDk0NzEKSzA5NDcyCkswOTQ3MwovLy8KTUdCMDIxCUdBQkEgc3ludGhlc2lzIElJCkswOTI1MQpLMDAxMzcKLy8vCk1HQjAyMglHQUJBIHN5bnRoZXNpcyBJSUkKSzAxNTgwCWJhY3ROT0cwMTc5MwovLy8KTUdCMDIzCURvcGFtaW5lIGRlZ3JhZGF0aW9uCkswMTAxNAlLMDA1NDUJSzAwMjc0CUswMDI3NglLMTMzNzEsSzEzMzcyCWJhY3ROT0cwODM0MAovLy8KTUdCMDI0CURPUEFDIHN5bnRoZXNpcwpLMDAyNzQJSzAwMjc2CUsxMzM3MSxLMTMzNzIJYmFjdE5PRzA4MzQwCkswMDEyOAlLMDAxMjkJQ09HMTAxMgovLy8KTUdCMDI1CU5pdHJpYyBveGlkZSBzeW50aGVzaXMgSSAoTk8gc3ludGhhc2UpCkswMDQ5MQliYWN0Tk9HMTI4MzYJYmFjdE5PRzUxNjc4Ci8vLwpNR0IwMjYJTml0cmljIG94aWRlIHN5bnRoZXNpcyBJSSAobml0cml0ZSByZWR1Y3Rhc2UpCkswMDM2OAlLMTU4NjQJYmFjdE5PRzEwMjQ3CUNPRzMzMDMKLy8vCk1HQjAyNwlOaXRyaWMgb3hpZGUgZGVncmFkYXRpb24gSSAoTk8gZGlveHlnZW5hc2UpCkswNTkxNgliYWN0Tk9HMDUwNTYJYmFjdE5PRzM2OTg0Ci8vLwpNR0IwMjgJTml0cmljIG94aWRlIGRlZ3JhZGF0aW9uIElJIChOTyByZWR1Y3Rhc2UpCkswNDU2MSxLMDIzMDUJYmFjdE5PRzAzNzY2LGJhY3ROT0cxODE5MglLMDQ1NjEsYmFjdE5PRzE4MTkyCWJhY3ROT0cwMzc2NixLMDIzMDUKLy8vCk1HQjAyOQlDbHBCIChBVFAtZGVwZW5kZW50IGNoYXBlcm9uZSBwcm90ZWluKQpLMDM2OTUJVElHUjAzMzQ2CUNPRzA1NDIKLy8vCk1HQjAzMAlQb2x5c2FjY2hhcmlkZSBBCmJhY3ROT0cwMDAxNApiYWN0Tk9HMTg1MTkKYmFjdE5PRzI0NzExCWJhY3ROT0cxODYzMAliYWN0Tk9HMzExMTYKLy8vCk1HQjAzMQkxNy1iZXRhLUVzdHJhZGlvbCBkZWdyYWRhdGlvbgpDT0cwMzAwCUswNzEyNAovLy8KTUdCMDMyCVF1aW5vbGluaWMgYWNpZCBzeW50aGVzaXMKSzAzNTE3LEswMDI3OAlUSUdSMDA1NTAsVElHUjAwNTUxCWJhY3ROT0cwMjIzOSxDT0cwMDI5CUswMzUxNyxUSUdSMDA1NTEJSzAzNTE3LENPRzAwMjkJVElHUjAwNTUwLEswMDI3OAlUSUdSMDA1NTAsQ09HMDAyOQliYWN0Tk9HMDIyMzksSzAwMjc4CWJhY3ROT0cwMjIzOSxUSUdSMDA1NTEKLy8vCk1HQjAzMwlRdWlub2xpbmljIGFjaWQgZGVncmFkYXRpb24KSzAwNzY3CVRJR1IwMDA3OAliYWN0Tk9HMDA0NDEKSzAwOTY5CWJhY3ROT0czMDQ2NAliYWN0Tk9HNjkzMzUKSzAxOTE2CUswMTk1MAlUSUdSMDA1NTIJQ09HMDE3MQovLy8KTUdCMDM0CUlzb3ZhbGVyaWMgYWNpZCBzeW50aGVzaXMgSSAoS0FESCBwYXRod2F5KQpLMDA4MjYJVElHUjAxMTIyCVRJR1IwMTEyMwlDT0cwMTE1CkswMDE2NixLMDAxNjcJSzExMzgxCkswOTY5OQpLMDAzODIJVElHUjAxMzUwCkswMTA3MwlLMTA4MDYJYmFjdE5PRzI3MDQzCi8vLwpNR0IwMzUJSXNvdmFsZXJpYyBhY2lkIHN5bnRoZXNpcyBJSSAoS0FEQyBwYXRod2F5KQpLMDA4MjYJVElHUjAxMTIyCVRJR1IwMTEyMwlDT0cwMTE1CkswMTU2OAlDT0czOTYxCkswMDEyOAliYWN0Tk9HMDA4NTIKLy8vCk1HQjAzNglTLUFkZW5vc3lsbWV0aGlvbmluZSAoU0FNKSBzeW50aGVzaXMKSzAwNzg5CVRJR1IwMTAzNAliYWN0Tk9HMDI5MzcKLy8vCk1HQjAzNwlJbm9zaXRvbCBzeW50aGVzaXMKSzAxODU4CWJhY3ROT0cwNjgxMQpLMDEwOTIJSzE4NjQ5CUNPRzA0ODMKLy8vCk1HQjAzOAlJbm9zaXRvbCBkZWdyYWRhdGlvbgpLMDAwMTAJVElHUjA0MzgwCWJhY3ROT0cwNzE1MwovLy8KTUdCMDM5CWctSHlkcm94eWJ1dHlyaWMgYWNpZCAoR0hCKSBkZWdyYWRhdGlvbgpLMDAwNDMJSzA4MzE4CUsxODEyMAliYWN0Tk9HMTA1NjUKLy8vCk1HQjA0MAlNZW5hcXVpbm9uZSBzeW50aGVzaXMgKHZpdGFtaW4gSzIpIEkKSzAxODUxCUswMjM2MQlLMDI1NTIJSzA0NzgxCkswMjU1MQpLMDg2ODAKSzAyNTQ5CkswMTkxMQpLMDE2NjEKSzEyMDczCUsxOTIyMgpLMDI1MjMJSzAwODA1CkswMjU0OApLMDMxODMKLy8vCk1HQjA0MQlNZW5hcXVpbm9uZSBzeW50aGVzaXMgKHZpdGFtaW4gSzIpIElJIChhbHRlcm5hdGl2ZSBwYXRod2F5OiBmdXRhbG9zaW5lIHBhdGh3YXkpCksxMTc4MgpLMTE3ODMKSzExNzg0CksxMTc4NQpLMDI1NDgKSzAzMTgzCi8vLwpNR0IwNDIJUFVGQXMgc3ludGhlc2lzIChBQSwgRVBBLCBESEEpCmJhY3ROT0cwNDU3OCxiYWN0Tk9HNjIyNDUsYmFjdE5PRzAyMTk0LGJhY3ROT0c0ODMwNwlDT0czMzIxLGJhY3ROT0cwMjE5NCxiYWN0Tk9HNDgzMDcKLy8vCk1HQjA0MwlBY2V0YXRlIHN5bnRoZXNpcyBJCkswMDY1NglLMDAxNjksSzAwMTcwLEswMDE3MSxLMDAxNzIJSzAzNzM3CkswMDYyNQlLMTM3ODgJSzE1MDI0CkswMDkyNQovLy8KTUdCMDQ0CUFjZXRhdGUgc3ludGhlc2lzIElJCkswMDY1NglLMDAxNjksSzAwMTcwLEswMDE3MSxLMDAxNzIJSzAzNzM3CkswMTkwNQovLy8KTUdCMDQ1CUFjZXRhdGUgc3ludGhlc2lzIElJSQpLMDA2NTYJSzAwMTY5LEswMDE3MCxLMDAxNzEsSzAwMTcyCUswMzczNwpLMTgxMTgKLy8vCk1HQjA0NglBY2V0YXRlIHN5bnRoZXNpcyBJVgpLMDE1NjgKSzE0NTE5Ci8vLwpNR0IwNDcJQWNldGF0ZSBkZWdyYWRhdGlvbgpLMDE4OTUJYmFjdE5PRzAwOTQ3CWJhY3ROT0cwMTI1OQovLy8KTUdCMDQ4CVByb3Bpb25hdGUgc3ludGhlc2lzIEkKSzEzOTIzCkswMDkzMglLMTk2OTcKLy8vCk1HQjA0OQlUcnlwdG9waGFuIGRlZ3JhZGF0aW9uCkswMTY2NwlUSUdSMDI2MTcJYmFjdE5PRzAyODgxCi8vLwpNR0IwNTAJR2x1dGFtYXRlIGRlZ3JhZGF0aW9uIEkKSzAwMjYwCUsxNTM3MQlLMDAyNjEKSzAwMTA5CkswMTAzOSxLMDEwNDAKSzAxNjE1Ci8vLwpNR0IwNTEJR2x1dGFtYXRlIGRlZ3JhZGF0aW9uIElJCkswMTg0NixLMTkyNjgKSzA0ODM1Ci8vLwpNR0IwNTIJQnV0eXJhdGUgc3ludGhlc2lzIEkKSzAwMjQ4CkswMDYzNApLMDA5MjkKLy8vCk1HQjA1MwlCdXR5cmF0ZSBzeW50aGVzaXMgSUkKSzAwMjQ4CkswMTAzNCxLMDEwMzUJSzE5NzA5Ci8vLwpNR0IwNTQJUHJvcGlvbmF0ZSBzeW50aGVzaXMgSUkKSzAxMDI2Ci8vLwpNR0IwNTUJUHJvcGlvbmF0ZSBzeW50aGVzaXMgSUlJCkswMTg0NwlLMDE4NDksSzAxODQ4CUsxMTk0MgpLMDU2MDYKSzAzNDE2LEsxNzQ4OSxLMTc0OTAKLy8vCk1HQjA1NglQcm9waW9uYXRlIGRlZ3JhZGF0aW9uIEkKSzAxOTA4CkswMTY1OQpLMDE3MjAKSzAxNjgyCkswMzQxNwo="
)

def _load_gbm_text():
    embedded_bytes = base64.b64decode(_GBM_EMBEDDED_B64)
    embedded_sha = hashlib.sha256(embedded_bytes).hexdigest()

    downloaded_bytes = None
    try:
        resp = requests.get(CONFIG["gbm_url"], timeout=15)
        resp.raise_for_status()
        downloaded_bytes = resp.content
    except Exception as e:
        print(f"[info] Could not download GBMs.v1.0.txt ({e}); using embedded offline copy instead.")

    if downloaded_bytes is not None:
        downloaded_sha = hashlib.sha256(downloaded_bytes).hexdigest()
        if downloaded_sha == embedded_sha:
            print("[ok] Downloaded GBMs.v1.0.txt is byte-identical to the pinned embedded copy.")
            print(f"     SHA-256: {downloaded_sha}")
            return downloaded_bytes.decode("utf-8"), downloaded_sha, "downloaded (verified)"
        else:
            print("[warning] Downloaded file does NOT match the pinned SHA-256 — falling back to the "
                  "embedded, verified copy to guarantee reproducibility. (This can happen if the "
                  "upstream repository has since changed; the embedded copy is still the exact file "
                  "this notebook was validated against.)")

    print(f"[ok] Using embedded GBMs.v1.0.txt. SHA-256: {embedded_sha}")
    return embedded_bytes.decode("utf-8"), embedded_sha, "embedded fallback"

GBM_RAW_TEXT, GBM_SHA256, GBM_SOURCE_MODE = _load_gbm_text()
CONFIG["gbm_sha256_expected"] = GBM_SHA256
print("GBM source used for this run:", GBM_SOURCE_MODE)


CELL4

# ============================================================
# 4. PARSE THE OFFICIAL GBM STEP LOGIC (AND across steps, OR within a step,
#    AND within a comma-joined complex)
# ============================================================
def parse_gbm_db(raw_text):
        # Parse GBMs.v1.0.txt into {module_id: {'name': str, 'steps': [[[ortholog,...], ...], ...]}}.
    #
    # File format (confirmed against the Nature Microbiology 2019 supplementary legend):
    # - each module block is terminated by a line containing only '///'
    # - the first line of a block is 'MGBxxx<TAB>Module name'
    # - every subsequent line = one obligatory functional STEP (steps are implicitly AND-ed)
    # - within a step-line, TAB-separated fields = alternative orthologs that can independently
    # satisfy the step (OR relationship)
    # - within one tab-separated field, COMMA-separated ortholog IDs are subunits of one enzyme
    # complex and are jointly required (AND relationship)
    modules = OrderedDict()
    blocks = [b for b in raw_text.strip().split("///") if b.strip()]
    for b in blocks:
        lines = [l for l in b.strip("\n").split("\n") if l.strip() != ""]
        if not lines:
            continue
        header = lines[0].split("\t")
        mod_id = header[0].strip()
        mod_name = header[1].strip() if len(header) > 1 else ""
        steps = []
        for line in lines[1:]:
            alt_fields = line.split("\t")
            alt_groups = []
            for field in alt_fields:
                members = [m.strip() for m in field.split(",") if m.strip()]
                if members:
                    alt_groups.append(members)
            if alt_groups:
                steps.append(alt_groups)
        if mod_id:
            modules[mod_id] = {"name": mod_name, "steps": steps}
    return modules

GBM = parse_gbm_db(GBM_RAW_TEXT)

# --- sanity checks: fail loudly and immediately if the parse looks wrong, rather than producing a
#     silently-broken downstream analysis ---
assert len(GBM) == 56, f"Expected 56 GBM modules, parsed {len(GBM)} — check GBM source file integrity."
assert "MGB005" in GBM and GBM["MGB005"]["name"] == "Tryptophan synthesis", "MGB005 sanity check failed."
assert len(GBM["MGB005"]["steps"]) == 5, "MGB005 should have 5 required steps."
assert GBM["MGB022"]["steps"][0][0] == ["K01580"], "MGB022 (GABA synthesis III) step-1 KO check failed."
print(f"[ok] Parsed {len(GBM)} GBM modules from a verified source ({GBM_SOURCE_MODE}).")
print(f"[ok] Sanity checks passed (MGB005 = Tryptophan synthesis, 5 steps; MGB022 step-1 = K01580).")

# Build a flat KO -> [(module_id, step_num, alt_group_num)] index, and the set of ALL KOs referenced
# anywhere in the GBM database (used later to report which KOs in your data are GBM-relevant at all).
KO_INDEX = defaultdict(list)
ALL_GBM_KOS = set()
for mid, d in GBM.items():
    for s_i, step in enumerate(d["steps"], 1):
        for g_i, grp in enumerate(step, 1):
            for member in grp:
                if re.fullmatch(r"K\d{5}", member):
                    KO_INDEX[member].append((mid, s_i, g_i))
                    ALL_GBM_KOS.add(member)
print(f"[ok] {len(ALL_GBM_KOS)} distinct KEGG Orthology (KO) identifiers are used across all 56 modules.")


CELL5

# ============================================================
# 5. UPLOAD / LOAD YOUR EGGNOG-MAPPER EXCEL FILE
# ============================================================
if IN_COLAB:
    from google.colab import files
    print("Select your eggNOG-mapper .xlsx export...")
    uploaded = files.upload()
    assert len(uploaded) >= 1, "No file uploaded."
    EGGNOG_PATH = list(uploaded.keys())[0]
else:
    # Non-Colab (local Jupyter) fallback: edit this path directly.
    EGGNOG_PATH = "eggnog_annotations.xlsx"

print("Using input file:", EGGNOG_PATH)
xl = pd.ExcelFile(EGGNOG_PATH)
print("Sheets found:", xl.sheet_names)


CELL6

# ============================================================
# 6. EXTRACT KO PRESENCE PER GENOME
# ============================================================
def _find_ko_column(df):
    lower_map = {c.lower(): c for c in df.columns}
    for cand in CONFIG["ko_column_candidates"]:
        if cand.lower() in lower_map:
            return lower_map[cand.lower()]
    return None

def _find_genome_column(df):
    lower_map = {c.lower(): c for c in df.columns}
    for cand in CONFIG["genome_id_candidates"]:
        if cand.lower() in lower_map:
            return lower_map[cand.lower()]
    return None

def _parse_ko_cell(val):
        # eggNOG-mapper KEGG_ko cell values look like 'ko:K00001,ko:K00002' or '-' or NaN.
    # Returns a (possibly empty) set of bare KO IDs (e.g. {'K00001','K00002'}).
    if pd.isna(val):
        return set()
    s = str(val).strip()
    if s in ("", "-", "nan", "NaN", "NA"):
        return set()
    kos = set()
    for tok in re.split(r"[;,\s]+", s):
        tok = tok.strip()
        if not tok:
            continue
        tok = tok.replace("ko:", "").replace("KO:", "")
        if re.fullmatch(r"K\d{5}", tok):
            kos.add(tok)
    return kos

def _clean_emapper_df(raw):
        # Drop eggNOG-mapper '##' comment rows if the sheet was exported with them still in place,
    # and use row 0 as header if the real header wasn't picked up as column names.
    df = raw.copy()
    # If the KO column isn't found, but a '#query' row appears inside the data (common when the
    # header got exported as data), try re-reading with the first non-comment row as header.
    if _find_ko_column(df) is None:
        first_col = df.columns[0]
        mask_comment = df[first_col].astype(str).str.startswith("##")
        df = df[~mask_comment]
        header_row_mask = df[first_col].astype(str).str.contains("query", case=False, na=False)
        if header_row_mask.any():
            header_idx = df.index[header_row_mask][0]
            new_header = df.loc[header_idx].tolist()
            df = df.loc[header_idx + 1:].copy()
            df.columns = new_header
    else:
        first_col = df.columns[0]
        mask_comment = df[first_col].astype(str).str.startswith("##")
        df = df[~mask_comment]
    return df.reset_index(drop=True)

genome_ko_sets = OrderedDict()   # genome_id -> set of KOs present
genome_gene_ko_rows = []         # audit trail: one row per (genome, query gene, KO)

mode = CONFIG["genome_split_mode"]
if mode == "auto":
    mode = "per_sheet" if len(xl.sheet_names) > 1 else "single_sheet_autodetect"

if mode == "per_sheet":
    for sheet in xl.sheet_names:
        raw = xl.parse(sheet, header=0)
        df = _clean_emapper_df(raw)
        ko_col = _find_ko_column(df)
        assert ko_col is not None, (
            f"Could not find a KEGG_ko-type column in sheet '{sheet}'. "
            f"Columns present: {list(df.columns)}. "
            f"Add the correct column name to CONFIG['ko_column_candidates'] and re-run from cell 6."
        )
        query_col = df.columns[0]
        kos_this_genome = set()
        for _, row in df.iterrows():
            kos = _parse_ko_cell(row[ko_col])
            kos_this_genome |= kos
            for ko in kos:
                genome_gene_ko_rows.append({"Genome": sheet, "Query_gene": row[query_col], "KO_ID": ko})
        genome_ko_sets[sheet] = kos_this_genome

elif mode == "single_genome":
    raw = xl.parse(xl.sheet_names[0], header=0)
    df = _clean_emapper_df(raw)
    ko_col = _find_ko_column(df)
    assert ko_col is not None, (
        f"Could not find a KEGG_ko-type column. Columns present: {list(df.columns)}. "
        f"Add the correct column name to CONFIG['ko_column_candidates'] and re-run from cell 6."
    )
    query_col = df.columns[0]
    genome_name = EGGNOG_PATH.rsplit("/", 1)[-1].rsplit(".", 1)[0]
    kos_this_genome = set()
    for _, row in df.iterrows():
        kos = _parse_ko_cell(row[ko_col])
        kos_this_genome |= kos
        for ko in kos:
            genome_gene_ko_rows.append({"Genome": genome_name, "Query_gene": row[query_col], "KO_ID": ko})
    genome_ko_sets[genome_name] = kos_this_genome

else:  # single_sheet_autodetect
    raw = xl.parse(xl.sheet_names[0], header=0)
    df = _clean_emapper_df(raw)
    ko_col = _find_ko_column(df)
    assert ko_col is not None, (
        f"Could not find a KEGG_ko-type column. Columns present: {list(df.columns)}. "
        f"Add the correct column name to CONFIG['ko_column_candidates'] and re-run from cell 6."
    )
    genome_col = _find_genome_column(df)
    query_col = df.columns[0]
    if genome_col is not None:
        for genome_id, sub in df.groupby(genome_col):
            kos_this_genome = set()
            for _, row in sub.iterrows():
                kos = _parse_ko_cell(row[ko_col])
                kos_this_genome |= kos
                for ko in kos:
                    genome_gene_ko_rows.append({"Genome": str(genome_id), "Query_gene": row[query_col], "KO_ID": ko})
            genome_ko_sets[str(genome_id)] = kos_this_genome
    else:
        genome_name = EGGNOG_PATH.rsplit("/", 1)[-1].rsplit(".", 1)[0]
        kos_this_genome = set()
        for _, row in df.iterrows():
            kos = _parse_ko_cell(row[ko_col])
            kos_this_genome |= kos
            for ko in kos:
                genome_gene_ko_rows.append({"Genome": genome_name, "Query_gene": row[query_col], "KO_ID": ko})
        genome_ko_sets[genome_name] = kos_this_genome
        print("[info] No sheet/genome-ID column found; the whole input was treated as ONE genome: "
              f"'{genome_name}'. If your workbook actually contains multiple genomes in one sheet, "
              "add a genome/sample-ID column and re-run, or set CONFIG['genome_split_mode']='per_sheet' "
              "with one sheet per genome instead.")

assert len(genome_ko_sets) > 0, "No genomes were parsed from the input file — check CONFIG and file layout."
for g, kos in genome_ko_sets.items():
    print(f"  {g}: {len(kos)} distinct KOs annotated, {len(kos & ALL_GBM_KOS)} of which are GBM-relevant")
GENOME_IDS = list(genome_ko_sets.keys())
print(f"\n[ok] {len(GENOME_IDS)} genome(s) loaded.")


CELL7

# ============================================================
# 7. RESOLVE AND/OR MODULE COMPLETENESS FOR EVERY GENOME AND EVERY MODULE
# ============================================================
def ko_present(genome_id, ko):
    return ko in genome_ko_sets[genome_id]

def group_satisfied(genome_id, alt_group):
        # An alternative-group is satisfied only if EVERY ortholog listed in it is present
    # (AND within a group = enzyme-complex subunits). Non-KO orthologs (COG/NOG/TIGR) cannot be
    # evaluated from eggNOG KEGG_ko data and are treated as 'not confirmed' (conservative; never
    # fabricated as present).
    for member in alt_group:
        if not re.fullmatch(r"K\d{5}", member):
            return False  # COG/NOG/TIGR-only alternative: cannot be confirmed from KO data alone
        if not ko_present(genome_id, member):
            return False
    return True

def step_present(genome_id, step):
        # A step is present if ANY of its OR-alternative groups is fully satisfied.
    return any(group_satisfied(genome_id, grp) for grp in step)

step_rows = []
module_rows = []
missing_ko_rows = []

for genome_id in GENOME_IDS:
    for mid, d in GBM.items():
        steps = d["steps"]
        n_steps = len(steps)
        n_present = 0
        for s_i, step in enumerate(steps, 1):
            present = step_present(genome_id, step)
            n_present += int(present)
            step_rows.append({
                "Genome": genome_id, "Module_ID": mid, "Module_name": d["name"],
                "Functional_step": s_i, "Step_present": int(present),
            })
            if not present:
                # record exactly which KOs would have satisfied this step, for troubleshooting
                for g_i, grp in enumerate(step, 1):
                    ko_members = [m for m in grp if re.fullmatch(r"K\d{5}", m)]
                    non_ko_members = [m for m in grp if not re.fullmatch(r"K\d{5}", m)]
                    missing_ko_rows.append({
                        "Genome": genome_id, "Module_ID": mid, "Module_name": d["name"],
                        "Functional_step": s_i, "Alternative_group": g_i,
                        "Required_KOs_for_this_alternative": ",".join(ko_members) if ko_members else "(none)",
                        "Required_non_KO_orthologs (uncheckable from KEGG_ko)": ",".join(non_ko_members) if non_ko_members else "(none)",
                        "KOs_present_in_this_genome": ",".join(sorted(k for k in ko_members if ko_present(genome_id, k))) or "(none)",
                        "KOs_missing_in_this_genome": ",".join(sorted(k for k in ko_members if not ko_present(genome_id, k))) or "(none)",
                    })
        completeness_pct = 100.0 * n_present / n_steps if n_steps > 0 else np.nan
        module_rows.append({
            "Genome": genome_id, "Module_ID": mid, "Module_name": d["name"],
            "N_required_steps": n_steps, "N_steps_present": n_present,
            "Completeness_%": round(completeness_pct, 2),
            "Status": "COMPLETE" if n_present == n_steps else "INCOMPLETE",
        })

steps_df = pd.DataFrame(step_rows)
module_df = pd.DataFrame(module_rows)
missing_df = pd.DataFrame(missing_ko_rows)

print(f"[ok] Computed completeness for {module_df['Module_ID'].nunique()} modules x {len(GENOME_IDS)} genome(s).")
print(module_df.groupby("Status").size())


CELL8

# ============================================================
# 8. SUMMARY TABLES
# ============================================================
thr = CONFIG["completeness_analytical_threshold_pct"]

module_summary = (
    module_df.groupby(["Module_ID", "Module_name", "N_required_steps"])
    .agg(
        N_genomes=("Genome", "nunique"),
        Mean_completeness_pct=("Completeness_%", "mean"),
        Median_completeness_pct=("Completeness_%", "median"),
        SD_pct=("Completeness_%", "std"),
        Min_pct=("Completeness_%", "min"),
        Max_pct=("Completeness_%", "max"),
        N_genomes_complete=("Status", lambda s: (s == "COMPLETE").sum()),
    )
    .reset_index()
)
module_summary["Pct_genomes_complete"] = (
    100.0 * module_summary["N_genomes_complete"] / module_summary["N_genomes"]
).round(2)
module_summary[f"N_genomes_ge{thr}pct (analytical threshold, not official GBM cutoff)"] = module_df.groupby(
    "Module_ID")["Completeness_%"].apply(lambda s: (s >= thr).sum()).values
for c in ["Mean_completeness_pct", "Median_completeness_pct", "SD_pct", "Min_pct", "Max_pct"]:
    module_summary[c] = module_summary[c].round(2)
module_summary = module_summary.sort_values(["Pct_genomes_complete", "Module_ID"], ascending=[False, True])

genome_profile = module_df.pivot(index="Genome", columns="Module_ID", values="Completeness_%")
genome_profile = genome_profile.loc[GENOME_IDS]

n_complete_modules = module_df[module_df.Status == "COMPLETE"]["Module_ID"].nunique()
n_total_modules = module_df["Module_ID"].nunique()
print(f"Across all genome(s), {n_complete_modules}/{n_total_modules} modules were COMPLETE in at least one genome.")
module_summary.head(15)


CELL9

# ============================================================
# 9. WRITE THE FINAL EXCEL REPORT
# ============================================================
gbm_def_rows = []
for mid, d in GBM.items():
    for s_i, step in enumerate(d["steps"], 1):
        for g_i, grp in enumerate(step, 1):
            gbm_def_rows.append({
                "Module_ID": mid, "Module_name": d["name"], "Functional_step": s_i,
                "Alternative_group (OR)": g_i,
                "Orthologs_required_in_group (AND if >1)": " + ".join(grp),
                "KO-type_members_only": " + ".join([m for m in grp if re.fullmatch(r"K\d{5}", m)]) or "(none - COG/NOG/TIGR only)",
            })
gbm_def_df = pd.DataFrame(gbm_def_rows)

ko_gene_df = pd.DataFrame(genome_gene_ko_rows) if genome_gene_ko_rows else pd.DataFrame(
    columns=["Genome", "Query_gene", "KO_ID"])

readme_lines = [
    ("GBM Module Completeness from eggNOG-mapper Annotations", True, 14),
    ("", False, 10),
    (f"Input file: {EGGNOG_PATH}", False, 10),
    (f"Genomes analysed: {len(GENOME_IDS)} ({', '.join(GENOME_IDS[:10])}{' ...' if len(GENOME_IDS) > 10 else ''})", False, 10),
    (f"GBM database: GBMs.v1.0 (Valles-Colomer et al., Nat Microbiol 4:623-632, 2019)", False, 10),
    (f"GBM source: {CONFIG['gbm_url']}", False, 10),
    (f"GBM source mode this run: {GBM_SOURCE_MODE}  |  SHA-256: {GBM_SHA256}", False, 10),
    ("", False, 10),
    ("METHOD", True, 12),
    ("For every module, a functional step is present if ANY official OR-alternative is fully", False, 10),
    ("satisfied (ALL KOs of that alternative present; comma-joined KOs = obligate enzyme-complex", False, 10),
    ("subunits, AND-required). A module is COMPLETE only if ALL of its required steps are present.", False, 10),
    ("Completeness % = (steps present) / (total required steps) x 100 - never a raw gene-count ratio.", False, 10),
    ("Alternatives whose ONLY orthologs are COG/NOG/TIGR (not KO) cannot be evaluated from eggNOG's", False, 10),
    ("KEGG_ko column and are conservatively treated as not confirmed (never assumed present).", False, 10),
    (f"The {thr}% threshold column is an analytical convenience, NOT an official GBM cutoff.", False, 10),
    ("", False, 10),
    ("SHEETS", True, 12),
    ("- GBM_Definitions: official step-level AND/OR logic for all 56 modules.", False, 10),
    ("- KO_Gene_Evidence: every (genome, query gene, KO) triple extracted from your eggNOG file.", False, 10),
    ("- Genome_Functional_Steps: per-genome, per-step present/absent (logic-resolved).", False, 10),
    ("- Module_Completeness_Summary: descriptive statistics per module across your genome(s).", False, 10),
    ("- Genome_Module_Profile: per-genome completeness %, one row per genome, wide format.", False, 10),
    ("- Missing_KOs_detail: for every incomplete step, exactly which KO(s) are missing.", False, 10),
    ("", False, 10),
    ("REPRODUCIBILITY", True, 12),
    ("This entire analysis is deterministic given the same input file and the same GBM database", False, 10),
    ("(verified by SHA-256 above). No thresholds were tuned to change the result; the only", False, 10),
    (f"user-editable analytical convenience is the {thr}% 'near-complete' label, clearly marked as such.", False, 10),
]

with pd.ExcelWriter(CONFIG["output_xlsx"], engine="xlsxwriter") as writer:
    wb = writer.book
    header_fmt = wb.add_format({"bold": True, "bg_color": "#1F4E78", "font_color": "white",
                                 "text_wrap": True, "valign": "vcenter", "border": 1})
    title_fmt = wb.add_format({"bold": True, "font_size": 14})
    sub_fmt = wb.add_format({"bold": True, "font_size": 12})
    normal_fmt = wb.add_format({"font_size": 10})

    ws = wb.add_worksheet("README")
    writer.sheets["README"] = ws
    ws.set_column(0, 0, 110)
    for r, (text, bold, size) in enumerate(readme_lines):
        fmt = title_fmt if (bold and size >= 14) else (sub_fmt if bold else normal_fmt)
        ws.write(r, 0, text, fmt)

    def dump(df, sheet_name, note=None):
        start = 0
        if note:
            wsx = wb.add_worksheet(sheet_name)
            writer.sheets[sheet_name] = wsx
            wsx.write(0, 0, note, wb.add_format({"italic": True, "font_size": 9}))
            start = 2
            df.to_excel(writer, sheet_name=sheet_name, startrow=start, index=False)
        else:
            df.to_excel(writer, sheet_name=sheet_name, index=False)
        wsx = writer.sheets[sheet_name]
        for c_i, col in enumerate(df.columns):
            wsx.write(start, c_i, col, header_fmt)
            try:
                width = min(60, max(12, df[col].astype(str).str.len().quantile(0.9) + 2))
            except Exception:
                width = 15
            wsx.set_column(c_i, c_i, width)
        wsx.freeze_panes(start + 1, 0)

    dump(gbm_def_df, "GBM_Definitions",
         "Official AND/OR step logic for all 56 GBMs.v1.0 modules (source verified above).")
    dump(ko_gene_df, "KO_Gene_Evidence",
         "Every (genome, query gene, KO) triple extracted from the KEGG_ko column of your eggNOG file.")
    dump(steps_df, "Genome_Functional_Steps",
         "Per-genome, per-step present(1)/absent(0), resolved using the official AND/OR logic.")
    dump(module_summary, "Module_Completeness_Summary",
         "Descriptive statistics of module completeness across your genome(s).")
    gp_out = genome_profile.reset_index()
    dump(gp_out, "Genome_Module_Profile",
         "Per-genome completeness (%) for every module; one row per genome.")
    dump(missing_df, "Missing_KOs_detail",
         "For every step NOT satisfied, the exact KO(s) required and which are missing in that genome.")

print(f"[ok] Workbook written: {CONFIG['output_xlsx']}")


CELL 10

# ============================================================
# 10. DOWNLOAD THE RESULT
# ============================================================
if IN_COLAB:
    from google.colab import files
    files.download(CONFIG["output_xlsx"])
else:
    print(f"Output written to: {CONFIG['output_xlsx']} (open it directly).")
