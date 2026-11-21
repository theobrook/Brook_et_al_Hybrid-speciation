#!/usr/bin/env python3
"""
009_fst_only_summary.py
=======================
Compiles and summarises FST-only outlier genes across all 15 pairwise
comparisons from the 008_tair_annotation outputs.

Focus metrics:
  FST_top5   — top 5% empirical FST outliers (primary)
  FST_top10  — top 10% empirical FST outliers (supplementary)

Outputs (all in 009_fst_only_summary/):
  master_table.tsv
      One row per gene × comparison × metric.

  gene_comparison_matrix.tsv
      Wide matrix: genes as rows, comparisons as columns.
      Values: "top5", "top10", "both", or "-".

  gene_annotation_lookup.tsv
      De-duplicated annotation table: one row per gene.

  universal_outliers.tsv
      Genes that are FST outliers in every pairwise comparison.

  exclusive_outliers.tsv
      Genes found in only one pairwise comparison.

  comparison_counts.tsv
      Per-comparison gene counts for top5 and top10.

  summary_stats.txt
      Full text report.
"""

import os
import csv
import re
from collections import defaultdict

# ── Paths ──────────────────────────────────────────────────────────────────

base_dir  = "results/outliers"
input_dir = os.path.join(base_dir, "008_tair_annotation", "pairwise_outliers")
output_dir = os.path.join(base_dir, "009_fst_only_summary")
os.makedirs(output_dir, exist_ok=True)

FOCUS_METRICS = {
    "FST_top5":  "top5",
    "FST_top10": "top10",
}

ANNOTATION_COLS = [
    "gene_id", "preferred_name",
    "blast_arab_id", "blast_pident", "blast_qcov", "blast_evalue",
    "blast_arab_description",
    "tair_symbol_aliases", "tair_short_description", "tair_curator_summary",
    "tair_phenotypes",
    "eggnog_description", "seed_ortholog",
    "GO_terms", "KEGG_ko", "KEGG_Pathway", "COG_category",
    "adaptation_flag", "adaptation_keywords",
    "reproduction_flag", "reproduction_keywords",
]

MASTER_COLS = ["comparison", "metric", "threshold_class"] + ANNOTATION_COLS

# ── Step 1: Load all data ──────────────────────────────────────────────────

print("Scanning 008_tair_annotation output directories...")

master_rows      = []
annotation_store = {}                      # gene_id -> best annotation dict
gene_occurrences = defaultdict(set)        # gene_id -> set of (comparison, metric)

for pair_name in sorted(os.listdir(input_dir)):
    pair_dir = os.path.join(input_dir, pair_name)
    if not os.path.isdir(pair_dir):
        continue

    for metric, threshold_class in FOCUS_METRICS.items():
        tsv = os.path.join(pair_dir, metric,
                           "all_outlier_genes_tair_annotated.tsv")
        if not os.path.exists(tsv):
            continue

        with open(tsv, newline="", encoding="utf-8") as f:
            rows = list(csv.DictReader(f, delimiter="\t"))

        for row in rows:
            gid = row["gene_id"]
            gene_occurrences[gid].add((pair_name, metric))

            # Keep best annotation (prefer rows with BLAST hit)
            if gid not in annotation_store or (
                row.get("blast_arab_id") and
                not annotation_store[gid].get("blast_arab_id")
            ):
                annotation_store[gid] = {c: row.get(c, "")
                                         for c in ANNOTATION_COLS}

            master_row = {
                "comparison":       pair_name,
                "metric":           metric,
                "threshold_class":  threshold_class,
            }
            master_row.update({c: row.get(c, "") for c in ANNOTATION_COLS})
            master_rows.append(master_row)

        print(f"  {pair_name} - {metric}: {len(rows)} genes")

print(f"\nTotal gene-comparison-metric records: {len(master_rows):,}")
print(f"Unique genes: {len(annotation_store):,}")

# ── Step 2: Re-run keyword flagging with full TAIR text ────────────────────

ADAPTATION_KEYWORDS = [
    ("abiotic", False), ("temperature stimulus", False),
    ("cold acclimation", False), ("cold stress", False), ("cold response", False),
    ("heat stress", False), ("heat shock", False),
    ("drought", False), ("water deprivation", False), ("desiccation", False),
    ("salt stress", False), ("salt tolerance", False), ("osmotic stress", False),
    ("oxidative stress", False), ("hypoxia", False), ("freezing", False),
    ("wounding", False), ("defense response", False), ("pathogen", False),
    ("immune", False), ("UV", True), ("ultraviolet", False),
    ("light stimulus", False), ("light response", False),
    ("photoperiod", False), ("circadian", False),
    ("stomata", False), ("transpiration", False),
    ("photosynthesis", False), ("chloroplast", False),
    ("abscisic acid", False), ("jasmonic", False), ("ethylene", False),
    ("auxin", False), ("reactive oxygen", False), ("antioxidant", False),
    ("dehydration", False), ("dehydrin", False),
    ("CBF", True), ("DREB", True), ("LEA protein", False),
    ("late embryogenesis abundant", False),
    ("HSP", True), ("HKT", True), ("NHX", True), ("SOS1", True),
    ("ICE1", True), ("CCA1", True), ("TOC1", True), ("LHY", True),
    ("iron-sulfur", False), ("metal cluster", False), ("metal binding", False),
    ("tocopherol", False), ("carotenoid", False),
    ("superoxide dismutase", False), ("peroxidase", False), ("catalase", False),
    ("aquaporin", False), ("water channel", False),
    ("cold-regulated", False), ("cold-responsive", False),
    ("senescence", False), ("stress response", False), ("stress tolerance", False),
    ("abiotic stress", False), ("biotic stress", False),
]

REPRO_KEYWORDS = [
    ("flower development", False), ("floral development", False),
    ("floral organ", False), ("floral meristem", False),
    ("pollen development", False), ("pollen tube", False),
    ("pollen germination", False), ("pollen wall", False),
    ("pistil", False), ("ovule", False),
    ("seed development", False), ("embryo development", False),
    ("meiosis", False), ("meiotic", False),
    ("gametophyte development", False), ("megagametogenesis", False),
    ("microsporogenesis", False), ("tapetum", False),
    ("self-incompatibility", False),
    ("flowering time", False), ("flowering locus", False),
    ("photoperiodism", False), ("vernalization", False),
    ("pollination", False),
    ("FLC", True), ("CONSTANS", False), ("SOC1", True),
    ("GIGANTEA", False), ("PHYTOCHROME", False), ("CRYPTOCHROME", False),
    ("FRIGIDA", False), ("AGAMOUS", False), ("APETALA", False),
    ("PISTILLATA", False), ("SEPALLATA", False),
    ("inflorescence", False),
    ("anther development", False), ("stamen development", False),
    ("carpel development", False),
    ("male gametophyte", False), ("female gametophyte", False),
]

def flag_keywords(text, kwlist):
    if not text: return []
    return [kw for kw, wb in kwlist if
            (re.search(r'\b' + re.escape(kw) + r'\b', text, re.I) if wb
             else kw.lower() in text.lower())]

print("\nRe-running keyword flagging with full TAIR text...")
n_adapt_new = n_repro_new = 0

for gid, ann in annotation_store.items():
    text = " ".join(filter(None, [
        ann.get("blast_arab_description", ""),
        ann.get("tair_short_description", ""),
        ann.get("tair_curator_summary", ""),
        ann.get("tair_phenotypes", ""),
        ann.get("eggnog_description", ""),
        ann.get("GO_terms", ""),
        ann.get("preferred_name", ""),
        ann.get("tair_symbol_aliases", ""),
    ]))
    af = flag_keywords(text, ADAPTATION_KEYWORDS)
    rf = flag_keywords(text, REPRO_KEYWORDS)
    old_a = ann.get("adaptation_flag", "")
    old_r = ann.get("reproduction_flag", "")
    ann["adaptation_flag"]     = "YES" if af else ""
    ann["adaptation_keywords"] = "; ".join(af) if af else ""
    ann["reproduction_flag"]     = "YES" if rf else ""
    ann["reproduction_keywords"] = "; ".join(rf) if rf else ""
    if af and not old_a: n_adapt_new += 1
    if rf and not old_r: n_repro_new += 1

# Propagate back to master rows
for row in master_rows:
    gid = row["gene_id"]
    if gid in annotation_store:
        for col in ["adaptation_flag", "adaptation_keywords",
                    "reproduction_flag", "reproduction_keywords"]:
            row[col] = annotation_store[gid].get(col, "")

total_adapt = sum(1 for a in annotation_store.values() if a.get("adaptation_flag"))
total_repro = sum(1 for a in annotation_store.values() if a.get("reproduction_flag"))
print(f"  Newly flagged — adaptation: {n_adapt_new}, reproduction: {n_repro_new}")
print(f"  Total flagged — adaptation: {total_adapt}, reproduction: {total_repro}")

# ── Step 3: Build occurrence structures ───────────────────────────────────

# metric -> comparison -> set of gene_ids
structures = {m: defaultdict(set) for m in FOCUS_METRICS}

for row in master_rows:
    structures[row["metric"]][row["comparison"]].add(row["gene_id"])

pairwise_comparisons = sorted(
    set().union(*[s.keys() for s in structures.values()])
)
n_pairwise = len(pairwise_comparisons)

# ── Step 4: Write master_table.tsv ────────────────────────────────────────

with open(os.path.join(output_dir, "master_table.tsv"),
          "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=MASTER_COLS, delimiter="\t",
                       extrasaction="ignore")
    w.writeheader(); w.writerows(master_rows)
print(f"\nWrote: master_table.tsv  ({len(master_rows):,} rows)")

# ── Step 5: Write gene_annotation_lookup.tsv ──────────────────────────────

with open(os.path.join(output_dir, "gene_annotation_lookup.tsv"),
          "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=ANNOTATION_COLS, delimiter="\t")
    w.writeheader()
    for gid in sorted(annotation_store):
        w.writerow(annotation_store[gid])
print(f"Wrote: gene_annotation_lookup.tsv  ({len(annotation_store):,} genes)")

# ── Step 6: Gene × comparison matrix ──────────────────────────────────────

all_genes = set().union(*[
    set().union(*s.values()) for s in structures.values()
])

def encode_cell(gid, comparison):
    in5  = gid in structures["FST_top5"].get(comparison, set())
    in10 = gid in structures["FST_top10"].get(comparison, set())
    if in5 and in10: return "both"
    if in5:          return "top5"
    if in10:         return "top10"
    return "-"

matrix_cols = (
    ["gene_id", "preferred_name", "blast_arab_id",
     "tair_short_description", "tair_curator_summary",
     "GO_terms", "adaptation_flag", "reproduction_flag",
     "n_comparisons_top5", "n_comparisons_top10"]
    + pairwise_comparisons
)

matrix_rows = []
for gid in sorted(all_genes):
    ann = annotation_store.get(gid, {})
    n5  = sum(1 for c in pairwise_comparisons
              if gid in structures["FST_top5"].get(c, set()))
    n10 = sum(1 for c in pairwise_comparisons
              if gid in structures["FST_top10"].get(c, set()))
    mrow = {
        "gene_id":                gid,
        "preferred_name":         ann.get("preferred_name", ""),
        "blast_arab_id":          ann.get("blast_arab_id", ""),
        "tair_short_description": ann.get("tair_short_description", ""),
        "tair_curator_summary":   ann.get("tair_curator_summary", ""),
        "GO_terms":               ann.get("GO_terms", ""),
        "adaptation_flag":        ann.get("adaptation_flag", ""),
        "reproduction_flag":      ann.get("reproduction_flag", ""),
        "n_comparisons_top5":     n5,
        "n_comparisons_top10":    n10,
    }
    for comp in pairwise_comparisons:
        mrow[comp] = encode_cell(gid, comp)
    matrix_rows.append(mrow)

# Sort by breadth (most widespread first)
matrix_rows.sort(key=lambda r: -(r["n_comparisons_top5"] +
                                  r["n_comparisons_top10"]))

with open(os.path.join(output_dir, "gene_comparison_matrix.tsv"),
          "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=matrix_cols, delimiter="\t")
    w.writeheader(); w.writerows(matrix_rows)
print(f"Wrote: gene_comparison_matrix.tsv  "
      f"({len(matrix_rows):,} genes × {n_pairwise} comparisons)")

# ── Step 7: Universal and exclusive outliers ───────────────────────────────

universal_rows = []
exclusive_rows = []

univ_cols = [
    "metric", "n_comparisons", "gene_id", "preferred_name",
    "blast_arab_id", "tair_short_description", "tair_curator_summary",
    "tair_phenotypes", "GO_terms",
    "adaptation_flag", "adaptation_keywords",
    "reproduction_flag", "reproduction_keywords", "eggnog_description",
]
excl_cols = [
    "metric", "comparison", "gene_id", "preferred_name",
    "blast_arab_id", "tair_short_description", "tair_curator_summary",
    "tair_phenotypes", "GO_terms",
    "adaptation_flag", "reproduction_flag", "eggnog_description",
]

for metric, occ in structures.items():
    comps_with_data = [c for c, gs in occ.items() if gs]
    if not comps_with_data:
        continue

    # Universal
    universal = set.intersection(*[occ[c] for c in comps_with_data])
    for gid in sorted(universal):
        ann = annotation_store.get(gid, {})
        universal_rows.append({
            "metric":             metric,
            "n_comparisons":      len(comps_with_data),
            "gene_id":            gid,
            "preferred_name":     ann.get("preferred_name", ""),
            "blast_arab_id":      ann.get("blast_arab_id", ""),
            "tair_short_description": ann.get("tair_short_description", ""),
            "tair_curator_summary":   ann.get("tair_curator_summary", ""),
            "tair_phenotypes":        ann.get("tair_phenotypes", ""),
            "GO_terms":           ann.get("GO_terms", ""),
            "adaptation_flag":    ann.get("adaptation_flag", ""),
            "adaptation_keywords": ann.get("adaptation_keywords", ""),
            "reproduction_flag":  ann.get("reproduction_flag", ""),
            "reproduction_keywords": ann.get("reproduction_keywords", ""),
            "eggnog_description": ann.get("eggnog_description", ""),
        })

    # Exclusive (in exactly one comparison)
    for comp, genes in occ.items():
        for gid in sorted(genes):
            if sum(1 for gs in occ.values() if gid in gs) == 1:
                ann = annotation_store.get(gid, {})
                exclusive_rows.append({
                    "metric":             metric,
                    "comparison":         comp,
                    "gene_id":            gid,
                    "preferred_name":     ann.get("preferred_name", ""),
                    "blast_arab_id":      ann.get("blast_arab_id", ""),
                    "tair_short_description": ann.get("tair_short_description", ""),
                    "tair_curator_summary":   ann.get("tair_curator_summary", ""),
                    "tair_phenotypes":        ann.get("tair_phenotypes", ""),
                    "GO_terms":           ann.get("GO_terms", ""),
                    "adaptation_flag":    ann.get("adaptation_flag", ""),
                    "reproduction_flag":  ann.get("reproduction_flag", ""),
                    "eggnog_description": ann.get("eggnog_description", ""),
                })

with open(os.path.join(output_dir, "universal_outliers.tsv"),
          "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=univ_cols, delimiter="\t")
    w.writeheader(); w.writerows(universal_rows)
print(f"Wrote: universal_outliers.tsv  ({len(universal_rows):,} records)")

with open(os.path.join(output_dir, "exclusive_outliers.tsv"),
          "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=excl_cols, delimiter="\t")
    w.writeheader(); w.writerows(exclusive_rows)
print(f"Wrote: exclusive_outliers.tsv  ({len(exclusive_rows):,} records)")

# ── Step 8: Comparison counts table ───────────────────────────────────────

count_rows = []
for comp in pairwise_comparisons:
    n5  = len(structures["FST_top5"].get(comp, set()))
    n10 = len(structures["FST_top10"].get(comp, set()))
    count_rows.append({
        "comparison":       comp,
        "n_top5_outliers":  n5,
        "n_top10_outliers": n10,
        "n_top10_only":     n10 - n5,   # in top10 but not top5
    })

with open(os.path.join(output_dir, "comparison_counts.tsv"),
          "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f,
        fieldnames=["comparison", "n_top5_outliers",
                    "n_top10_outliers", "n_top10_only"],
        delimiter="\t")
    w.writeheader(); w.writerows(count_rows)
print(f"Wrote: comparison_counts.tsv")

# ── Step 9: Summary report ─────────────────────────────────────────────────

sep  = "=" * 70
sep2 = "-" * 70
half = n_pairwise // 2

lines = [
    sep,
    "009 FST-ONLY OUTLIER SUMMARY REPORT",
    "Metrics: FST_top5 (primary) and FST_top10 (supplementary)",
    sep, "",
    f"Pairwise comparisons: {n_pairwise}",
    "  " + "\n  ".join(pairwise_comparisons),
    "",
]

# Gene counts by metric
lines += [sep2, "GENE COUNTS BY METRIC", sep2, ""]
for metric, occ in structures.items():
    all_g    = set().union(*occ.values()) if occ else set()
    per_comp = {c: len(gs) for c, gs in occ.items()}
    mn = min(per_comp.values()) if per_comp else 0
    mx = max(per_comp.values()) if per_comp else 0
    av = sum(per_comp.values()) / len(per_comp) if per_comp else 0
    lines += [
        f"  {metric}:",
        f"    Total unique genes (across all pairs): {len(all_g)}",
        f"    Per-comparison range: {mn}–{mx}  (mean {av:.1f})",
        "",
    ]

# Per-comparison counts
lines += [sep2, "PER-COMPARISON GENE COUNTS", sep2, ""]
header = f"{'Comparison':<40} {'top5':>6} {'top10':>6} {'top10_only':>10}"
lines += [header, "-" * len(header)]
for row in count_rows:
    lines.append(
        f"{row['comparison']:<40} "
        f"{row['n_top5_outliers']:>6} "
        f"{row['n_top10_outliers']:>6} "
        f"{row['n_top10_only']:>10}"
    )
lines.append("")

# Universal outliers
lines += [sep2, "UNIVERSAL OUTLIERS (in ALL pairwise comparisons)", sep2, ""]
for metric, occ in structures.items():
    comps_with_data = [c for c, gs in occ.items() if gs]
    universal = (set.intersection(*[occ[c] for c in comps_with_data])
                 if comps_with_data else set())
    lines.append(f"  {metric}: {len(universal)} universal genes "
                 f"(across {len(comps_with_data)} comparisons)")
    for gid in sorted(universal):
        ann  = annotation_store.get(gid, {})
        name = ann.get("preferred_name") or gid
        arab = ann.get("blast_arab_id", "")
        desc = (ann.get("tair_short_description") or
                ann.get("eggnog_description", ""))
        desc = (desc[:80] + "…") if len(desc) > 80 else desc
        lines.append(f"    {name}  [{arab}]  {desc}")
    lines.append("")

# Broadly shared outliers (>= half of comparisons)
lines += [sep2,
          f"BROADLY SHARED OUTLIERS (>= {half} of {n_pairwise} comparisons)",
          sep2, ""]
for metric, occ in structures.items():
    all_g = set().union(*occ.values()) if occ else set()
    broad = {
        gid: sum(1 for gs in occ.values() if gid in gs)
        for gid in all_g
    }
    broad = {g: n for g, n in broad.items() if n >= half}
    lines.append(f"  {metric}: {len(broad)} genes in >= {half} comparisons")
    for gid, n in sorted(broad.items(), key=lambda x: -x[1]):
        ann  = annotation_store.get(gid, {})
        name = ann.get("preferred_name") or gid
        arab = ann.get("blast_arab_id", "")
        desc = (ann.get("tair_short_description") or
                ann.get("eggnog_description", ""))
        desc = (desc[:70] + "…") if len(desc) > 70 else desc
        lines.append(f"    [{n:2d}/{n_pairwise}] {name}  [{arab}]  {desc}")
    lines.append("")

# Exclusive outliers
lines += [sep2,
          "EXCLUSIVE OUTLIERS (in exactly ONE comparison)",
          sep2, ""]
for metric, occ in structures.items():
    excl_per_comp = defaultdict(list)
    for gid in set().union(*occ.values()) if occ else []:
        if sum(1 for gs in occ.values() if gid in gs) == 1:
            comp = next(c for c, gs in occ.items() if gid in gs)
            excl_per_comp[comp].append(gid)
    total_excl = sum(len(v) for v in excl_per_comp.values())
    lines.append(f"  {metric}: {total_excl} exclusive genes")
    for comp in sorted(excl_per_comp):
        lines.append(f"    {comp}: {len(excl_per_comp[comp])} exclusive genes")
    lines.append("")

# Adaptation / reproduction flags
lines += [sep2, "FLAGGED GENES (adaptation / reproduction)", sep2, ""]
for metric, occ in structures.items():
    all_g  = set().union(*occ.values()) if occ else set()
    adapt  = [g for g in all_g
              if annotation_store.get(g, {}).get("adaptation_flag")]
    repro  = [g for g in all_g
              if annotation_store.get(g, {}).get("reproduction_flag")]
    lines.append(f"  {metric}:  adaptation-flagged: {len(adapt)},  "
                 f"reproduction-flagged: {len(repro)}")
    for gid in sorted(adapt):
        ann  = annotation_store.get(gid, {})
        name = ann.get("preferred_name") or gid
        kw   = ann.get("adaptation_keywords", "")
        lines.append(f"    [ADAPT] {name} — {kw}")
    for gid in sorted(repro):
        ann  = annotation_store.get(gid, {})
        name = ann.get("preferred_name") or gid
        kw   = ann.get("reproduction_keywords", "")
        lines.append(f"    [REPRO] {name} — {kw}")
    lines.append("")

# GO / BLAST / TAIR coverage
lines += [sep2, "ANNOTATION COVERAGE", sep2, ""]
for metric, occ in structures.items():
    all_g = set().union(*occ.values()) if occ else set()
    n = len(all_g)
    with_go    = sum(1 for g in all_g if annotation_store.get(g, {}).get("GO_terms"))
    with_blast = sum(1 for g in all_g if annotation_store.get(g, {}).get("blast_arab_id"))
    with_tair  = sum(1 for g in all_g if annotation_store.get(g, {}).get("tair_short_description"))
    lines.append(
        f"  {metric} (n={n}):  "
        f"BLAST hits: {with_blast} ({100*with_blast//max(n,1)}%)  "
        f"TAIR desc: {with_tair} ({100*with_tair//max(n,1)}%)  "
        f"GO terms: {with_go} ({100*with_go//max(n,1)}%)"
    )
lines.append("")

lines += [sep, "END OF REPORT", sep, ""]

with open(os.path.join(output_dir, "summary_stats.txt"),
          "w", encoding="utf-8") as f:
    f.write("\n".join(lines))
print("Wrote: summary_stats.txt")

print(f"\n{'='*70}")
print("009 FST-ONLY SUMMARY COMPLETE")
print(f"All outputs saved to: {output_dir}")
print(f"{'='*70}")
print("\nOutput files:")
for fname in ["master_table.tsv", "gene_comparison_matrix.tsv",
              "gene_annotation_lookup.tsv", "universal_outliers.tsv",
              "exclusive_outliers.tsv", "comparison_counts.tsv",
              "summary_stats.txt"]:
    print(f"  {fname}")
