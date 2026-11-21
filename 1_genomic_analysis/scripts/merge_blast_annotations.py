#!/usr/bin/env python3
"""
007_merge_blast_annotations.py
===============================
Merges BLAST (006) with eggNOG (005) annotations for all outlier genes
per metric. Mirrors 010_merge_blast_annotations.py exactly in logic and
output format, with paths updated for the 002_fst_only pipeline.

Input:
  006_blast_arabidopsis/best_arabidopsis_hits.tsv
  005_arabidopsis_orthologs/pairwise_outliers/{pair}/all_outlier_genes_all_metrics.tsv

Output:
  007_merged_annotations/pairwise_outliers/{pair}/{metric}/
    all_outlier_genes_annotated.tsv
    local_adaptation_candidates.tsv
    reproductive_isolation_candidates.tsv
    arabidopsis_ids_for_TAIR.txt
    summary.txt
  007_merged_annotations/pairwise_outliers/{pair}/
    overall_summary.txt
"""

import os, csv, re
from collections import defaultdict

base_dir   = "results/outliers"
blast_file = os.path.join(base_dir, "006_blast_arabidopsis", "best_arabidopsis_hits.tsv")
annot_dir  = os.path.join(base_dir, "005_arabidopsis_orthologs", "pairwise_outliers")
output_dir = os.path.join(base_dir, "007_merged_annotations")
os.makedirs(output_dir, exist_ok=True)

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
]

def flag_keywords(text, kwlist):
    if not text: return []
    return [kw for kw, wb in kwlist if
            (re.search(r'\b' + re.escape(kw) + r'\b', text, re.I) if wb
             else kw.lower() in text.lower())]

OUTPUT_FIELDS = [
    "gene_id", "preferred_name",
    "blast_arab_id", "blast_pident", "blast_qcov", "blast_evalue",
    "blast_bitscore", "blast_arab_description",
    "eggnog_arab_id", "arab_lookup_method",
    "seed_ortholog", "eggnog_description",
    "GO_terms", "COG_category", "KEGG_ko", "KEGG_Pathway",
    "adaptation_flag", "adaptation_keywords",
    "reproduction_flag", "reproduction_keywords",
    "metrics",
]

def build_merged_row(row, hit):
    text = " ".join([row.get("description", ""), row.get("preferred_name", ""),
                     row.get("GO_terms", ""), hit.get("blast_arab_description", "")])
    af = flag_keywords(text, ADAPTATION_KEYWORDS)
    rf = flag_keywords(text, REPRO_KEYWORDS)
    return {
        "gene_id":                row["gene_id"],
        "preferred_name":         row.get("preferred_name", ""),
        "blast_arab_id":          hit.get("blast_arab_id", ""),
        "blast_pident":           hit.get("blast_pident", ""),
        "blast_qcov":             hit.get("blast_qcov", ""),
        "blast_evalue":           hit.get("blast_evalue", ""),
        "blast_bitscore":         hit.get("blast_bitscore", ""),
        "blast_arab_description": hit.get("blast_arab_description", ""),
        "eggnog_arab_id":         row.get("arabidopsis_ortholog", ""),
        "arab_lookup_method":     row.get("arab_lookup_method", ""),
        "seed_ortholog":          row.get("seed_ortholog", ""),
        "eggnog_description":     row.get("description", ""),
        "GO_terms":               row.get("GO_terms", ""),
        "COG_category":           row.get("COG_category", ""),
        "KEGG_ko":                row.get("KEGG_ko", ""),
        "KEGG_Pathway":           row.get("KEGG_Pathway", ""),
        "adaptation_flag":        "YES" if af else "",
        "adaptation_keywords":    "; ".join(af) if af else "",
        "reproduction_flag":      "YES" if rf else "",
        "reproduction_keywords":  "; ".join(rf) if rf else "",
        "metrics":                row.get("metric", row.get("metrics", "")),
    }

# Load BLAST hits
print("Loading BLAST results...")
blast_hits = {}
with open(blast_file, newline="") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        gene = row.get("query", "")
        if gene:
            blast_hits[gene] = {
                "blast_arab_id":          row.get("arabidopsis_id", ""),
                "blast_pident":           row.get("pident", ""),
                "blast_qcov":             row.get("qcov", ""),
                "blast_evalue":           row.get("evalue", ""),
                "blast_bitscore":         row.get("bitscore", ""),
                "blast_arab_description": row.get("description", ""),
            }
print(f"Loaded BLAST hits for {len(blast_hits)} genes")

# Main loop
pair_folders = sorted([d for d in os.listdir(annot_dir)
                       if os.path.isdir(os.path.join(annot_dir, d))])
print(f"Found {len(pair_folders)} pairwise comparisons")

for pair_name in pair_folders:
    annot_file = os.path.join(annot_dir, pair_name,
                              "all_outlier_genes_all_metrics.tsv")
    out_base   = os.path.join(output_dir, "pairwise_outliers", pair_name)
    print(f"\n{'='*60}\nProcessing: {pair_name}\n{'='*60}")

    if not os.path.exists(annot_file):
        print(f"  WARNING: not found: {annot_file}"); continue
    os.makedirs(out_base, exist_ok=True)

    with open(annot_file, newline="") as f:
        all_rows = list(csv.DictReader(f, delimiter="\t"))
    print(f"  Loaded {len(all_rows)} genes")

    all_metrics = sorted(set(
        m.strip()
        for row in all_rows
        for m in row.get("metric", row.get("metrics", "")).split(",")
        if m.strip()
    ))

    rows_by_metric = defaultdict(list)
    for row in all_rows:
        hit    = blast_hits.get(row["gene_id"], {})
        merged = build_merged_row(row, hit)
        for m in [m.strip() for m in
                  row.get("metric", row.get("metrics", "")).split(",")
                  if m.strip()]:
            rows_by_metric[m].append(merged)

    overall = [f"Analysis: {pair_name}", "=" * 60, ""]

    for metric in all_metrics:
        metric_rows = rows_by_metric.get(metric, [])
        if not metric_rows:
            overall.append(f"{metric}: 0 genes"); continue

        mdir = os.path.join(out_base, metric)
        os.makedirs(mdir, exist_ok=True)

        with open(os.path.join(mdir, "all_outlier_genes_annotated.tsv"),
                  "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=OUTPUT_FIELDS, delimiter="\t")
            w.writeheader(); w.writerows(metric_rows)

        for fname, flag in [
            ("local_adaptation_candidates.tsv",      "adaptation_flag"),
            ("reproductive_isolation_candidates.tsv", "reproduction_flag"),
        ]:
            sub = [r for r in metric_rows if r[flag] == "YES"]
            if sub:
                with open(os.path.join(mdir, fname), "w", newline="") as f:
                    w = csv.DictWriter(f, fieldnames=OUTPUT_FIELDS, delimiter="\t")
                    w.writeheader(); w.writerows(sub)

        aids = sorted(set(r["blast_arab_id"] for r in metric_rows
                          if r["blast_arab_id"]))
        if aids:
            with open(os.path.join(mdir, "arabidopsis_ids_for_TAIR.txt"), "w") as f:
                f.write("\n".join(aids) + "\n")

        n_blast = sum(1 for r in metric_rows if r["blast_arab_id"])
        n_ad    = sum(1 for r in metric_rows if r["adaptation_flag"] == "YES")
        n_re    = sum(1 for r in metric_rows if r["reproduction_flag"] == "YES")

        with open(os.path.join(mdir, "summary.txt"), "w") as f:
            f.write(f"Metric: {metric}\nPair: {pair_name}\n{'='*60}\n\n"
                    f"Total outlier genes: {len(metric_rows)}\n"
                    f"BLAST hits: {n_blast}\nArab orthologs: {len(aids)}\n"
                    f"Adaptation candidates: {n_ad}\n"
                    f"Reproduction candidates: {n_re}\n")

        overall.append(f"{metric}: {len(metric_rows)} genes, "
                       f"{n_blast} BLAST, {n_ad} adapt, {n_re} repro")
        print(f"    {metric}: {len(metric_rows)} genes, "
              f"{n_blast} BLAST, {n_ad} adapt, {n_re} repro")

    with open(os.path.join(out_base, "overall_summary.txt"), "w") as f:
        f.write("\n".join(overall) + "\n")

print(f"\n{'='*60}\nMERGING COMPLETED!\nResults: {output_dir}\n{'='*60}")
