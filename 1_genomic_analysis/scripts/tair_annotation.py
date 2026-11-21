#!/usr/bin/env python3
"""
008_tair_annotation.py
=======================
Joins TAIR data with merged annotations from 007_merge_blast_annotations.py.
Mirrors 011_tair_annotation.py exactly, with paths updated for 002_fst_only.

Requires TAIR data files in 008_tair_annotation/tair_data/:
  Araport11_functional_descriptions_20241231.txt
  gene_aliases_20241231.txt
  Locus_Germplasm_Phenotype_20241231.txt

These are the same files used by 011_tair_annotation.py — just symlink or
copy from the main pipeline's tair_data folder rather than re-downloading.

Input:
  007_merged_annotations/pairwise_outliers/{pair}/{metric}/
    all_outlier_genes_annotated.tsv

Output:
  008_tair_annotation/pairwise_outliers/{pair}/{metric}/
    all_outlier_genes_tair_annotated.tsv
    summary.txt
  008_tair_annotation/pairwise_outliers/{pair}/
    overall_summary.txt
"""

import os, csv, re
from collections import defaultdict

base_dir   = "results/outliers"
tair_dir   = os.path.join(base_dir, "008_tair_annotation", "tair_data")
merged_dir = os.path.join(base_dir, "007_merged_annotations", "pairwise_outliers")
output_dir = os.path.join(base_dir, "008_tair_annotation")
os.makedirs(output_dir, exist_ok=True)

# ── TAIR data ──────────────────────────────────────────────────────────────
# NOTE: The TAIR files are large. Rather than copying them, create symlinks:
#   cd 008_tair_annotation && mkdir -p tair_data && cd tair_data
#   ln -s ../../.../001_pairwise_outliers/011_tair_annotation/tair_data/* .

func_desc_file = os.path.join(tair_dir, "Araport11_functional_descriptions_20241231.txt")
aliases_file   = os.path.join(tair_dir, "gene_aliases_20241231.txt")
phenotype_file = os.path.join(tair_dir, "Locus_Germplasm_Phenotype_20241231.txt")

tair_descriptions = {}
tair_aliases      = {}
tair_phenotypes   = defaultdict(set)

print("Loading TAIR functional descriptions...")
if os.path.exists(func_desc_file):
    with open(func_desc_file, encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith(("Locus", "locus", "#")): continue
            fields = line.strip().split("\t")
            if len(fields) < 4: continue
            lb = re.sub(r'\.\d+$', '', fields[0].strip())
            sd = fields[3].strip() if len(fields) > 3 else ""
            cs = fields[4].strip() if len(fields) > 4 else ""
            cd = fields[5].strip() if len(fields) > 5 else ""
            if lb not in tair_descriptions:
                tair_descriptions[lb] = {
                    "tair_short_description":       sd,
                    "tair_curator_summary":         cs,
                    "tair_computational_description": cd,
                }
            else:
                e = tair_descriptions[lb]
                if cs and not e["tair_curator_summary"]:   e["tair_curator_summary"] = cs
                if sd and not e["tair_short_description"]: e["tair_short_description"] = sd
    print(f"  Loaded {len(tair_descriptions)} loci")
else:
    print(f"  WARNING: not found: {func_desc_file}")
    print("  Tip: symlink from main pipeline's 011_tair_annotation/tair_data/")

print("Loading gene aliases...")
if os.path.exists(aliases_file):
    with open(aliases_file, encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith(("locus", "#")): continue
            fields = line.strip().split("\t")
            if len(fields) < 2: continue
            lb  = re.sub(r'\.\d+$', '', fields[0].strip())
            sym = fields[1].strip()
            fn  = fields[2].strip() if len(fields) > 2 else ""
            als = f"{sym} ({fn})" if fn and fn != sym else sym
            tair_aliases.setdefault(lb, [])
            if als not in tair_aliases[lb]: tair_aliases[lb].append(als)
    print(f"  Loaded {len(tair_aliases)} loci")
else:
    print(f"  WARNING: not found: {aliases_file}")

print("Loading phenotype data...")
if os.path.exists(phenotype_file):
    with open(phenotype_file, encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith(("Locus", "locus", "#")): continue
            fields = line.strip().split("\t")
            if len(fields) < 5: continue
            lb  = re.sub(r'\.\d+$', '', fields[0].strip())
            phe = fields[4].strip() if len(fields) > 4 else ""
            if phe and phe != "no phenotype observed":
                if len(phe) > 200: phe = phe[:200] + "..."
                tair_phenotypes[lb].add(phe)
    print(f"  Loaded {len(tair_phenotypes)} loci")
else:
    print(f"  WARNING: not found: {phenotype_file}")

OUTPUT_FIELDS = [
    "gene_id", "preferred_name",
    "blast_arab_id", "blast_pident", "blast_qcov", "blast_evalue",
    "blast_arab_description",
    "tair_symbol_aliases", "tair_short_description",
    "tair_curator_summary", "tair_phenotypes",
    "eggnog_description", "seed_ortholog",
    "GO_terms", "KEGG_ko", "KEGG_Pathway", "COG_category",
    "adaptation_flag", "adaptation_keywords",
    "reproduction_flag", "reproduction_keywords",
    "metrics",
]

def get_tair(arab_id):
    if not arab_id: return {}, [], set()
    lb = re.sub(r'\.\d+$', '', arab_id.upper())
    return (tair_descriptions.get(lb, {}),
            tair_aliases.get(lb, []),
            tair_phenotypes.get(lb, set()))

# Main loop
pair_folders = sorted([d for d in os.listdir(merged_dir)
                       if os.path.isdir(os.path.join(merged_dir, d))])
print(f"\nFound {len(pair_folders)} pairwise comparisons")

for pair_name in pair_folders:
    merged_base = os.path.join(merged_dir, pair_name)
    out_base    = os.path.join(output_dir, "pairwise_outliers", pair_name)
    print(f"\n{'='*60}\nProcessing: {pair_name}\n{'='*60}")
    os.makedirs(out_base, exist_ok=True)

    metric_folders = sorted([
        d for d in os.listdir(merged_base)
        if os.path.isdir(os.path.join(merged_base, d))
    ])

    overall = [f"Analysis: {pair_name}", "=" * 60, ""]

    for metric in metric_folders:
        input_file = os.path.join(merged_base, metric,
                                  "all_outlier_genes_annotated.tsv")
        if not os.path.exists(input_file):
            overall.append(f"{metric}: no input file"); continue

        with open(input_file, newline="") as f:
            rows = list(csv.DictReader(f, delimiter="\t"))
        if not rows:
            overall.append(f"{metric}: 0 genes"); continue

        mdir = os.path.join(out_base, metric)
        os.makedirs(mdir, exist_ok=True)

        out_rows = []
        n_tair   = 0

        for row in rows:
            arab_id        = row.get("blast_arab_id", "")
            desc, als, phe = get_tair(arab_id)
            trow = {
                "gene_id":                row.get("gene_id", ""),
                "preferred_name":         row.get("preferred_name", ""),
                "blast_arab_id":          arab_id,
                "blast_pident":           row.get("blast_pident", ""),
                "blast_qcov":             row.get("blast_qcov", ""),
                "blast_evalue":           row.get("blast_evalue", ""),
                "blast_arab_description": row.get("blast_arab_description", ""),
                "tair_symbol_aliases":    "; ".join(als) if als else "",
                "tair_short_description": desc.get("tair_short_description", ""),
                "tair_curator_summary":   desc.get("tair_curator_summary", ""),
                "tair_phenotypes":        " | ".join(sorted(phe)) if phe else "",
                "eggnog_description":     row.get("eggnog_description", ""),
                "seed_ortholog":          row.get("seed_ortholog", ""),
                "GO_terms":               row.get("GO_terms", ""),
                "KEGG_ko":                row.get("KEGG_ko", ""),
                "KEGG_Pathway":           row.get("KEGG_Pathway", ""),
                "COG_category":           row.get("COG_category", ""),
                "adaptation_flag":        row.get("adaptation_flag", ""),
                "adaptation_keywords":    row.get("adaptation_keywords", ""),
                "reproduction_flag":      row.get("reproduction_flag", ""),
                "reproduction_keywords":  row.get("reproduction_keywords", ""),
                "metrics":                row.get("metrics", ""),
            }
            if desc.get("tair_short_description") or als or phe: n_tair += 1
            out_rows.append(trow)

        with open(os.path.join(mdir, "all_outlier_genes_tair_annotated.tsv"),
                  "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=OUTPUT_FIELDS, delimiter="\t")
            w.writeheader(); w.writerows(out_rows)

        n_blast = sum(1 for r in out_rows if r["blast_arab_id"])

        with open(os.path.join(mdir, "summary.txt"), "w") as f:
            f.write(f"Metric: {metric}\nPair: {pair_name}\n{'='*60}\n\n"
                    f"Total outlier genes: {len(out_rows)}\n"
                    f"BLAST hits: {n_blast}\nTAIR annotated: {n_tair}\n\n")
            curated = [r for r in out_rows if r["tair_curator_summary"]]
            if curated:
                f.write(f"Genes with TAIR curator summaries ({len(curated)}):\n")
                for r in curated:
                    nm = r["preferred_name"] or r["gene_id"]
                    f.write(f"  {nm} ({r['blast_arab_id']}): "
                            f"{r['tair_short_description']}\n"
                            f"    {r['tair_curator_summary'][:300]}\n")
                    if r["tair_phenotypes"]:
                        f.write(f"    Phenotypes: {r['tair_phenotypes'][:300]}\n")
                    f.write("\n")

        overall.append(f"{metric}: {len(out_rows)} genes, "
                       f"{n_blast} BLAST, {n_tair} TAIR")
        print(f"    {metric}: {len(out_rows)} genes, "
              f"{n_blast} BLAST, {n_tair} TAIR")

    with open(os.path.join(out_base, "overall_summary.txt"), "w") as f:
        f.write("\n".join(overall) + "\n")

print(f"\n{'='*60}\nTAIR ANNOTATION COMPLETED!\nResults: {output_dir}\n{'='*60}")
print("\nNOTE: If TAIR data files were missing, symlink them from the main pipeline:")
print("  cd 008_tair_annotation && mkdir -p tair_data && cd tair_data")
print("  ln -s ../../../../001_pairwise_outliers/011_tair_annotation/tair_data/* .")
