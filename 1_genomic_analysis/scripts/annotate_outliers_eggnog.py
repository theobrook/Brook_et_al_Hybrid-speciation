#!/usr/bin/env python3
"""
005_arabidopsis_orthologs.py
============================
Annotates FST-only outlier transcripts with eggNOG data and flags
adaptation/reproduction candidates.

Mirrors 008_arabidopsis_orthologs.py exactly in logic and output format,
but reads from the 002_fst_only/003_annotation/ tree.

Processes both FST_top10 and FST_top5 metric folders for all 15 pairs.

Output: 002_fst_only/005_arabidopsis_orthologs/pairwise_outliers/{pair}/{metric}/
          all_outlier_genes_annotated.tsv
          local_adaptation_candidates.tsv
          reproductive_isolation_candidates.tsv
          arabidopsis_ids_for_TAIR.txt
        002_fst_only/005_arabidopsis_orthologs/pairwise_outliers/{pair}/
          all_outlier_genes_all_metrics.tsv
"""

import os
import csv
import re
from collections import defaultdict

# ── Paths ──────────────────────────────────────────────────────────────────

base_dir       = "results/outliers"
eggnog_file    = "results/eggnog/coprosma_annotation.emapper.annotations"
annotation_dir = os.path.join(base_dir, "003_annotation", "pairwise_outliers")
output_dir     = os.path.join(base_dir, "005_arabidopsis_orthologs")
os.makedirs(output_dir, exist_ok=True)

# ── Keyword lists (identical to 008_arabidopsis_orthologs.py) ──────────────

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
    ("hydroxyphenylpyruvate dioxygenase", False),
    ("NRT1", True), ("nitrate transport", False),
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

GENE_NAME_TO_ARAB = {
    "ADR1":"AT1G33560","HSP90":"AT5G56030","HSP70":"AT3G12580",
    "HSP101":"AT1G74310","HSP81":"AT5G56030",
    "CBF1":"AT4G25490","CBF2":"AT4G25470","CBF3":"AT4G25480",
    "DREB1A":"AT4G25480","DREB2A":"AT5G05410",
    "SOS1":"AT2G01980","NHX1":"AT5G27150","HKT1":"AT4G10310",
    "ICE1":"AT3G26744","COR15A":"AT2G42540",
    "RD29A":"AT5G52310","RD29B":"AT5G52300",
    "FLC":"AT5G10140","FT":"AT1G65480","CO":"AT5G15840",
    "SOC1":"AT2G45660","GI":"AT1G22770","SVP":"AT2G22540",
    "FRI":"AT4G00650","FCA":"AT4G16280","FPA":"AT2G43410",
    "AP1":"AT1G69120","AP2":"AT4G36920","AP3":"AT3G54340",
    "AG":"AT4G18960","PI":"AT5G20240","SEP1":"AT5G15800","SEP3":"AT1G24260",
    "CCA1":"AT2G46830","LHY":"AT1G01060","TOC1":"AT5G61380",
    "PHYA":"AT1G09570","PHYB":"AT2G18790",
    "CRY1":"AT4G08920","CRY2":"AT1G04400",
}

def flag_keywords(text, kwlist):
    if not text: return []
    return [kw for kw, wb in kwlist if
            (re.search(r'\b' + re.escape(kw) + r'\b', text, re.I) if wb
             else kw.lower() in text.lower())]

def extract_arab_id(seed):
    if not seed or seed == "-": return ""
    if "3702." in seed:
        parts = seed.split(".")
        if len(parts) >= 2: return ".".join(parts[1:])
    m = re.search(r'(AT[1-5CM]G\d{5})', seed, re.I)
    return m.group(1).upper() if m else ""

def lookup_arab(seed, pref):
    aid = extract_arab_id(seed)
    if aid: return aid, "seed_ortholog"
    if pref and pref != "-":
        u = pref.upper().strip()
        if u in GENE_NAME_TO_ARAB: return GENE_NAME_TO_ARAB[u], "gene_name_lookup"
        b = re.sub(r'[\.\-]?\d+$', '', u)
        if b in GENE_NAME_TO_ARAB: return GENE_NAME_TO_ARAB[b], "gene_name_lookup"
    return "", ""

FIELDNAMES = [
    "gene_id", "preferred_name", "arabidopsis_ortholog", "arab_lookup_method",
    "seed_ortholog", "description", "category_flag",
    "adaptation_keywords", "reproduction_keywords",
    "metric", "GO_terms", "COG_category", "KEGG_ko", "KEGG_Pathway",
]

# ── Load eggNOG ────────────────────────────────────────────────────────────

print("Loading eggNOG annotations...")
eggnog = {}
with open(eggnog_file) as f:
    for line in f:
        if line.startswith("#"): continue
        fields = line.strip().split("\t")
        if len(fields) < 13: continue
        eggnog[fields[0]] = {
            "seed_ortholog":  fields[1]  if len(fields) > 1  else "",
            "COG_category":   fields[6]  if len(fields) > 6  else "",
            "Description":    fields[7]  if len(fields) > 7  else "",
            "Preferred_name": fields[8]  if len(fields) > 8  else "",
            "GOs":            fields[9]  if len(fields) > 9  else "",
            "KEGG_ko":        fields[11] if len(fields) > 11 else "",
            "KEGG_Pathway":   fields[12] if len(fields) > 12 else "",
        }
print(f"Loaded {len(eggnog)} eggNOG entries")

# ── Annotation function ────────────────────────────────────────────────────

def annotate(gene_list, metric):
    rows = []
    for gene in gene_list:
        a    = eggnog.get(gene, {})
        seed = a.get("seed_ortholog", "")
        desc = a.get("Description", "")
        pref = a.get("Preferred_name", "")
        gos  = a.get("GOs", "")
        ko   = a.get("KEGG_ko", "")
        kp   = a.get("KEGG_Pathway", "")
        cog  = a.get("COG_category", "")
        aid, am = lookup_arab(seed, pref)
        text = " ".join([desc, pref, gos])
        af   = flag_keywords(text, ADAPTATION_KEYWORDS)
        rf   = flag_keywords(text, REPRO_KEYWORDS)
        cat  = ("ADAPTATION + REPRODUCTION" if af and rf else
                "LOCAL_ADAPTATION"          if af else
                "REPRODUCTIVE_ISOLATION"    if rf else "")
        rows.append({
            "gene_id":                gene,
            "preferred_name":         pref if pref != "-" else "",
            "arabidopsis_ortholog":   aid,
            "arab_lookup_method":     am,
            "seed_ortholog":          seed if seed != "-" else "",
            "description":            desc if desc != "-" else "",
            "category_flag":          cat,
            "adaptation_keywords":    "; ".join(af) if af else "",
            "reproduction_keywords":  "; ".join(rf) if rf else "",
            "metric":                 metric,
            "GO_terms":               gos  if gos  != "-" else "",
            "COG_category":           cog  if cog  != "-" else "",
            "KEGG_ko":                ko   if ko   != "-" else "",
            "KEGG_Pathway":           kp   if kp   != "-" else "",
        })
    return rows

def write_outputs(rows, out_dir, metric, label):
    os.makedirs(out_dir, exist_ok=True)
    with open(os.path.join(out_dir, "all_outlier_genes_annotated.tsv"),
              "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDNAMES, delimiter="\t")
        w.writeheader(); w.writerows(rows)

    for fname, flag in [
        ("local_adaptation_candidates.tsv",       "ADAPTATION"),
        ("reproductive_isolation_candidates.tsv",  "REPRODUCT"),
    ]:
        sub = [r for r in rows if flag in r["category_flag"]]
        if sub:
            with open(os.path.join(out_dir, fname), "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=FIELDNAMES, delimiter="\t")
                w.writeheader(); w.writerows(sub)

    aids = sorted(set(r["arabidopsis_ortholog"] for r in rows
                      if r["arabidopsis_ortholog"]))
    if aids:
        with open(os.path.join(out_dir, "arabidopsis_ids_for_TAIR.txt"), "w") as f:
            f.write("\n".join(aids) + "\n")

    n_ad = sum(1 for r in rows if "ADAPTATION" in r["category_flag"])
    n_re = sum(1 for r in rows if "REPRODUCT"  in r["category_flag"])
    print(f"    {metric}: {len(rows)} genes, {n_ad} adapt, {n_re} repro")

# ── Main loop ──────────────────────────────────────────────────────────────

pair_folders = sorted([d for d in os.listdir(annotation_dir)
                       if os.path.isdir(os.path.join(annotation_dir, d))])
print(f"Found {len(pair_folders)} pairwise comparison folders")

for pair_name in pair_folders:
    print(f"\n{'='*60}\nProcessing: {pair_name}\n{'='*60}")
    annot_base = os.path.join(annotation_dir, pair_name)
    out_base   = os.path.join(output_dir, "pairwise_outliers", pair_name)
    os.makedirs(out_base, exist_ok=True)

    all_genes_metrics = defaultdict(set)

    for metric in sorted(os.listdir(annot_base)):
        metric_dir   = os.path.join(annot_base, metric)
        outlier_file = os.path.join(metric_dir, "outlier_transcripts.txt")
        if not os.path.isdir(metric_dir) or not os.path.exists(outlier_file):
            continue
        with open(outlier_file) as f:
            genes = [l.strip() for l in f if l.strip()]
        if not genes:
            print(f"    {metric}: 0 genes, skipping"); continue
        for g in genes:
            all_genes_metrics[g].add(metric)
        rows = annotate(genes, metric)
        write_outputs(rows, os.path.join(out_base, metric), metric, pair_name)

    # Combined file across metrics
    if all_genes_metrics:
        combined_path = os.path.join(out_base, "all_outlier_genes_all_metrics.tsv")
        with open(combined_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=FIELDNAMES, delimiter="\t")
            w.writeheader()
            for gene in sorted(all_genes_metrics):
                a    = eggnog.get(gene, {})
                seed = a.get("seed_ortholog", "")
                desc = a.get("Description", "")
                pref = a.get("Preferred_name", "")
                gos  = a.get("GOs", "")
                ko   = a.get("KEGG_ko", "")
                kp   = a.get("KEGG_Pathway", "")
                cog  = a.get("COG_category", "")
                aid, am = lookup_arab(seed, pref)
                text = " ".join([desc, pref, gos])
                af   = flag_keywords(text, ADAPTATION_KEYWORDS)
                rf   = flag_keywords(text, REPRO_KEYWORDS)
                cat  = ("ADAPTATION + REPRODUCTION" if af and rf else
                        "LOCAL_ADAPTATION"          if af else
                        "REPRODUCTIVE_ISOLATION"    if rf else "")
                w.writerow({
                    "gene_id":               gene,
                    "preferred_name":        pref if pref != "-" else "",
                    "arabidopsis_ortholog":  aid,
                    "arab_lookup_method":    am,
                    "seed_ortholog":         seed if seed != "-" else "",
                    "description":           desc if desc != "-" else "",
                    "category_flag":         cat,
                    "adaptation_keywords":   "; ".join(af) if af else "",
                    "reproduction_keywords": "; ".join(rf) if rf else "",
                    "metric": ", ".join(sorted(all_genes_metrics[gene])),
                    "GO_terms":     gos if gos != "-" else "",
                    "COG_category": cog if cog != "-" else "",
                    "KEGG_ko":      ko  if ko  != "-" else "",
                    "KEGG_Pathway": kp  if kp  != "-" else "",
                })
        print(f"  Combined: {len(all_genes_metrics)} unique genes")

print(f"\n{'='*60}\nALL ANNOTATIONS COMPLETED!\nResults: {output_dir}\n{'='*60}")
