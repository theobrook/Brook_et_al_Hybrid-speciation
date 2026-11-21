#!/usr/bin/env python3
"""
admixture_q_to_indfile.py

Combines ADMIXTURE replicate .Q files for one K into a single CLUMPP
indfile (DATATYPE 0), in the same format written by STRUCTURE HARVESTER:

    Label  Pop  (Miss)  Ploidy  :  Q1  Q2  ...  QK

with a blank line between replicate blocks.

Usage:
    python3 admixture_q_to_indfile.py K2_rep1.Q K2_rep2.Q ... > K2.indfile
"""

import sys

for qfile in sys.argv[1:]:
    with open(qfile) as f:
        rows = [line.split() for line in f if line.strip()]
    for i, qs in enumerate(rows, 1):
        q_str = " ".join(f"{float(q):.3f}" for q in qs)
        print(f"  {i}   {i} (0) 1 :  {q_str}")
    print()
