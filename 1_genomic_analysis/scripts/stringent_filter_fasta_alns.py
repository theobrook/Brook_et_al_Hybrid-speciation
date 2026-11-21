#!/usr/bin/env python3

"""
Filter FASTA alignments for phylogenomic analysis

This script filters aligned FASTA files based on missing data and heterozygosity.
It performs three stages of filtering:
1. Remove columns with excessive missing data
2. Remove sequences with excessive missing data
3. Remove sequences with excessive heterozygosity

Author: Owen G Osborne (modified)
"""

from Bio import AlignIO
import glob
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import numpy as np
import argparse
import sys
import os

parser = argparse.ArgumentParser(
    description='Filter aligned FASTA files for missing data and heterozygous sites.',
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog="""
Examples:
  # Use default parameters (moderate filtering)
  python3 filter_fasta_alns_updated.py -o filtered_fastas/

  # Conservative filtering
  python3 filter_fasta_alns_updated.py -o filtered_fastas/ -g 30 -G 50 -H 10 -l 500 -s 35

  # Relaxed filtering
  python3 filter_fasta_alns_updated.py -o filtered_fastas/ -g 70 -G 80 -H 20 -l 200 -s 25
"""
)

parser.add_argument(
    "-i", "--input_extension",
    help="Input file extension (default: fas)",
    default="fas"
)
parser.add_argument(
    "-o", "--outdir",
    required=True,
    help="Path to output folder (will be created if it doesn't exist)"
)
parser.add_argument(
    "-l", "--min_length",
    help="Minimum alignment length in bp (default: 300)",
    default=300,
    type=int
)
parser.add_argument(
    "-g", "--max_gaps_per_site",
    help="Maximum percentage of missing data/gaps per column (default: 50)",
    default=50,
    type=int
)
parser.add_argument(
    "-G", "--max_gaps_per_sequence",
    help="Maximum percentage of missing data/gaps per sequence (default: 75)",
    default=75,
    type=int
)
parser.add_argument(
    "-H", "--max_het_per_sequence",
    help="Maximum percentage of heterozygous sites per sequence (default: 15)",
    default=15,
    type=int
)
parser.add_argument(
    "-s", "--min_seqs",
    help="Minimum number of sequences to retain alignment (default: 30)",
    default=30,
    type=int
)

args = parser.parse_args()

# Setup
files = f'*.{args.input_extension}'
outdir = args.outdir
if not outdir.endswith(os.path.sep):
    outdir += os.path.sep

# Create output directory if it doesn't exist
os.makedirs(outdir, exist_ok=True)

minLen = args.min_length
maxpcGap = args.max_gaps_per_sequence
minSeqs = args.min_seqs
maxpcGap_s = args.max_gaps_per_site
maxpcHet = args.max_het_per_sequence
numfiles = len(glob.glob(files))
logpath = f'{outdir}fasta_alignment_filtering.log'

# Print settings
print(f"\n{'='*70}")
print(f"ALIGNMENT FILTERING SETTINGS")
print(f"{'='*70}")
print(f"Input files found: {numfiles} (*.{args.input_extension})")
print(f"Output directory: {outdir}")
print(f"\nFiltering parameters:")
print(f"  Max missing data per column: {maxpcGap_s}%")
print(f"  Max missing data per sequence: {maxpcGap}%")
print(f"  Max heterozygosity per sequence: {maxpcHet}%")
print(f"  Minimum alignment length: {minLen} bp")
print(f"  Minimum sequences per alignment: {minSeqs}")
print(f"\nLog file: {logpath}")
print(f"{'='*70}\n")

# Initialize log file
log = open(logpath, 'w')
log.write("filename,prefilter_len,postfilter_len,passed_len_filter,prefilter_nSeqs,"
          "nSeqs_removed_gap_filt,nSeqs_removed_het_filt,postfilter_nSeqs,passed_nSeq_filter\n")

if numfiles > 0:
    nPass = 0
    nFail_length = 0
    nFail_nseqs = 0

    for fasta in sorted(glob.glob(files)):
        fastaname = fasta.split('/')[-1]

        try:
            aln = AlignIO.read(fasta, "fasta")
        except Exception as e:
            print(f"Warning: Skipping {fastaname} - {e}")
            continue

        seqLen1 = aln.get_alignment_length()
        nSeq1 = len(aln)

        # Stage 1: Column-wise filtering
        aln_a = np.array([list(rec) for rec in aln], dtype=np.str_)
        delsites = []

        for i in range(seqLen1):
            site = list(aln_a[:, i])
            nSeq = len(site)
            nGap_s = float(site.count('N') + site.count('n') + site.count('-'))
            pcGap_s = nGap_s / nSeq * 100
            if pcGap_s > maxpcGap_s:
                delsites.append(i)

        aln_a = np.delete(aln_a, delsites, 1)

        # Reconstruct alignment after column filtering
        filt1_aln = MultipleSeqAlignment([
            SeqRecord(Seq(''.join(aln_a[c,])), id=aln[c].id)
            for c in range(len(aln))
        ])

        length = filt1_aln.get_alignment_length()

        # Check minimum length
        if length >= minLen:
            # Stage 2: Sequence-wise gap filtering
            filt2_aln = MultipleSeqAlignment([])
            for current_seq in filt1_aln:
                seq = str(current_seq.seq)
                nGap = float(seq.count('N') + seq.count('n') + seq.count('-'))
                pcGap = nGap / length * 100
                if pcGap < maxpcGap:
                    filt2_aln.append(SeqRecord(Seq(seq), id=current_seq.id))

            # Stage 3: Heterozygosity filtering
            filt3_aln = MultipleSeqAlignment([])
            for current_seq in filt2_aln:
                seq = str(current_seq.seq)
                nHet = float(sum(seq.count(base) for base in 'WwSsMmKkRrYy'))
                nATGC = float(sum(seq.count(base) for base in 'ATGCatgc'))

                if nATGC + nHet > 0:
                    pcHet = nHet / (nHet + nATGC) * 100
                else:
                    pcHet = 0

                if pcHet < maxpcHet:
                    filt3_aln.append(SeqRecord(Seq(seq), id=current_seq.id))

            postFiltSeqs = len(filt3_aln)

            # Check minimum sequences
            if postFiltSeqs >= minSeqs:
                outpath = f'{outdir}{fastaname}'
                with open(outpath, 'w') as f:
                    AlignIO.write(filt3_aln, f, "fasta")
                pasN = "Y"
                nPass += 1
            else:
                pasN = "N"
                nFail_nseqs += 1

            pasLen = "Y"
            remGap = len(aln) - len(filt2_aln)
            remHet = len(filt2_aln) - len(filt3_aln)
        else:
            pasLen = "N"
            nFail_length += 1
            remGap = "NA"
            remHet = "NA"
            postFiltSeqs = "NA"
            pasN = "NA"

        # Write to log
        log.write(f"{fastaname},{seqLen1},{length},{pasLen},{nSeq1},"
                 f"{remGap},{remHet},{postFiltSeqs},{pasN}\n")

    log.close()

    # Summary
    print(f"\n{'='*70}")
    print(f"FILTERING SUMMARY")
    print(f"{'='*70}")
    print(f"Total alignments processed: {numfiles}")
    print(f"Alignments passed all filters: {nPass} ({nPass/numfiles*100:.1f}%)")
    print(f"Failed length filter: {nFail_length}")
    print(f"Failed minimum sequences filter: {nFail_nseqs}")
    print(f"\nDetailed log: {logpath}")
    print(f"{'='*70}\n")

else:
    print(f"\nNo fasta files found with extension '.{args.input_extension}'")
    print("Check that the working directory contains fasta files and --input_extension is set correctly\n")
    sys.exit(1)
