#!/usr/bin/env python3

from Bio import AlignIO
import glob
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
import argparse
import sys

parser = argparse.ArgumentParser(description='Concatenate multiple fasta alignments with missing sequences into a single alignment. If taxa are missing, they are replaced by Ns.')
parser.add_argument("-i", "--input_extension", help="Input file extension (default: 'fas'). Files must be in current working directory.", default="fas")
parser.add_argument("-o", "--outfile", required=True, help="Path to output file.")
parser.add_argument("-t", "--taxa", required=True, help="Comma-delimited taxon names.", type=str)

args = parser.parse_args()

files = f'*.{args.input_extension}'
outfile = args.outfile
taxa = args.taxa.split(",")
numfiles = len(glob.glob(files))
numtaxa = len(taxa)

print(f'{numfiles} fasta files and {numtaxa} taxa found, alignments will be concatenated and written to {outfile}')

if numfiles > 0:
    # Create an initial alignment with empty sequences for each taxon
    cataln = MultipleSeqAlignment([SeqRecord(Seq(""), id=taxon) for taxon in taxa])

    for fasta in glob.glob(files):
        try:
            aln = AlignIO.read(fasta, "fasta")
        except ValueError as e:
            print(f"Warning: Skipping file {fasta} due to an error: {e}")
            continue  # Skip this file and move to the next one

        seqLen = aln.get_alignment_length()
        newaln = MultipleSeqAlignment([])

        for catrec in cataln:
            catid = str(catrec.id)
            seq = "N" * seqLen  # Default to N if no match is found
            for rec in aln:
                if str(rec.id) == catid:
                    seq = str(rec.seq)
                    break

            # Concatenate sequences
            newseq = str(catrec.seq) + seq
            newaln.append(SeqRecord(Seq(newseq), id=catid))

        # Update the concatenated alignment with the new sequences
        cataln = newaln

    # Write the concatenated alignment to the output file
    with open(outfile, 'w') as f:
        AlignIO.write(cataln, f, "fasta")

    print("Finished")
    sys.exit()

else:
    print("\nNo fasta files found. Check that the working directory contains fasta files and --input_extension is set correctly.")
    sys.exit()
