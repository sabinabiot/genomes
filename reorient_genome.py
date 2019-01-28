#!/usr/bin/python
import os
import sys
import argparse
from Bio import SeqIO
from Bio import SeqRecord

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", metavar="FILE", help="Input Genbank file", required=True, action="store", type=str, dest="input")
parser.add_argument("-o", "--output", metavar="FILE", help="Output FASTA file", required=True, action="store", type=str, dest="output")
parser.add_argument("-g", "--gene", metavar="FILE", help="First gene of the reoriented genome", required=True, action="store", type=str, dest="gene")
args = parser.parse_args()

seqrecord = SeqIO.read(args.input, "genbank")
gene = [feature for feature in seqrecord.features if feature.type == "CDS" and feature.qualifiers.get("gene", [""])[0] == args.gene][0]
seq = None
if gene.location.strand == 1:
    start = gene.location.start
    seq = seqrecord.seq[start:] + seqrecord.seq[:start]
if gene.location.strand == -1:
    start = gene.location.end
    seq = seqrecord.seq[start:] + seqrecord.seq[:start]
    seq = seq.reverse_complement()

outfile = open(args.output, "w")
seqrecord = SeqRecord.SeqRecord(id=seqrecord.id, seq=seq, description="")
SeqIO.write(seqrecord, outfile, "fasta")
outfile.close()
