#!/usr/bin/env python
import os
import sys
import urllib
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("--genbank", metavar="FILE", help="Input Genbank file", required=True, action="store", type=str, dest="genbank")
parser.add_argument("--fasta", metavar="FILE", help="Input FASTA file", required=True, action="store", type=str, dest="fasta")
args = parser.parse_args()

print "##gff-version 3"
for genbank_record in SeqIO.parse(args.genbank, "genbank"):
    for feature in genbank_record.features:
        seq_region = genbank_record.id
        start = int(feature.location.start)      
        end = int(feature.location.end)
        strand =  "+"
        if feature.location.strand == -1:
            strand = "-"
        gff3_attributes = "" 
        if feature.type != "CDS": continue
        locus_tag = feature.qualifiers["locus_tag"][0]
        gene_name = feature.qualifiers.get("gene", [""])[0]
        description = feature.qualifiers.get("product", [""])[0]
        protein_id = feature.qualifiers.get("protein_id", [""])[0]
        gff3_attributes = "ID=%s;locus_tag=%s" % (locus_tag, locus_tag)
        if gene_name != "":
            gff3_attributes += ";gene=%s" % (gene_name)
        if description != "":
            gff3_attributes += ";product=%s" % (urllib.quote(description))
        if protein_id != "":
            gff3_attributes += ";protein_id=%s" % (protein_id)
        print "\t".join(( seq_region, "prokka", feature.type, str(start + 1), str(end), ".", strand, ".", gff3_attributes ))

print "##FASTA"
for line in open(args.fasta):
    print line,




