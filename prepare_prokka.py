#!/usr/bin/env python
import os
import sys
import urllib
import argparse
import collections
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


parser = argparse.ArgumentParser()
parser.add_argument("--prokka-dir", metavar="DIR", help="Prokka result directory", required=True, action="store", type=str, dest="prokka_dir")
parser.add_argument("--prokka-prefix", metavar="PREFIX", help="Prokka file prefix", required=True, action="store", type=str, dest="prokka_prefix")
parser.add_argument("--output-dir", metavar="DIR", help="Output file directory", required=True, action="store", type=str, dest="output_dir")
parser.add_argument("--output-prefix", metavar="PREFIX", help="Output file prefix", required=True, action="store", type=str, dest="output_prefix")
args = parser.parse_args()

if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)

genbank_record = SeqIO.read("%s/%s.gbk" % (args.prokka_dir, args.prokka_prefix), "genbank")
genbank_record.id = args.output_prefix
genbank_record.name = args.output_prefix
genbank_record.features = [feature for feature in genbank_record.features if feature.type != "gene"]
for feature in genbank_record.features:
    if feature.qualifiers.has_key("locus_tag"):
        locus_tag = feature.qualifiers["locus_tag"][0].replace("PROKKA", args.output_prefix)
        feature.qualifiers["locus_tag"] = [locus_tag]
    if feature.qualifiers.has_key("protein_id"):
        feature.qualifiers["protein_id"] = feature.qualifiers["locus_tag"]
    if feature.type == "misc_RNA": 
        feature.qualifiers["gene"] = feature.qualifiers["product"]

outfile = open("%s/%s.gbk" % (args.output_dir, args.output_prefix), "w")
SeqIO.write(genbank_record, outfile, "genbank")
outfile.close()

outfile = open("%s/%s.fasta" % (args.output_dir, args.output_prefix), "w")
SeqIO.write(genbank_record, outfile, "fasta")
outfile.close()

outfile = open("%s/%s.ffn" % (args.output_dir, args.output_prefix), "w")
for feature in genbank_record.features:
    if not feature.type in ["CDS", "rRNA", "tRNA", "tmRNA", "misc_RNA"]: continue
    seq = feature.extract(genbank_record.seq)
    seqrecord = SeqRecord(id=feature.qualifiers["locus_tag"][0], description=feature.qualifiers.get("product", [""])[0], seq=seq)
    SeqIO.write(seqrecord, outfile, "fasta")
outfile.close()

outfile = open("%s/%s.faa" % (args.output_dir, args.output_prefix), "w")
for feature in genbank_record.features:
    if feature.type != "CDS": continue
    seq = Seq(feature.qualifiers["translation"][0])
    seqrecord = SeqRecord(id=feature.qualifiers["protein_id"][0], description=feature.qualifiers.get("product", [""])[0], seq=seq)
    SeqIO.write(seqrecord, outfile, "fasta")
outfile.close()

repeat_counter = collections.defaultdict(int)
outfile = open("%s/%s.gff3" % (args.output_dir, args.output_prefix), "w")
for feature in genbank_record.features:
    seq_region = genbank_record.id
    start = int(feature.location.start)      
    end = int(feature.location.end)
    strand =  "+"
    if feature.location.strand == -1:
        strand = "-"
    gff3_attributes = "" 
    if feature.type == "source": continue
    if feature.type in ["CDS", "rRNA", "tRNA", "tmRNA", "misc_RNA"]: 
        locus_tag = feature.qualifiers["locus_tag"][0]
        gene_name = feature.qualifiers.get("gene", [""])[0]
        description = feature.qualifiers.get("product", [""])[0]
        gff3_attributes = "ID=%s;Locus_tag=%s;Name=%s;Description=%s" % (locus_tag, locus_tag, gene_name, urllib.quote(description))
        if feature.type == "CDS":
            protein_id = feature.qualifiers.get("protein_id", [""])[0]
            gff3_attributes += ";Protein_ID=%s" % (protein_id,)
    if feature.type == "repeat_region":
        repeat_counter[args.output_prefix] += 1
        repeat_id = "%s_CRISPR_%s" % (args.output_prefix, str(repeat_counter[args.output_prefix]).zfill(3))
        repeat_family = feature.qualifiers.get("rpt_family", [""])[0]
        gff3_attributes = "ID=%s;Repeat_family=%s" % (repeat_id, repeat_family)
    outfile.write( "\t".join(( seq_region, "Prokka", feature.type, str(start + 1), str(end), ".", strand, ".", gff3_attributes )) + "\n" )
outfile.close()




