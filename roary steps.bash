#roary piplene script
#genomes assembled and annotated with prokka, input files gff
#in case prokka gives you gff3 not gff - use 

#roary_gff.py -> usage: python roary_gff.py --genbank FILE_NAME.gbk --fasta FILE_NAME.fna > FILE_NAME.roary.gff


python roary_gff.py --genbank IFB5408.gb --fasta IFB5408.fasta > IFB5408_roary.gff
python roary_gff.py --genbank IFB5427.gb --fasta IFB5427.fasta > IFB5427_roary.gff
python roary_gff.py --genbank IFB5432.gb --fasta IFB5432.fasta > IFB5432_roary.gff
python roary_gff.py --genbank IFB5441.gb --fasta IFB5441.fasta > IFB5441_roary.gff
python roary_gff.py --genbank IFB5485.gb --fasta IFB5485.fasta > IFB5485_roary.gff
python roary_gff.py --genbank IFB5486.gb --fasta IFB5486.fasta > IFB5486_roary.gff
python roary_gff.py --genbank IFB5597.gb --fasta IFB5597.fasta > IFB5597_roary.gff
python roary_gff.py --genbank IFB5604.gb --fasta IFB5604.fasta > IFB5604_roary.gff
python roary_gff.py --genbank IFB5605.gb --fasta IFB5605.fasta > IFB5605_roary.gff
python roary_gff.py --genbank IFB5619.gb --fasta IFB5619.fasta > IFB5619_roary.gff
python roary_gff.py --genbank IFB5623.gb --fasta IFB5623.fasta > IFB5623_roary.gff
python roary_gff.py --genbank IFB5626.gb --fasta IFB5626.fasta > IFB5626_roary.gff

#roary start, for more options -h will give you more stuff

roary -f ./roary -e -n -r *.gff

#you need to make a tree so you can use it to make roary_plots.py

fasttreeMP -nt -gtr core_gene_alignment.aln > CdiffTree.newick
python roary_plots.py CdiffTee.newick gene_presence_absence.csv

#or to use roary made tree, this is better
python roary_plots.py accessory_binary_genes.fa.newick gene_presence_absence.csv

#he gives all needed files, for additional stuff

create_pan_genome_plots.R


#COMMON ISSUES
#always make folders without SPACE, always UNDERSCORE
#when he is popin error 'Could not extract any protein sequences from /home/HCGS/Desktop/for_Roary/IFB5408.gff. Does the file contain the assembly as well as the annotation?'
#look into gff file, fasta file has to have the same name as the replicon number
