#Install BioPython:
sudo apt-get install python-dev python-pip python-setuptools python-reportlab python-matplotlib python-numpy python-scipy python-mysqldb
sudo pip install BioPython -U
sudo apt-get install ncbi-blast+ blast2

#Create separate catalog for analysis

#Place raw data from catalogs pacbio/ (content of Analysis_Results/) and illumina/ (zipped IFB_XXXX.fastq) 

#Merging of FASTQ from PacBio:

cat pacbio/*.fastq > pacbio.fastq

#Prepare pacbio.fofn ('file of file names') that has the directory to all *.bas.h5 in pacbio catalog e.g..
pacbio/m150717_104246_42149_c100818232550000001823171311031517_s1_p0.bas.h5
pacbio/m150725_013821_42149_c100854972550000001823186312311503_s1_p0.bas.h5

#Unzip Illumina:
gunzip illumina/*.gz
bunzip2 illumina/*.bz2

#QC of Illumina:
docker run --rm -u `stat -c "%u:%g" .` -v `pwd`:/home   genome_assembly:0.1   java -Xms1g -Xmx1g -jar /opt/trimmomatic/trimmomatic.jar PE -threads 12 -phred33    illumina/IFB_XXXX_R1.fastq illumina/IFB_XXXX_R2.fastq    illumina_R1.fastq /dev/null illumina_R2.fastq /dev/null ILLUMINACLIP:/opt/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:5:20 MINLEN:100 2> illumina.trimming.txt

#De novo assemlby with SPAdes
docker run --rm -u `stat -c "%u:%g" .` -v `pwd`:/home   genome_assembly:0.1   spades.py -t 12 -m 20 --careful -o spades/ --pe1-1 illumina_R1.fastq --pe1-2 illumina_R2.fastq --pacbio pacbio.fastq

#File spades/scaffolds.fasta has the final genome version. Delete any conting below 500 bp and save as genome_spades.fasta

#Prokka
docker run --rm -u `stat -c "%u:%g" .` -v `pwd`:/home   prokka:1.11   prokka --cpus 12 --compliant --centre XXX --rfam --force --outdir prokka_1/ --prefix IFB_XXXX --genus Pectobacterium --species parmentieri --strain IFB_XXXX genome_spades.fasta

#If you have serveral contigs:
./prepare_prokka_contig.py --prokka-dir prokka_1 --prokka-prefix IFB_XXXX --output-dir prokka --output-prefix IFB_XXXX
docker run --rm -u `stat -c "%u:%g" .` -v `pwd`:/home   prokka:1.11   prokka --cpus 12 --compliant --centre XXX --rfam --force --outdir prokka_2/ --prefix IFB_XXXX --genus Pectobacterium --species parmentieri --strain IFB_XXXX genome.fasta

#In case of full chromosome:
./prepare_prokka.py --prokka-dir prokka_2 --prokka-prefix IFB_XXXX --output-dir prokka --output-prefix IFB_XXXX

#Reorient the genome
./reorient_genome.py -i prokka_1/IFB_XXXX.gbk -g dnaA -o genome.fasta
