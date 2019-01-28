install.packages("micropan")
library(micropan)

?panPrep #Preparing a FASTA file before starting comparisons of sequences in a pan-genome study.
#we make a folder data -> fasta, proteins,prepped folders and genome_table with basic genome stats, load genome_table
#into R
#look up - case study

for( i in 1:dim(genome_table)[1]){
    cat("Preparing", genome_table$File[i], "...\n")
    in.file <- file.path("data/proteins", genome_table$File[i])
    gid <- genome_table$GID.tag[i]
    out.file <- file.path("data/prepped", genome_table$File[i])
    panPrep(in.file, gid, out.file)
    }
#checkin the data
fdta <- readFasta("data/prepped/IFB5408_GID1.faa")
print(substring(fdta$Header[1:3],1,50))

#performing blast
in.files <- file.path("data/prepped", dir("data/prepped"))
out.folder <- "blast"
blastAllAll(in.files, out.folder)
# readBlastTable('blast/GID1_vs_GID2.txt')

#clutering basing on blasts
blast.files <- file.path("blast", dir("blast"))
blast.distances <- bDist(blast.files)
save(blast.distances, file="res/blast_distances.RData")

#The variable blast.distances is now a data.frame with 3 columns. The
#third column contains the distances, and it is always a good idea to make a
#histogram of these distances to verify that it looks reasonable:
hist(blast.distances[,3], breaks=50, col="tan4",
     xlab="BLAST distance", ylab="Number of distances")
#The shape of this histogram will vary somewhat from study to study, but there should always be a
#large number of very small distances. If not, it means all proteins in all genomes
#are quite different, which is really strange for a pan-genome.

cluster.blast <- bClust(blast.distances, linkage="complete", threshold=0.85)
length(unique(cluster.blast))
print(cluster.blast[1:7])

#creating panMartix
pm.blast <- panMatrix(cluster.blast)
plot(pm.blast)
summary(pm.blast)

#making a Tree, The panTree function will perform an average linkage hierarchical 
#clustering of the genomes based on the computed distances, and return a Pantree object.
blast.tree <- panTree(pm.blast, nboot=100) # tree with bootstrapping
my.lab <- genome_table$Strain
names( my.lab ) <- genome_table$GID.tag
my.col <- genome_table$Color
names( my.col ) <- genome_table$GID.tag
plot(blast.tree, 
     leaf.lab=my.lab, 
     col=my.col,
     xlab="Manhattan distances",
     main="Pectobacterium parmentieri Phylogenetic Tree")

#pan genome - Binomial mixture models
binomix <- binomixEstimate(pm.blast, K.range=2:7)
print(binomix$BIC.table)
summary(binomix)
#An alternative estimate of pan-genome size is obtained by the Chao lower bound estimator:
pan.size <- chao(pm.blast)
print(pan.size)
plot(binomix)

#other
par(mar=c(4,2,2,5), oma=c(2,2,1,1))
r <- rarefaction(pm.blast, n.perm=100)
plot(r)

h <- heaps(pm.blast, n.perm=100)
print(h) # Only alpha is of interest here, if its is above 1.0 - closed pan-genome, 0 - 0.99 - open.

f <- fluidity(pm.blast, n.sim=100)
print(f)

J <- distJaccard(pm.blast)
print(mean(J))
hist(J, breaks=10, col="tan")


