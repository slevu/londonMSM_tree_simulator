### Blast
### download from LANL
### All subtype B from 2253 to 3798
### with label: accession - patient ID - country - sample year
### N =5607 seq (10mb)
### including HXB2

rm(list=ls())
library(ape)

B <- read.dna(file = "./data/blast/lanl_b.fasta", format = "fasta")
head(rownames(B))
dim(B)
tail(rownames(B))
## UK = GB
rownames(B)[grep(rownames(B), pattern = 'GB')]
as.character(B[1:6, 1:6])

##- make a fasta file as the database for BLAST
write.dna(as.character(B[1:6, ]), file = "../softwares/ncbi-blast-2.3.0+/db/db.fas", format = "fasta")

# cmd <- "blastn -help"
# cmd <- "echo $PATH"
 cmd1 <- "cd ../softwares/ncbi-blast-2.3.0+/db/; ls" 
# cmd2 <- "makeblastdb -in db.fas -dbtype nucl -parse_seqids"
system(cmd1)
