###--- generate sequences
### phangorn simSeq function
### need a substitution rate:
  ### either a consensus estimate from the litterature:
    ### In Lythgoe & Fraser, between host synonymous = 1.3e-3,
    ### non-synonymous = 3.4e-3 
  ### or an output from LSD: 10e-5 as default it seems ?
### from Erik:
### Here is a code snipped where I simulated 
### a sequence alignment with 1000 bases and 
### recompute the distance matrix with ape: 
  ### phydat <- simSeq(rvtree_outtree, rate = MU,  
  ### l = 1e3, type="DNA", bf=rep(.25, 4 ), Q=NULL)
  ### msa <- as.DNAbin(phydat)
  ### D <- dist.dna( msa, model = 'JC69' )
### rvtree_outtree is a ape::phylo, 
### and you can swap in the simulated time trees. 

library(phangorn)
?simSeq