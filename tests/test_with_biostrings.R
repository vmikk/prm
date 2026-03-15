## Test with Biostrings package
## Based on https://benjjneb.github.io/dada2/ITS_workflow.html

library(Biostrings)
library(ShortRead)

## FASTQ files
fnFs <- "./testdata/R1.fq.gz"
fnRs <- "./testdata/R2.fq.gz"

## Primer sequences
FWD <- "ASCYGYGGTAAYWCCAGC"
REV <- "TCHNHGNATTTCACCNCT"

## Create all orientations of the input sequence
allOrients <- function(primer) {

    dna <- DNAString(primer)
    orients <- c(
        Forward    = dna,
        Complement = Biostrings::complement(dna),
        Reverse    = Biostrings::reverse(dna),
        RevComp    = Biostrings::reverseComplement(dna))
    
    ## Convert back to character vector
    return(sapply(orients, toString))
}

FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)

FWD.orients
#  Forward     "ASCYGYGGTAAYWCCAGC"
#  Complement  "TSGRCRCCATTRWGGTCG" 
#  Reverse     "CGACCWYAATGGYGYCSA" 
#  RevComp     "GCTGGWRTTACCRCRGST" 

REV.orients
# Forward      "TCHNHGNATTTCACCNCT"
# Complement   "AGDNDCNTAAAGTGGNGA"
# Reverse      "TCNCCACTTTANGHNHCT" 
# RevComp      "AGNGGTGAAATNCDNDGA" 


## Function to counts number of reads in which the primer is found
primerHits <- function(primer, fn) {
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}


## Count the number of times the primers appear in the forward and reverse read,
## while considering all possible primer orientations
rbind(
    FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs),
    FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs),
    REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs),
    REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs))

#                  Forward Complement Reverse RevComp
# FWD.ForwardReads     481          0       0       0
# FWD.ReverseReads     461          0       0       3
# REV.ForwardReads     496          0       0       4
# REV.ReverseReads     479          0       0       0

