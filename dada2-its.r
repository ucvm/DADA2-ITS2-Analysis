################################################################################
#	dada2 ITS pipeline for nemabiome characterization		       #
################################################################################

#	Load libraries
library(dada2);packageVersion("dada2")

library(ShortRead);packageVersion("ShortRead")

library(Biostrings); packageVersion("Biostrings")

wd <- "/home/sgavriliuk/nemabiome/pipelines/dada2"
path <- "/home/sgavriliuk/nemabiome/original_miseq_files"
list.files(path)

fnFs <- sort(list.files(path, pattern = "_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2_001.fastq", full.names = TRUE))

FWD <- "ACGTCTGGTTCAGGGTTGTT"
REV <- "TTAGTTTCTTTTCCTCCGCT"

allOrients <- function(primer) {
	require(Biostrings)
	dna <- DNAString(primer)
	orients <- c(Forward = dna, Complement = complement(dna),
		     Reverse = reverse(dna), RevComp = reverseComplement(dna))
	return(sapply(orients, toString))
}

FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)

fnFs.filtN <- file.path(wd, "filtN", basename(fnFs))
fnRs.filtN <- file.path(wd, "filtN", basename(fnRs))

filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, 
	compress = FALSE, multithread = TRUE)

primerHits <- function(primer, fn) {
	nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
	return(sum(nhits > 0))
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
    FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
    REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
    REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

cutadapt <- "/home/sgavriliuk/miniconda3/envs/qiime2-2018.11/bin/cutadapt"
system2(cutadapt, args = "--version")

path.cut <- file.path(wd, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)

R1.flags <- paste("-g", FWD, "-a", REV.RC) 
R2.flags <- paste("-G", REV, "-A", FWD.RC) 

for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
    FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
    REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
    REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

cutFs <- sort(list.files(path.cut, pattern = "_R1_001.fastq", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_R2_001.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

jpeg('cutFs_quality_profile.jpg')
plotQualityProfile(cutFs[1:2])
dev.off()

jpeg('cutRs_quality_profile.jpg')
plotQualityProfile(cutRs[1:2])
dev.off()

filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 2, maxEE = c(2,2), 
    rm.phix = TRUE, compress = TRUE, multithread = TRUE)  
head(out)

errF <- learnErrors(filtFs, multithread = TRUE); errR <- learnErrors(filtRs, multithread = TRUE)

jpeg("error-f.jpg")
plotErrors(errF, nominalQ = TRUE)
dev.off()

jpeg("error-r.jpg")
plotErrors(errR, nominalQ = TRUE)
dev.off()

derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)

names(derepFs) <- sample.names
names(derepRs) <- sample.names

#	Increasing banded alignment heuristic
setDadaOpt(BAND_SIZE=32)

dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, 
    getN))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged")
rownames(track) <- sample.names
head(track)


################################################################################
#			Assigning Taxonomy				       #
################################################################################
ref_db <- "/home/sgavriliuk/nemabiome/nematode_ITS2_rdp.fasta"
taxa <- assignTaxonomy(seqtab, ref_db, multithread=TRUE, tryRC=TRUE)

taxa.print <- taxa
rownames(taxa.print)
head(taxa.print)

################################################################################
#		Extracting the goods from R				       #
#	giving our seq headers more manageable names (ASV_1, ASV_2...)	       #
################################################################################
asv_seqs <- colnames(seqtab)
asv_headers <- vector(dim(seqtab)[2], mode="character")

for (i in 1:dim(seqtab)[2]) 
{
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

#	making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs_bandsize_32.fa")

#	count table:
asv_tab <- t(seqtab)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASVs_counts_bandsize_32.txt", sep="\t", quote=F)

#	tax table:
asv_tax <- taxa
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax, "ASVs_taxonomy_bandsize_32.txt", sep="\t", quote=F)
