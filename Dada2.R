library(dada2)
library(ggplot2); packageVersion("ggplot2")


# define filepath
path = "./demux"
list.files(path)

# define filenames
fnFs = sort(list.files(path, pattern="_R1_001.fastq.gz", full.names=TRUE))
sample.names = sapply(strsplit(basename(fnFs),"_"),`[`,1)

# initialize filtered files
filt_path <- file.path(path, "filtered")
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))

# trim to 250 base pairs
out <- filterAndTrim(fwd=fnFs,filt=filtFs,truncLen=250,truncQ=0,trimLef=0,maxLen=Inf,minLen=0,maxN=0,minQ=0,maxEE=Inf)

# learn erros
errF = learnErrors(filtFs, multithread=TRUE)

# dereplicate
derepFs = derepFastq(filtFs, verbose=TRUE)
names(derepFs) = sample.names

# sample inference
dadaFs = dada(derepFs, err=errF, multithread=FALSE)
dadaFs[[1]]

rm(derepFs,errF)

# model error plot
p = plotErrors(dadaFs[[1]], nominalQ=TRUE)
ggsave("dada_errors_F.png", plot=p)

# make tables
seqtab = makeSequenceTable(dadaFs)
dim(seqtab)
saveRDS(seqtab,"seqtab.rds")
write.csv(seqtab,"seqtab.csv")

table(nchar(getSequences(seqtab)))

# remove chimeras
seqtab.nochim = removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
saveRDS(seqtab,"seqtab.nochim.rds")
write.csv(seqtab,"seqtab.nochim.csv")

# summarize 
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoised", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)

rm(dadaFs)

# assign taxonomy
taxtab <- assignTaxonomy(seqtab.nochim, "../Training/silva_nr_v123_train_set.fa.gz", multithread=TRUE)
taxtab <- addSpecies(taxtab, "../Training/silva_species_assignment_v123.fa.gz")
taxtab.print <- taxtab # Removing sequence rownames for display only
rownames(taxtab.print) <- NULL
head(taxtab.print)
saveRDS(taxtab, "taxtab.rds")
write.csv(taxtab,"taxtab.csv")
