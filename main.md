# Microbiome Stress Project processing pipeline: <br>
### &nbsp; &nbsp; &nbsp; &nbsp; from raw interleaved fastq’s from SRA to merged biomes, 
### &nbsp; &nbsp; &nbsp; &nbsp; representative dada2 ESVs, and corresponding phylogeny.

---

## Download fastq’s from NCBI SRA

1. Enter the accession number provided on the SRA website: https://www.ncbi.nlm.nih.gov/Traces/study/?go=home

2. Select all the runs (click green plus sign) and download the "RunInfo Table". This can provide some metadata information (treatments…).

3. Copy all the "Experiment" accession numbers (starting with SRX) from the "RunInfo Table" - usually there is 1 experiment number per sample.

4. Paste these "experiment" numbers separated by a comma at the link below. Ensure no extra-space between the comma and the number or SRA will not recognize it.
https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=search_seq_name

5. Download the Fastq file. All the samples will be grouped in one file.

---

## Join Paired Ends (if the reads are not already paired)**

*Note: using “Author_Year” and “Author” generically to explain our nomenclature used for each study. Anywhere these terms appear, the code is run once for each study.*

**Rename Files so they can be de-interleaved (i.e separate the forward and reverse reads into 2 distinct files). In the Terminal:**<br>

```bash
sed '/^@/ s/\./_/' Author_Year_sra_data.fastq > Author2.fastq
```
```bash
sed '/^+/ s/\./_/' Author2.fastq > Author_raw_renamed.fastq 
```

**Extract Forward and Reverse reads, common in SRA datasets that they are deposited together in the same file (not necessary for unpublished). In Qiime 1:**<br>

```bash
extract_reads_from_interleaved_file.py -i Author_raw_renamed.fastq -o Author_ForRev --forward_read_identifier .1 --reverse_read_identifier .2
```

**Join Paired Ends from what was just de-interleaved, here in Qiime1 (possible to join PE in Qiime 2 using “qiime vsearch join-pairs” function):**<br>

```bash
join_paired_ends.py -f Author_ForRev/forward_reads.fastq -r Author_ForRev/reverse_reads.fastq -o Author_PE
```
---

**Trimming the 515F Forward primers to just retain a specific region of the V4 hyervariable region of the 16S rRNA gene, in terminal:**<br>

(install cutadapt with: pip install --user --upgrade cutadapt)
~/.local/bin/cutadapt -g GTGYCAGCMGCCGCGGTAA Author_PE.fastq > Author_trimmed.fastq

---

## Split Fastq by sample for import to Qiime2. Done in Qiime1, here:

split_sequence_file_on_sample_ids.py -i Author_PE/fastqjoin.join.fastq -o Author_q2_importable --file_type fastq

---

## Trim 250 bp out from 515F primer in R package “DADA2”

Author_trimmed <- file.path(filt_path, paste0(sample.names, "Author_trimmed.fastq.gz")) 
Author_RevTrim <- filterAndTrim(fwd=fnFs,filt=Author_trimmed,truncLen=250,truncQ=0,trimLef=0,maxLen=Inf,minLen=0,maxN=0,minQ=0,maxEE=Inf) 

--

## DADA2 to identify ESVs. We ran this analysis on each independent study separately and merged the ESV tables and representative sequences afterwards (see merging steps below in Qiime 2)

### Infer Dada2 error model 
errFs = learnErrors(Author_RevTrim, multithread=TRUE) 

### Plot error model 
p = plotErrors(errFs, nominalQ=TRUE) 
ggsave("dada_errors_F.png", plot=p) 

### Dereplicate reads 
derepFs = derepFastq(Author_RevTrim, verbose=TRUE) 

### Infer samples
dadaFs = dada(derepFs, err=errFs, multithread=TRUE) 

### Make tables 
seqtab = makeSequenceTable(dadaFs) 
saveRDS(seqtab,"seqtab.rds”) 

### Remove chimera:
seqtab.nochim = removeBimeraDenovo(seqtab, method='consensus', verbose=TRUE) 
dim(seqtab.nochim) 
sum(seqtab.nochim)/sum(seqtab) 
saveRDS(seqtab.nochim,"seqtab.nochim.rds”) 

### Convert to biom file in R
library(devtools) 
library(phyloseq) 
library(biomformat) 
library(ggplot2) 

seqtab = readRDS(‘seqtab.rds') 
seqtab.nochim = readRDS('seqtab.nochim.rds’) 
taxtab = readRDS(‘taxtab.rds’) 
ps = phyloseq(otu_table(seqtab.nochim,taxa_are_rows=FALSE),tax_table(taxtab))  
phyloseq-class experiment-level object 

stn.biom = make_biom(t(seqtab)) 
write_biom(stn.biom, './seqtab.nochim.biom’) 
 
---

## In Qiime2 merging biom files and representative ESVs from each individually processed dataset:

### Import biom for each study:
qiime tools import \
  --input-path Author_dada2.biom \
  --type 'FeatureTable[Frequency]' \
  --source-format BIOMV210Format \
  --output-path Author_dada2.qza
…same for each study…

### Merge biom files into one:
qiime feature-table merge \
  --i-tables Ernakovich_dada2.qza \
  --i-tables Jurburg_dada2.qza \
  --i-tables Zhang_dada2.qza \
  --i-tables Fuentes_dada2.qza \
  --i-tables Lin_dada2.qza \
  --i-tables Nunes_dada2.qza \
  --i-tables Sun_dada2.qza \
  --i-tables Zhai_dada2.qza \
  --i-tables Simonin_dada2.qza \
  --o-merged-table MBSP_Dada2_merged.qza

### Import biom for each study:
qiime tools import \
  --input-path Author.fna \
  --output-path Author_seqs.qza \
  --type 'FeatureData[Sequence]'
…same for each study…

### Merge ESVs into one representative ESV fasta file.
qiime feature-table merge-seqs \
  --i-data Ernakovich_seqs.qza \
  --i-data Jurburg_seqs.qza \
  --i-data Zhang_seqs.qza \
  --i-data Fuentes_seqs.qza \
  --i-data Lin_seqs.qza \
  --i-data Nunes_seqs.qza \
  --i-data Sun_seqs.qza \
  --i-data Zhai_seqs.qza \
  --i-data Simonin_seqs.qza \
  --o-merged-data MBSP_repseqs.qza

---

## All archaeal, fungal, chloroplast, and mitochondrial ESVs identified with assign_taxonomy against the SILVA128 99% taxonomy (Qiime2 did not have classify-consensus with Silva at the time we processed the data). These ESVs were removed from the merged biom file, and the tree file (below), aside from several archaea as outgroups for rooting the tree.

assign_taxonomy.py -i MBSP.fasta -t consensus_taxonomy_7_levels.txt -r 99_otus_16S.fasta

#Filtering merged biom file:
filter_otus_from_otu_table.py -i Dada2_merged_table/feature-table.biom -o MBSP_NoUnassigned.biom -e Dada2_seqs/UnassignedTipNames.txt
filter_otus_from_otu_table.py -i MBSP_NoUnassigned.biom -o MBSP_NoChloroplast.biom -e Dada2_seqs/ChloroplastTipNames.txt
filter_otus_from_otu_table.py -i MBSP_NoChloroplast.biom -o MBSP_BactOnly.biom -e Dada2_seqs/ArchaealTipNames.txt

#########
### Qiime1 generate alignment using SILVA128 reference.

align_seqs.py -i dna-sequences.fasta -m pynast -t 99_otus_aligned.fasta -o MBSP_Qiime1_alignment

#filter alignment
filter_alignment.py -i MBSP_Qiime1_alignment/dna-sequences_aligned.fasta -o Filtered -e 0.10 -g 0.80

#########
### PASTA for improved alignment and building tree

#GUI interface (detailed software instructions available here: https://github.com/smirarab/pasta) but these are the settings used for the Microbiome Stress Project tree and improved alignment:
Aligner: MAFFT
Merger: Muscle
Tree Estimator: Fasttree
Model GTR+G20

Output: MBSP.tre
