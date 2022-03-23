# 02. DADA2 - Prepare sequencing files
# Following tutorial https://benjjneb.github.io/dada2/tutorial.html
# Table S6
# Ref for DADA2: Callahan, B. J. et al. DADA2: High-resolution sample inference from Illumina amplicon data. Nat. Methods 13, 581–583 (2016).
# Ref for tab2fasta.py: https://github.com/dariober/bioinformatics-cafe
# Ref for QIIME2: Bolyen, E. et al. Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. Nat. Biotechnol. 37, 852–857 (2019).
# Ref for phyloseq: McMurdie, P. J. & Holmes, S. phyloseq: An R Package for Reproducible Interactive Analysis and Graphics of Microbiome Census Data. PLOS ONE 8, e61217 (2013).
# Water and sediment sequences were analyzed separately; this script also was repeated for each sequence batch

#### Prepare files ####
setwd(path)
fns <- list.files(input)
fns

fastqs <- fns[grepl(".fastq$", fns)]
fastqs <- sort(fastqs) # Sort ensures forward/reverse reads are in same order
fnFs <- fastqs[grepl("_R1", fastqs)] # Just the forward read files
fnRs <- fastqs[grepl("_R2", fastqs)] # Just the reverse read files

# Get sample names from the first part of the forward read filenames
sample.names <- sapply(strsplit(fnFs, "_R"), `[`, 1)

# Fully specify the path for the fnFs and fnRs
fnFs <- file.path(input, fnFs)
fnRs <- file.path(input, fnRs)

#### Filtering and Trimming ####
cairo_pdf(file.path(quality, "Rplot.qual_profile.raw.pdf"), onefile = T)
for(i in seq_along(fnFs)) {
    print(plotQualityProfile(fnFs[[i]]))
    print(plotQualityProfile(fnRs[[i]]))
}
dev.off()

## Change these values based on quality profile ##
# SeqBatch1
truncLen <- 240
trimLeft <- 10
trimRight <- 0

# SeqBatch2
truncLen <- 215
trimLeft <- 10
trimRight <- 0

# SeqBatch3
truncLen <- 220
trimLeft <- 0
trimRight <- 0

# SeqBatch4
truncLen <- 220
trimLeft <- 10
trimRight <- 0

# Sediment
truncLen <- 220
trimLeft <- 10
trimRight <- 0
##################################################

n <- 0
ee <- 2
lowcomplex <- 8
maxlength <- Inf
str_trimLeft <- paste(trimLeft, collapse = '_')
str_truncLen <- paste(truncLen, collapse = '_')
str_ee <- paste(ee, collapse = '_')
base <- paste("filt_trim", ".tlen", str_truncLen, ".tleft", str_trimLeft, ".ee", str_ee, ".n", n, sep = "")
filt_path <- file.path(filter, base)
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt_trim.fastq"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt_trim.fastq"))

filtered_out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = truncLen, trimLeft = trimLeft,
maxN=n, maxEE=ee, truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE, rm.lowcomplex = lowcomplex)

filtered_out

cairo_pdf(file.path(quality, "Rplot.qual_profile.filt_trim.pdf"), onefile = T)
for(i in seq_along(fnFs)) {
    print(plotQualityProfile(filtFs[[i]]))
    print(plotQualityProfile(filtRs[[i]]))
}
dev.off()

#### Learn Error Rates ####
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

cairo_pdf(file.path(error, "Rplot.estimated_error_rates.pdf"), onefile = T)
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
dev.off()

cairo_pdf(file.path(error, "Rplot.error_profileFor.filt_trim.pdf"), onefile = T)
plot(colSums(errF$trans))
dev.off()

cairo_pdf(file.path(error, "Rplot.error_profileRev.filt_trim.pdf"), onefile = T)
plot(colSums(errR$trans))
dev.off()

#### Dereplication ####
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

f_summary_derep <- file.path(dereplication, "summary.derepFastq.txt")
for(i in seq_along(fnFs)) {
    strF <- summary(derepFs[[i]])
    strR <- summary(derepRs[[i]])
    write(strF, file = f_summary_derep, append = T)
    write(strR, file = f_summary_derep, append = T)
}

#### Sample Inference ASVs ####
dadaFs <- dada(derepFs, err=errF, multithread=TRUE, HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE, HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32)

#### Merge ####
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
head(mergers[[1]])

#sample.names
for(i in seq_along(fnFs)) {
    sample <- sample.names[i]
    names <- rownames(mergers[[i]])
    mergers1 <- cbind(name = c(names), mergers[[i]][c(2:length(mergers[[i]][1,]),1)])
    f_mergers <- file.path(merging, paste("mergers", ".", sample, ".txt", sep = ""))
    write.table(mergers1, file = f_mergers, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}

#### Construct sequence table ####
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

table(nchar(getSequences(seqtab)))

# Remove sequences that are too long or too short
seqtab1 <- seqtab[,nchar(colnames(seqtab)) %in% seq(227,237)] #Batch1: 262-267; Batch2: 262-266; Batch3: 251-257; Batch4: 227-237; Sediment: 230-237
table(nchar(getSequences(seqtab1)))

seqtab <- seqtab1

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

# percent of reads that did not have chimeras. (Note: 100-%% = % of reads that were chimeras)
sum(seqtab.nochim)/sum(seqtab)

# Inspect distribution of sequence lengths
table(nchar(colnames(seqtab.nochim)))

#### Create data tables ####
getN <- function(x) sum(getUniques(x))
track <- cbind(filtered_out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

# create txt file of track table
dir.create(paste0(path,"/",batch,"_dada2_data"))
track_table <- file.path(paste0(path,"/",batch,"_dada2_data"), "summary_reads_table.txt")
write.table(track, file = track_table, append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

# Summary of # sequences, chimeras, paired reads (Txt file created)
f_summary_seqtab <- file.path(paste0(path,"/",batch,"_dada2_data"), "summary.seqtabh.txt")
es <- c()
e1 <- c(name = 'seqtab', nseq = sum(seqtab), nrep = dim(seqtab)[2])
e2 <- c(name = 'seqtab.nochim', nseq = sum(seqtab.nochim), nrep = dim(seqtab.nochim)[2])
es <- rbind(es, e1,e2)
write.table(es, file = f_summary_seqtab, append = F, quote = F, sep = "\t", row.names = F, col.names = T)

# Summary of Sequence length (Txt file created)
f_summary_seqlen <- file.path(paste0(path,"/",batch,"_dada2_data"), "summary.seqlength.txt")
con21 <- table(nchar(colnames(seqtab)))
es21 <- as.data.frame(con21)

con22 <- table(nchar(colnames(seqtab.nochim)))
es22 <- as.data.frame(con22)
es2 <- plyr::join(es21,es22,by=c('Var1'))
colnames(es2) <- c('length','seqtab','seqtab.nochim')
write.table(es2, file = f_summary_seqlen, append = F, quote = F, sep = "\t", row.names = F, col.names = T)

# Transform seqtab_nochim
t.seqtab <- t(seqtab)
t.seqtab1 <- as.data.frame(cbind(t.seqtab, seq = rownames(t.seqtab)))
f_seqtab <- file.path(paste0(path,"/",batch,"_dada2_data"), "seqtab_nochim.txt")
write.table(t.seqtab1, file = f_seqtab, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# Transform seqtab_nochim.txt
t.seqtab.nochim <- t(seqtab.nochim)
t.seqtab.nochim1 <- as.data.frame(cbind(t.seqtab.nochim, seq = rownames(t.seqtab.nochim)))
f_seqtab_nochim <- file.path(paste0(path,"/",batch,"_dada2_data"), "seqtab_nochim.txt")
write.table(t.seqtab.nochim1, file = f_seqtab_nochim, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# create fasta file
f_seqtab_nochim_fa <- file.path(paste0(path,"/",batch,"_dada2_data"), "seqtab_nochim.uchime.fa")
uniquesToFasta(seqtab.nochim, f_seqtab_nochim_fa, ids = NULL, mode = "w")

# Label and number each "ASV"
labels <- paste('ASV', seq(0,length(t.seqtab.nochim[,1])-1), sep = "_")
labels2 <- paste(seq(0,length(t.seqtab.nochim[,1])-1))

# create fasta file with labels
f_seqtab_nochim_fa <- file.path(paste0(path,"/",batch,"_dada2_data"), "seqtab_nochim.fa")
uniquesToFasta(seqtab.nochim, f_seqtab_nochim_fa, ids = labels, mode = "w")

# create txt file of seqtab_nochim
f_seqtab_nochim_asis <- file.path(paste0(path,"/",batch,"_dada2_data"), "seqtab_nochim.asis.txt")
write.table(seqtab.nochim, file = f_seqtab_nochim_asis, append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

# create txt file of seqtab_nochim as "ASV" table
t.seqtab.nochim.otu <- as.data.frame(cbind('names' = labels, t.seqtab.nochim))
f_seqtab_nochim_otu <- file.path(paste0(path,"/",batch,"_dada2_data"), "seqtab_nochim.fa.asv_table.txt")
write.table(t.seqtab.nochim.otu, file = f_seqtab_nochim_otu, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

#### Repeat the above commands for different sequence batches ####

#### Clustering ASVs by QIIME2 (using bash) ####
# Prepare the files
wd=<path to QIIME2 folder>
cd ${wd}/input

# Rename the seqtab_nochim.txt in the dada2 batch data folders, copy to folder for qiime
# Rename sequence identifiers, retrieve sequences, and combine to "asv_table_wSeq.txt" ; make separate file "asv_table.txt" which does not have sequences
# Manually combine ASV tables (add 0 for ASVs not in other samples)

input=asv_table_wSeq.txt
output=seqtab_nochim_concatenated.fa
python3 tab2fasta.py ${input} 2 1 > ${output} #tab2fasta.py: https://github.com/dariober/bioinformatics-cafe
sed -i -e '1,2d' ${output} #delete first two lines from file
sed 's/_NA//g' ${output} > ${output}2 #delete all the NAs
rm ${output}-e ${output}
mv ${output}2 ${output}

conda activate qiime2-2021.4

# import sequences to qiime
input=${wd}/input/seqtab_nochim_concatenated.fa
output=${wd}/input/seqtab_nochim_concatenated.qza

qiime tools import \
  --input-path ${input} \
  --output-path ${output} \
  --type 'FeatureData[Sequence]'

# import otu table to qiime
biom convert -i asv_table.txt -o asv_biom.biom --table-type="OTU table" --to-hdf5

qiime tools import \
--input-path asv_biom.biom \
--type 'FeatureTable[Frequency]' \
--input-format BIOMV210Format \
--output-path feature-table.qza

# Cluster the sequences at 0.97 level
table=${wd}/qiime2_clustering/feature-table.qza
input=${wd}/qiime2_clustering/seqtab_nochim.qza
i=0.97

qiime vsearch cluster-features-de-novo \
  --i-table ${table} \
  --i-sequences ${input} \
  --p-perc-identity ${i} \
  --o-clustered-table ${wd}/qiime2_clustering/table-dn-${i}.qza \
  --o-clustered-sequences ${wd}/qiime2_clustering/rep-seqs-dn-${i}.qza

qiime feature-table summarize \
    --i-table ${wd}/qiime2_clustering/table-dn-${i}.qza \
    --o-visualization ${wd}/qiime2_clustering/table-dn-${i}-summary.qzv

qiime tools view ${wd}/qiime2_clustering/table-dn-${i}-summary.qzv

# export the table
qiime tools export \
  --input-path table-dn-${i}.qza \
  --output-path table-dn-${i}.biom

biom convert -i table-dn-${i}.biom/feature-table.biom -o table-dn-${i}.tsv --to-tsv

# export sequences
qiime tools export --input-path ${wd}/rep-seqs-dn-${i}.qza  --output-path ${wd}/export-${i}

# convert fasta to tabulated
input=${wd}/export-${i}/dna-sequences.fasta
output=${wd}/export-${i}/dna-sequences.txt

seqkit fx2tab ${input} > ${output}

#### add to otu table "table-dn-0.97.xlsx" ####

#### Assign Taxonomy using with QIIME2 ####
classifier=<path to silva-138-99-515-806-nb-classifier.qza>
reads=rep-seqs-dn-${i}.qza
output=rep-seqs-dn-${i}_silva138-99.qza
viz_out=rep-seqs-dn-${i}_silva138-99.qzv

qiime feature-classifier classify-sklearn \
  --i-classifier ${classifier} \
  --i-reads ${reads} \
  --o-classification ${output}

qiime metadata tabulate \
  --m-input-file ${output} \
  --o-visualization ${viz_out}

qiime tools export \
    --input-path rep-seqs-dn-0.97_silva138-99.qza \
    --output-path rep-seqs-dn-0.97_silva138-99

#### Import QIIME2 data to R ####
library(qiime2R)
library(phyloseq)

qiime_path <- "<path to qiime working dir>"
setwd(qiime_path)

### metadata
library(readxl)
DATA_PHYLOSEQ <- read_excel("<path to metadata.xlsx>", na = "NA")
DATA_PHYLOSEQ_FIXED <- data.frame(DATA_PHYLOSEQ, row.names = 1)
head(DATA_PHYLOSEQ_FIXED)

### taxa
taxonomy <- read_qza(paste(qiime_path,"/rep-seqs-dn-0.97_silva138-99.qza",sep=""))
taxonomy <- as.matrix(parse_taxonomy(taxonomy$data))
head(taxonomy)

### features
features <- read_qza(paste(qiime_path,"/table-dn-0.97.qza",sep=""))
head(features)

write.table(features,paste(qiime_path,"feature_table_0.97.csv", sep=""), sep=",", row.names = T)
write.table(taxonomy,paste(qiime_path,"taxonomy_table_0.97.csv", sep=""), sep=",", row.names = T)

#### Change ASV to OTU, combine feature table, taxonomy table, and sequences together ####

#### import new taxonomy and otu table ####
combo_table <- read_excel("<path to combined table.xlsx>", sheet=1, na = "NA")
combo_table.df <- as.data.frame(as.matrix(combo_table))
tax_table <- combo_table.df[,1:10]
rownames(tax_table) <- tax_table$OTU_ID
tax_table$OTU_ID <- NULL
tax_table$ASV_ID <- NULL
tax_table$Seq <- NULL
tax_table$Kingdom <- gsub("d__", "", tax_table$Kingdom)
tax_table <- as.matrix(tax_table)
head(tax_table)

otu_table <- combo_table.df[, -c(1:10)]
otu_table <- as.data.frame(sapply(otu_table, as.numeric))
rownames(otu_table) <- combo_table.df$OTU_ID
otu_table <- as.matrix(otu_table)
head(otu_table)

### Make Phyloseq object
ps_use <- phyloseq(
    otu_table(otu_table, taxa_are_rows = TRUE),
    sample_data(DATA_PHYLOSEQ_FIXED),
    tax_table(tax_table)
)

ps_use
