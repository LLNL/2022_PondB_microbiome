# 03. Filtering OTUs
# Table S6, Table S9, Table S12, Figure S1
# Ref for checking control samples: https://benjjneb.github.io/decontam/vignettes/decontam_intro.html
# Ref for prevalence checking: https://ucdavis-bioinformatics-training.github.io/2017-September-Microbial-Community-Analysis-Workshop/friday/MCA_Workshop_R/phyloseq.html)
# Ref for rarefaction curve: ranacapa v0.1.0: Kandlikar GS, Gold ZJ, Cowen MC et al. ranacapa: An R package and Shiny web app to explore environmental DNA data with exploratory statistics and interactive visualizations [version 1; peer review: 1 approved, 2 approved with reservations]. F1000Research 2018, 7:1734 (https://doi.org/10.12688/f1000research.16680.1)
# Ref for phyloseq: McMurdie, P. J. & Holmes, S. phyloseq: An R Package for Reproducible Interactive Analysis and Graphics of Microbiome Census Data. PLOS ONE 8, e61217 (2013).
# Water and sediment sequences were analyzed separately

#### Phyloseq object ####
ps_use

#### Check library sizes ####
df <- as.data.frame(sample_data(ps_use)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(ps_use)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_Control, shape=Sequence_batch)) + geom_point()

#### Check OTUs in Control samples ####
sample_data(ps_use)$is.neg <- sample_data(ps_use)$Sample_Control == "Control"
contamdf.prev <- isContaminant(ps_use, conc="DNA_conc_ng_uL", method="combined", neg="is.neg", threshold=0.1, batch="Sequence_batch")

cairo_pdf(file.path(path_phy, "threshold_remove_contaminants_decontam.pdf"), onefile = T)
hist(contamdf.prev$p, 100, xlim = c(0,1))
dev.off()

table(contamdf.prev$contaminant)

# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa <- transform_sample_counts(ps_use, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_Control == "Control", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_Control == "Sample", ps.pa)

write.table(contamdf.prev97,paste("contamdf.prev.freq.csv", sep=""), sep=",", row.names = T)

#### Manually inspect the table; assign OTUs as contaminants if are Eukaryota, Mitochondria, Chloroplast, Kingdom/Phyla without assignment; Control samples seem to not impact Pond B water/sediment samples so will not do decontam filtering, but following the decontam code to remove Eukaryota, Mitochondria, Chloroplast, Kingdom/Phyla without assignment ####

#### Load in new table and remove contaminants ####
fix_cont <- read.csv(paste0(path_phy,"/contamdf.prev.freq.csv"), row.names=1)
fix_cont$contaminant <- as.logical(fix_cont$contaminant_manual)
head(fix_cont)
nrow(fix_cont)
nrow(contamdf.prev)

head(fix_cont[ which(fix_cont$contaminant_manual=='TRUE'), ])
table(fix_cont$contaminant_manual)

# Make new phyloseq object
ps.pa <- transform_sample_counts(ps_use, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_Control == "Control", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_Control == "Sample", ps.pa)

# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg), contaminant=fix_cont$contaminant_manual)

# prune out the contaminants
ps.noncontam <- prune_taxa(!fix_cont$contaminant_manual, ps_use)
ps_use
ps.noncontam

# remove control samples
ps.noncontam.sub = subset_samples(ps.noncontam, Sample_Control != "Control")

# remove samples with only 0 counts
ps.noncontam.sub <- prune_samples(sample_sums(ps.noncontam.sub) > 0, ps.noncontam.sub)
ps.noncontam.sub <- filter_taxa(ps.noncontam.sub, function(x) sum(x) > 0, TRUE)
ps.noncontam.sub

#### Examine the prevalence of each taxa ####
prev0 = apply(X = otu_table(ps.noncontam.sub),
        MARGIN = 1,
        FUN = function(x){sum(x > 0)})
        prevdf = data.frame(Prevalence = prev0,
        TotalAbundance = taxa_sums(ps.noncontam.sub),
        tax_table(ps.noncontam.sub))

prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(ps.noncontam.sub, taxonomic.rank = "Phylum"))

# Plot the prevalence
plot_prev <- ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps.noncontam.sub),color=Phylum)) +
                    geom_hline(yintercept = 0.02, alpha = 0.5, linetype = 2) +
                    geom_vline(xintercept = 5, alpha = 0.5, linetype = 2) +
                    geom_point(size = 2, alpha = 0.7) +
                    scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
                    facet_wrap(~Phylum) + theme(legend.position="none")

ggsave("prevalencethreshold.all.noncontam.pdf", path = path_phy, scale = 1, width = 15, height = 10, units = c("in"), dpi = 300)

# Zoom in to the plot
plot_prev + ylim(0,0.05) + xlim(0,20)
ggsave("prevalencethreshold.zoom.all.noncontam.pdf", path = path_phy, scale = 1, width = 15, height = 10, units = c("in"), dpi = 300)

#### Check the plots ####

# Keep phyla > 0
keepPhyla = table(prevdf$Phylum)[(table(prevdf$Phylum) > 0)]
prevdf1 = subset(prevdf, Phylum %in% names(keepPhyla))

# Track the kept phyla
prevdf1_track <- plyr::ddply(prevdf1, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
prevdf1_table <- file.path(path_phy, "phyloseq_prevalence_phylum_kept.txt")
write.table(prevdf1_track, file = prevdf1_table, append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

# Skip prevalence threshold
prevalenceThreshold = 0

# Execute prevalence filter, using `prune_taxa()` function
ps.noncontam.prev = prune_taxa((prev0 > prevalenceThreshold), ps.noncontam.sub)
ps.noncontam.prev

#### Plot the number of reads ####
sample_sum_df <- data.frame(sum = sample_sums(ps.noncontam.prev))

# Histogram of sample read counts
ggplot(sample_sum_df, aes(x = sum)) +
    geom_histogram(color = "black", fill = "indianred", binwidth = 2500) +
    ggtitle("Distribution of sample sequencing depth") +
    xlab("Read counts") +
    ylab("Number of samples") +
    #plot_nomargins_x +
    scale_y_continuous(expand = expansion(mult = c(0, .1))) +
theme(panel.background = element_rect(fill = "white", colour = "black", size = 1, linetype = "solid"),
    panel.border = element_rect(colour="black", size=1, fill=NA),
    panel.grid.major = element_line(size = 0),
    panel.grid.minor = element_line(size = 0),
    axis.text = element_text(size=10, colour="black"),
    axis.title = element_text(face="bold", size=10, colour="black"),
    plot.title = element_text(face="bold", size=12),
    legend.position="none")

ggsave("sample_sequencing_depth.pdf", path = path_phy, scale = 1, width = 5, height = 5, units = c("in"), dpi = 300)

#### Rarefaction curve: Figure S1 ####
ps.noncontam <- phyloseq(
    otu_table(otu_table(ps.noncontam.prev), taxa_are_rows = TRUE),
    sample_data(DATA_PHYLOSEQ_FIXED),
    tax_table(tax_table(ps.noncontam.prev))
)
ps.noncontam

rarefraction.curve <- ggrare(ps.noncontam, step = 500, se = TRUE, label="SampleID")

rarefraction.curve.plot <- rarefraction.curve +
facet_grid(Strat_group~Location_description, scales="free") +
scale_x_continuous(limits=c(0,500000), breaks=c(0,200000,400000)) +
scale_y_continuous(limits=c(0,3000), breaks=c(0,1000,2000,3000)) +
theme(panel.background = element_rect(fill = "white", colour = "black", size = 1, linetype = "solid"),
    panel.border = element_rect(colour="black", size=1, fill=NA),
    strip.background=element_rect(fill='white', colour='white'),
    strip.text = element_text(face="bold", size=10),
    panel.grid.major = element_line(size = 0),
    panel.grid.minor = element_line(size = 0),
    axis.text = element_text(face="bold", size=10, colour="black"),
    axis.title = element_text(face="bold", size=10),
    legend.position="bottom",
    legend.key = element_rect(fill = "white"),
    legend.title = element_blank()) +
guides(color="none")

rarefraction.curve.plot
ggsave("Rarefraction_curve.svg", path = path_phy, scale = 1, width = 10, height = 10, units = c("in"), dpi = 300)

#### Filter the samples ≥ 5 read counts ####

ps_filt5 <- phyloseq_filter_sample_wise_abund_trim(ps.noncontam, minabund = 5, rm_zero_OTUs = TRUE)
ps_filt5

write.table(tax_table(ps_filt5), file = "TAX_PHYLOSEQ_ps_filt5.txt", append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(otu_table(ps_filt5), file = "OTU_PHYLOSEQ_ps_filt5.txt", append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
