# 08.Songbird_OTUs.r
# Figure 4, Figure 7, Figure S7, Figure S8, Table S10, Table S11
# Ref for Songbird: Morton, J. T. et al. Establishing microbial composition measurement standards with reference frames. Nat. Commun. 10, 2719 (2019).
# Ref for Qurro: Fedarko, M. W. et al. Visualizing ’omic feature rankings and log-ratios using Qurro. NAR Genomics Bioinforma. 2, (2020).
# Water and sediment sequences were analyzed separately

#### Phyloseq Object ####
ps.noncontam.tree

#### Set Directory ####
songbird <- file.path(paste(path_phy, "Songbird", sep=""))
dir.create(songbird, showWarnings=FALSE)
setwd(songbird)

#### Prepare files for Songbird ####
# format taxonomy table
tax <- as(tax_table(ps.noncontam.tree),"matrix")
tax <- as.data.frame(tax)

tax$Kingdom <- paste("k", tax$Kingdom, sep="__")
tax$Phylum <- paste("p", tax$Phylum, sep="__")
tax$Class <- paste("c", tax$Class, sep="__")
tax$Order <- paste("o", tax$Order, sep="__")
tax$Family <- paste("f", tax$Family, sep="__")
tax$Genus <- paste("g", tax$Genus, sep="__")
tax$Species <- paste("s", tax$Species, sep="__")

tax_cols <- c("Kingdom", "Phylum", "Class","Order","Family","Genus","Species")
tax$taxonomy <- do.call(paste, c(tax[tax_cols], sep=";"))
for(co in tax_cols) tax[co]<-NULL

write.table(tax, "tax_for_qiime2.txt", quote=FALSE, col.names=FALSE, sep="\t")

# make a biomformat otu table
otu <- as(otu_table(ps.noncontam.tree),"matrix")
otu_biom <- make_biom(data=otu)
write_biom(otu_biom,"otu_biom.biom")
write.table(otu_table(ps.noncontam.tree), file = "otu_table.txt", sep = "\t", row.names = TRUE, col.names = NA)

# export metadata table
write.table(sample_data(ps.noncontam.tree), file = "metadata_for_qiime2.txt", sep = "\t", row.names = TRUE, col.names = NA)

#### Import data to QIIME2 (on bash) ####
conda activate qiime2-2020.6
wd=<path to Songbird working dir>
cd $wd

sed 's/"//g' metadata_for_qiime2.txt > metadata_for_qiime2_fixed.txt
# also add #SampleID to header
# add Train column for Songbird

biom convert -i otu_biom.biom -o otu_biom_HDF5.biom --to-hdf5

biom add-metadata -i otu_biom_HDF5.biom -o otu_wTax_metadata.biom --observation-metadata-fp tax_for_qiime2.txt --sc-separated taxonomy --observation-header OTUID,taxonomy --sample-metadata-fp metadata_for_qiime2_fixed.txt

# import to QIIME2
qiime tools import \
    --input-path otu_biom_HDF5.biom \
    --type 'FeatureTable[Frequency]' \
    --input-format BIOMV210Format \
    --output-path feature-table.qza

# import tax table to QIIME2
qiime tools import \
    --type 'FeatureData[Taxonomy]' \
    --input-format HeaderlessTSVTaxonomyFormat \
    --input-path tax_for_qiime2.txt \
    --output-path taxonomy.qza

# check import
qiime feature-table summarize \
    --i-table feature-table.qza \
    --m-sample-metadata-file metadata_for_qiime2_fixed.txt \
    --o-visualization summary_vis.qzv

qiime tools view summary_vis.qzv

#### Songbird ####

# Make the model; parameter of interest
dir=<name of folder for formula>
mkdir ${dir}

qiime songbird multinomial \
    --i-table feature-table.qza \
    --m-metadata-file metadata_for_qiime2_fixed.txt \
    --p-formula "<parameter of interest>" \
    --p-epochs 10000 \
    --p-differential-prior 0.5 \
    --p-summary-interval 1 \
    --p-num-random-test-examples 10 \ # For sediment samples: --p-num-random-test-examples 4
    --o-differentials ${dir}/differentials.qza \
    --o-regression-stats ${dir}/regression-stats.qza \
    --o-regression-biplot ${dir}/regression-biplot.qza \
    --p-training-column "Test_Train" \ #include if using
    --verbose

# Make the null model
null_dir=null_model
mkdir ${null_dir}

qiime songbird multinomial \
    --i-table feature-table.qza \
    --m-metadata-file metadata_for_qiime2_fixed.txt \
    --p-formula "1" \
    --p-epochs 10000 \
    --p-differential-prior 0.5 \
    --p-summary-interval 1 \
    --p-num-random-test-examples 10 \ # For sediment samples: --p-num-random-test-examples 4
    --o-differentials ${null_dir}/null-diff.qza \
    --o-regression-stats ${null_dir}/null-stats.qza \
    --o-regression-biplot ${null_dir}/null-biplot.qza \
    --p-training-column "Test_Train" \ #include if using
    --verbose

# Visualize the first model's regression stats and the null model
qiime songbird summarize-paired \
    --i-regression-stats ${dir}/regression-stats.qza \
    --i-baseline-stats ${null_dir}/null-stats.qza \
    --o-visualization ${dir}/paired-summary.qzv

qiime tools view ${dir}/paired-summary.qzv

# Plot the OTU rankings
qiime qurro differential-plot \
    --i-ranks ${dir}/differentials.qza \
    --i-table feature-table.qza \
    --m-sample-metadata-file metadata_for_qiime2_fixed.txt \
    --m-feature-metadata-file tax_for_qiime2.txt \
    --verbose \
    --o-visualization ${dir}/qurro_plot_q2.qzv

qiime tools view ${dir}/qurro_plot_q2.qzv

# Export the Songbird differentials
qiime metadata tabulate \
    --m-input-file ${dir}/differentials.qza \
    --o-visualization ${dir}/differentials-viz.qzv

qiime tools export \
  --input-path ${dir}/differentials-viz.qzv \
  --output-path ${dir}/differentials

#### Plot Songbird figures in R (Repeat for plotting differences between incipient stratification/stratification, Fe concentrations, and DO concentrations; and Sediments) ####
mypath = <path to Songbird data>
setwd(mypath)

# load bugs to focus on [tsv file with columns Family, Genus, Bug_type (iron oxidizer/reducer, methanogen, etc), Top_5perc_rank, Include (whether to include bug in plots), color (using hex color code), Top_10perc_rank]
bug_list <- read.delim(<path to selected bugs>, sep = "\t", stringsAsFactors=FALSE, fileEncoding="latin1")
bug_list$Family_Genus <- paste(bug_list$Family, bug_list$Genus) # use Family_Genus to match tax IDs; can change this to other combination for matching
rownames(bug_list) <- bug_list$Family_Genus
bug_list <- bug_list[ which(bug_list$Include=='Y'), ]
bug_list

# load otu/tax table
otu_tax_table <- read_excel("<path to file with OTU table and tax ID>")
otu_tax_table$Family_Genus <- paste(otu_tax_table$Family, otu_tax_table$Genus)
head(otu_tax_table)
write.csv(otu_tax_table, "otu_tax_table_ps_filt5.csv")

# subset OTUs from main otu/tax table (numerator)
focus_bugs <- otu_tax_table %>% filter(otu_tax_table$Family_Genus %in% bug_list$Family_Genus)
unique(focus_bugs$Family_Genus)
nrow(focus_bugs)
nrow(bug_list)

# get columns of interest only (Family_Genus and counts for all samples)
focus_bugs <- focus_bugs[,c(92, 11:91)]

# sum counts for each Family_Genus
focus_bugs <- focus_bugs %>%
      group_by(Family_Genus) %>%
      summarise_if(is.numeric, funs(sum))

# organize
focus_bugs <- as.data.frame(focus_bugs)
rownames(focus_bugs) <- focus_bugs$Family_Genus
focus_bugs[,1] <- NULL
focus_bugs <- as.data.frame(t(focus_bugs))
focus_bugs <- tibble::rownames_to_column(focus_bugs, "SampleID")

# setup denominator
bug_denominator <- otu_tax_table %>% filter(!otu_tax_table$Family_Genus %in% bug_list$Family_Genus)

# sum denominator
bug_sum_denominator_num <- bug_denominator %>% dplyr::select(where(is.numeric))
bug_sum_denominator_num <- as.data.frame(colSums(bug_sum_denominator_num))
bug_sum_denominator_num <- tibble::rownames_to_column(bug_sum_denominator_num, "SampleID")
head(bug_sum_denominator_num)

# calculate natural log ratios (Note: log in R is natural log)
df1_bugs <- focus_bugs
df2_bugs <- bug_sum_denominator_num
df3_bugs <- cbind(df1_bugs[1], log(df1_bugs[, -1] / df2_bugs[match(df1_bugs$SampleID, df2_bugs$SampleID), -1]))
head(df3_bugs)

# add descriptions
rownames(df3_bugs) <- df3_bugs$SampleID
df3_bugs$SampleID <- NULL
df3_bugs <- as.data.frame(t(df3_bugs))
df3_bugs$Family_Genus <- rownames(df3_bugs)
df3_bugs$Bug_type <- bug_list$Bug_type[match(df3_bugs$Family_Genus, bug_list$Family_Genus)]
df3_bugs$Top_5perc_rank <- bug_list$Top_5perc_rank[match(df3_bugs$Family_Genus, bug_list$Family_Genus)]
df3_bugs$color <- bug_list$color[match(df3_bugs$Family_Genus, bug_list$Family_Genus)]
head(df3_bugs)

# add stratificaiton group layer
df3_bugs_melt <- reshape2::melt(df3_bugs, id=c("Bug_type", "Top_5perc_rank", "Family_Genus", "color"))
colnames(df3_bugs_melt) <- c("Bug_type", "Top_5perc_rank", "Family_Genus", "color", "SampleID", "Natural_Log_Ratio")
df3_bugs_melt$SampleID <- gsub("\\.", "-", df3_bugs_melt$SampleID)
df3_bugs_melt$Strat_group <- DATA_PHYLOSEQ_FIXED$Strat_group[match(df3_bugs_melt$SampleID, DATA_PHYLOSEQ_FIXED$SampleID)]
df3_bugs_melt$strat_time <- ordi_data$strat_time[match(df3_bugs_melt$SampleID, ordi_data$SampleID)]
df3_bugs_melt$Depth_meters <- ordi_data$Depth_meters[match(df3_bugs_melt$SampleID, ordi_data$SampleID)]

# order the strat group
df3_bugs_melt$Strat_group <- gsub('Unstratified', 'Incipient Stratification', df3_bugs_melt$Strat_group)
df3_bugs_melt$Strat_group <- factor(df3_bugs_melt$Strat_group, levels=c("Incipient Stratification","Epilimnion","Thermocline", "Hypolimnion"), ordered = TRUE)

# remove Inf values
df3_bugs_melt <- df3_bugs_melt[is.finite(df3_bugs_melt$Natural_Log_Ratio),]

head(df3_bugs_melt)
write.csv(df3_bugs_melt, "Select_bugs_natural_log_ratio_v3.csv")

#### Plot Figure 7 ####

# Subset strat 3 months
df3_bugs_melt_sub <- subset(df3_bugs_melt, strat_time == "Stratified (~3 months)")

# Order the Bug_type
df3_bugs_melt_sub$Bug_type <- factor(df3_bugs_melt_sub$Bug_type, levels=c("Magnetotactic", "Iron oxidizer", "Iron oxidizer/reducer", "Photoferrotroph/iron reducer", "Iron reducer", "Sulfate/iron reducer", "Sulfate reducer", "Sulfide oxidizer", "Methylotroph", "Methanogen"), ordered=TRUE)

# Set up theme
plot_theme <- theme(panel.background = element_rect(fill = "white", colour = "black", size = 1, linetype = "solid"),
    panel.border = element_rect(colour="black", size=1, fill=NA),
    strip.background=element_rect(fill='white', colour='white', size = 0),
    strip.text = element_text(face="bold", size=15),
    panel.spacing.x=unit(1, "lines"),
    panel.spacing.y=unit(1, "lines"),
    panel.grid.major = element_line(size = 0),
    panel.grid.minor = element_line(size = 0),
    axis.text = element_text(size=15, colour="black"),
    axis.title = element_text(face="bold", size=15),
    legend.position="right",
    legend.key = element_rect(fill = "transparent"),
    legend.title = element_text(face="bold", size=15),
    legend.text = element_text(size=15, colour="black"),
    legend.background=element_blank())
plot_guide <- guides(
    fill = guide_legend(ncol=1, title="Taxa (Family_Genus)", override.aes=list(color="black", linetype=0)),
    color="none",
    linetype = "none")

# Set up coloring
colors <- distinct(df3_bugs_melt_sub, Family_Genus, color)
pal <- colors$color
names(pal) <- colors$Family_Genus
pal

figure_strat3months <- ggplot(df3_bugs_melt_sub, aes(x=Depth_meters, y=Natural_Log_Ratio)) +
    geom_point(aes(fill=Family_Genus), shape=21, colour="black", size=3) +
    geom_smooth(method = 'loess', size=1, alpha=0.25, linetype="dashed", aes(color=Family_Genus, fill=Family_Genus)) +
    facet_wrap(. ~ Bug_type, scales="free", nrow=2, labeller = labeller(Bug_type = label_wrap_gen(10))) +
    geom_vline(xintercept=2, linetype="dashed", color = "black", size=0.5) +
    geom_vline(xintercept=6.5, linetype="dashed", color = "black", size=0.5) +
    geom_vline(xintercept=0, linetype="solid", color = "black", size=0.5) +
    labs(x = "Depth (m)", y="Natural Log Ratio") +
    scale_y_continuous(lim=c(-12,0), breaks = seq(-12, 0, by = 4), expand = expansion(mult = c(0, 0))) +
    coord_flip() +
    scale_x_reverse(lim=c(10,-0.75), breaks = seq(10, 0, by = -2), expand = expansion(mult = c(0, 0))) +
    scale_fill_manual(values=pal) +
    scale_color_manual(values=pal) +
    plot_theme + plot_guide

figure_strat3months
ggsave("Selected_Bugs_strat3months.pdf", path = mypath, scale = 1, width = 18, height = 9, units = c("in"), dpi = 300)

#### For Fe combined plot (Figure S7) ####
Both_HighFe <- plot_grid(figure_Both + theme(legend.position="none", axis.title=element_blank()),
                 figure_HighFe + theme(legend.position="none", axis.title=element_blank()),
                 ncol=2, align = "v", axis="b")

Fe_hypo_ratio_All <- plot_grid(Both_HighFe,NULL,
                 figure_MidFe + theme(legend.position="none", axis.title=element_blank()),
                 nrow=3, align = "v", axis="b", rel_heights=c(1,0.1,1))

Fe_hypo_ratio_All

save_file <- paste("Combo_Fe_log_ratio.pdf", sep="")
ggsave(save_file, path = mypath, plot = Fe_hypo_ratio_All, scale = 1, width = 10, height = 8, units = c("in"), dpi = 300)

#### For O2 combined plot (Figure S8) ####
Both_O2 <- plot_grid(figure_obl + theme(legend.position="none", axis.title=element_blank()),
                 figure_O2 + theme(legend.position="none", axis.title=element_blank()),
                 ncol=2, align = "v", axis="b", rel_widths=c(1,0.75))

save_file <- paste("Combo_O2_log_ratio.pdf", sep="")
ggsave(save_file, path = mypath, plot = Both_O2, scale = 1, width = 15, height = 4, units = c("in"), dpi = 300)

#### For Fe and O2 log ratio plot ranking (Figure 4) ####
logratio <- read.delim("<path to data from Qurro>", sep = "\t", stringsAsFactors=FALSE, fileEncoding="latin1")
logratio$SampleID <- DATA_PHYLOSEQ_FIXED$SampleID[match(rownames(DATA_PHYLOSEQ_FIXED),logratio$Sample.ID)]
logratio$Season <- factor(logratio$Season, ordered=TRUE, levels=c("Spring", "Summer", "Fall"))
head(logratio)

# Set up theme
plot_theme <- theme(panel.background = element_rect(fill = "white", colour = "black", size = 1, linetype = "solid"),
    panel.border = element_rect(colour="black", size=1, fill=NA),
    panel.grid.major = element_line(size = 0),
    panel.grid.minor = element_line(size = 0),
    axis.text = element_text(size=15, colour="black"),
    axis.title = element_text(face="bold", size=20),
    legend.position="right",
    legend.key = element_rect(fill = "transparent"),
    legend.title = element_text(face="bold", size=20),
    legend.text = element_text(size=15, colour="black"),
    legend.background=element_blank())
plot_guide <- guides(
    fill = guide_legend(ncol=1, title="strat_time", override.aes=list(color="black", linetype=0)),
    color="none",
    linetype = "none")

figure <- ggplot(logratio, aes(x=Depth_meters, y=Current_Natural_Log_Ratio)) +
    geom_point(aes(fill=Season), shape=21, colour="black", size=3) +
    geom_smooth(method = 'loess', size=0.5, alpha=0.25, aes(color=Season, fill=Season, linetype=Season)) +
    geom_vline(xintercept=2, linetype="dashed", color = "black", size=0.5) +
    geom_vline(xintercept=6.5, linetype="dashed", color = "black", size=0.5) +
    geom_vline(xintercept=0, linetype="solid", color = "black", size=0.5) +
    labs(x = "Water Column Depth (m)", y="Natural Log Ratio") +
    coord_flip() +
    scale_x_reverse(lim=c(10,-0.75), breaks = seq(10, 0, by = -2), expand = expansion(mult = c(0, 0))) +
    scale_fill_manual(name="Stratification?", values=c("#a10702","#688e26", "#e9c46a"), labels=c("Incipient Stratification", "Stratified (~3 months)", "Stratified (~7 months)")) +
    scale_color_manual(name="Stratification?", values=c("#a10702","#688e26", "#e9c46a"), labels=c("Incipient Stratification", "Stratified (~3 months)", "Stratified (~7 months)")) +
    scale_linetype_manual(name="Stratification?", values=c("dotted", "longdash", "dotdash"), labels=c("Incipient Stratification", "Stratified (~3 months)", "Stratified (~7 months)")) +
    plot_theme + plot_guide
figure
ggsave("Log_Ratio_Qurro.pdf", path = mypath, scale = 1, width = 6, height = 5, units = c("in"), dpi = 300)

#### Plot OTU ranking for Fe and DO (Figure 4) ####
## For Fe
Fe_ranks <- read_excel("<path to Fe differentials.xlsx>", sheet="metadata")
colnames(Fe_ranks) <- c('featureid','Kingdom','Phylum','Class','Order','Family','Genus','Species','Bug_type','Intercept',
                        'Mid_High_Fe', 'Mid_Low_Fe','Mid_SuperLow_Fe')

## For DO
colnames(DO_ranks) <- c('featureid','Intercept','Aerobic_Anoxic', 'Aerobic_LowO2')
DO_ranks$Aerobic_LowO2 <- as.numeric(DO_ranks$Aerobic_LowO2)

# Sort by ranking (Fe)
Fe_ranks <- Fe_ranks[order(Fe_ranks$Mid_High_Fe),]

# Sort by ranking (DO)
DO_ranks <- DO_ranks[order(DO_ranks$Aerobic_LowO2),]

# Number the ranks (Fe or DO)
Fe_ranks$rank <- seq.int(nrow(Fe_ranks))
DO_ranks$rank <- seq.int(nrow(DO_ranks))

# Set up theme
plot_theme <- theme(panel.background = element_rect(fill = "white", colour = "black", size = 1, linetype = "solid"),
    panel.border = element_rect(colour="black", size=1, fill=NA),
    panel.spacing.x=unit(1, "lines"),
    panel.spacing.y=unit(1, "lines"),
    panel.grid.major = element_line(size = 0),
    panel.grid.minor = element_line(size = 0),
    axis.text = element_text(size=15, colour="black"),
    axis.title = element_text(face="bold", size=20),
    axis.title.y = element_text(angle = 90))
plot_guide <- guides(fill= guide_legend(order=1, override.aes = list(size=0.2, color="black", alpha=1)),
                     color = "none", linetype = "none")
                                   
# Change y variables to DO or Fe
ggplot(DO_ranks, aes(x=rank, y=Aerobic_LowO2)) +
    geom_line(color="black") +
    geom_ribbon(data=subset(DO_ranks, rank>1877 & rank<2085),aes(ymax=Aerobic_LowO2),ymin=0, fill="#941b0c",colour=NA,alpha=1)+
    geom_ribbon(data=subset(DO_ranks, rank>1117 & rank<1878),aes(ymax=Aerobic_LowO2),ymin=0, fill="grey",colour=NA,alpha=0.5)+
    geom_ribbon(data=subset(DO_ranks, rank>208 & rank<1118),aes(ymax=Aerobic_LowO2),ymin=0, fill="grey",colour=NA,alpha=0.5)+
    geom_ribbon(data=subset(DO_ranks, rank>0 & rank<209),aes(ymax=Aerobic_LowO2),ymin=0, fill="#f6aa1c",colour=NA,alpha=1)+
    geom_hline(yintercept=0, linetype="longdash", color = "black", size=0.3) +
    #geom_area(mapping = aes(x = ifelse(rank>1877 & rank<2085, rank, 0)), fill = "red") +
    ylab("log(LowO2/Aerobic)") +
    xlab("Taxa Rank") +
    scale_y_continuous(limits=c(-5.5,5.5), breaks=c(-4,-2,0,2,4), expand = c(0, 0)) +
    scale_x_continuous(limits=c(1,2085), breaks=c(1,500,1000,1500,2000), expand = c(0, 0)) +
    plot_theme + plot_guide
ggsave("Log_Ratio_ranks.pdf", path = mypath, scale = 1, width = 4.5, height = 4, units = c("in"), dpi = 300)
