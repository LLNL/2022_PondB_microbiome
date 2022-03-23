# 07. Beta Diversity
# Figure 3c, Figure 5c, Figure S6, Figure S10
# Ref for CLR PCA: Gloor, G. B.; Macklaim, J. M.; Pawlowsky-Glahn, V.; Egozcue, J. J. Microbiome Datasets Are Compositional: And This Is Not Optional. Front. Microbiol. 2017, 8. https://doi.org/10.3389/fmicb.2017.02224.
# Ref for PhILR PCA: Silverman, J. D.; Washburne, A. D.; Mukherjee, S.; David, L. A. A Phylogenetic Transform Enhances Analysis of Compositional Microbiota Data. eLife 2017, 6, e21887. https://doi.org/10.7554/eLife.21887.
# Water and sediment sequences were analyzed separately

#### Phyloseq Object ####
ps.noncontam.tree

#### Set Directory ####
beta <- file.path(paste(path_phy, "/beta_diversity/Ordination", sep=""))
dir.create(beta, showWarnings=FALSE, recursive=TRUE)
setwd(beta)

#### CLR transformation ####
set.seed(2021)

# Obtain the OTU table
ps_filt_otu <- as.data.frame(otu_table(ps.noncontam.tree))
metadata <- as.data.frame(as.matrix(sample_data(ps.noncontam.tree)))

# Transform the OTU table and replace 0 values with an estimate
f.n0 <- cmultRepl(t(ps_filt_otu), method="CZM", label=0, output="p-counts")

# CLR Transformation
f.clr <- codaSeq.clr(f.n0, IQLR=FALSE)
head(f.clr)

# Check distribution
qplot(rowSums(f.clr),bins=50) + xlab("clr transformed counts-per-sample")
qplot(colSums(f.clr),bins=50) + xlab("clr transformed counts-per-taxa")

# Add CLR to phyloseq object
ps_clr <- phyloseq(
    otu_table(t(f.clr), taxa_are_rows = TRUE),
    sample_data(sample_data(ps.noncontam.tree)),
    phy_tree(filt_tree),
    tax_table(tax_table(ps.noncontam.tree))
)
ps_clr

# get PCA
f.pcx <- prcomp(f.clr)

# extract results from PCA
var <- get_pca_var(f.pcx)

# get eigenvalue
eig.val <- get_eigenvalue(f.pcx)
head(eig.val)

#### PhILR transformation ####
library("philr")

# Add clr to phyloseq object
ps_philr <- phyloseq(
    otu_table(t(f.n0), taxa_are_rows = TRUE),
    sample_data(sample_data(ps.noncontam.tree)),
    phy_tree(filt_tree),
    tax_table(tax_table(ps.noncontam.tree))
)
ps_philr

# check the object
is.rooted(phy_tree(ps_philr))
is.binary.tree(phy_tree(ps_philr))

# labels and name balance
phy_tree(ps_philr) <- makeNodeLabel(phy_tree(ps_philr), method="number", prefix='n')

name.balance(phy_tree(ps_philr), tax_table(ps_philr), 'n1')

# calculate philr
otu.table <- as.data.frame(otu_table(ps_philr))
otu.table <- as.matrix(t(otu.table))
tree <- phy_tree(ps_philr)
metadata <- sample_data(ps_philr)
tax <- tax_table(ps_philr)

gp.philr <- philr(otu.table, tree, part.weights='enorm.x.gm.counts', ilr.weights='blw.sqrt', return.all=TRUE)

# get philr info
gp.dist <- dist(gp.philr$df.ilrp, method="euclidean")
gp.pcoa <- ordinate(ps_philr, 'PCoA', distance=gp.dist)

#### Prepare for combining CLR and PhILR data to plot (Water Samples) ####
ordi_clr <- as.data.frame(f.pcx$x[,1:2])
ordi_clr$ID_keep <- rownames(ordi_clr)
head(ordi_clr)

DATA_PHYLOSEQ_FIXED$ID_keep <- rownames(DATA_PHYLOSEQ_FIXED)
ordi_data <- merge(DATA_PHYLOSEQ_FIXED, ordi_clr, by = c('ID_keep'))
ordi_data$strat_time <- ifelse(ordi_data$Season == "Fall", "Stratified (~7 months)",
                  ifelse(ordi_data$Season == "Spring", "Incipient Stratification",
                         ifelse(ordi_data$Season == "Summer", "Stratified (~3 months)", "")))
ordi_data$strat_time <- factor(ordi_data$strat_time, ordered = TRUE, levels = c("Incipient Stratification", "Stratified (~3 months)", "Stratified (~7 months)"))

ordi.scores_philr <- as.data.frame(gp.pcoa$vectors)
colnames(ordi.scores_philr) <- c("Axis.1philr", "Axis.2philr")
ordi.scores_philr <- ordi.scores_philr[c("Axis.1philr","Axis.2philr")]
ordi.scores_philr$ID_keep <- rownames(ordi.scores_philr)

# Merge data with ordi_data
ordi_data <- merge(ordi_data, ordi.scores_philr, by = c('ID_keep'))
head(ordi_data)
colnames(ordi_data)

# remove bad samples
ordi_data_sub <- subset(ordi_data, SampleID != "W64")
ordi_data_sub <- subset(ordi_data_sub, SampleID != "W66")
ordi_data_sub <- subset(ordi_data_sub, SampleID != "W67")

# remove batch 4 (same as batch 1 and 3)
ordi_data_sub_NoBatch4 <- subset(ordi_data_sub, Sequence_batch != "Batch4")

# data for plotting batches
ordi_data_sub_batch <- subset(ordi_data_sub, Plot_batch == "Y")

#### Prepare for combining CLR and PhILR data to plot (Sediment Samples) ####
# CLR
ordi_clr = ordinate(ps_clr, "PCoA", "euclidean")
head(ordi_clr$values)
head(ordi_clr$vectors)

ordi.scores_clr <- as.data.frame(ordi_clr$vectors)
ordi.scores_clr[,3:20] <- NULL
colnames(ordi.scores_clr) <- c("Axis.1clr", "Axis2.clr")
ordi.scores_clr$ID_keep <- rownames(ordi.scores_clr)

# PhILR
ordi.scores_philr <- as.data.frame(gp.pcoa$vectors)
ordi.scores_philr[,3:20] <- NULL
colnames(ordi.scores_philr) <- c("Axis.1philr", "Axis.2philr")
ordi.scores_philr$ID_keep <- rownames(ordi.scores_philr)

# Merge data with ordi_data
DATA_PHYLOSEQ_FIXED$ID_keep <- rownames(DATA_PHYLOSEQ_FIXED)
ordi_data <- merge(DATA_PHYLOSEQ_FIXED, ordi.scores_clr, by = c('ID_keep'))
ordi_data <- merge(ordi_data, ordi.scores_philr, by = c('ID_keep'))
head(ordi_data)

#### Plot CLR and PhILR PCA ####
# Set up theme
plot_theme <- theme(panel.background = element_rect(fill = "white", colour = "black", size = 1, linetype = "solid"),
    panel.border = element_rect(colour="black", size=1, fill=NA),
    strip.background=element_rect(fill='white', colour='white', size = 0),
    strip.text = element_text(face="bold", size=20),
    panel.spacing.x=unit(0.5, "lines"),
    panel.grid.major = element_line(size = 0),
    panel.grid.minor = element_line(size = 0),
    axis.text = element_text(size=15, colour="black"),
    axis.title = element_text(face="bold", size=20),
    legend.position="right",
    legend.key = element_rect(fill = "transparent"),
    legend.title = element_text(face="bold", size=20),
    legend.text = element_text(size=20),
    legend.background = element_rect(fill = "transparent"),
    legend.box.background = element_rect(fill = "transparent", colour = NA))
plot_guide <- guides(fill = guide_colourbar(frame.linewidth = 1, frame.colour = "black", ticks = TRUE, ticks.colour = "black", ticks.linewidth = 1, reverse=T),
    shape = guide_legend(order=1, override.aes = list(size = 5, color="black", alpha=1)),
    color = "none")

### CLR PCA
## Water Samples
clr_plot_final <- ggplot(ordi_data_sub_NoBatch4, aes(x=PC1, y = PC2)) +
    geom_hline(yintercept=0, size=.2, linetype = "dashed", color="black") + geom_vline(xintercept=0, size=.2, linetype = "dashed", color="black") +
    geom_point(color="black", size=4, aes(fill=Depth_meters, shape=strat_time)) +
    xlab(paste("PC1: 29.2%")) +
    ylab(paste("PC2: 15.3%")) +
    scale_shape_manual(name = "Stratification?", values=c(21, 23, 24)) +
    scale_fill_gradient2(name = "Depth (m)", breaks = c(0, 2, 6, 10), limits = c(0, 10), low="#fffab1", mid="#ed105f", high="#023e73", midpoint=5) +
     plot_theme + plot_guide
clr_plot_final

save_file <- paste("clr.pdf", sep="")
ggsave(save_file, path = beta, scale = 1, width = 8, height = 5, units = c("in"), dpi = 300)

# label CLR PCA
clr_plot_label <- clr_plot_final +
    geom_label_repel(aes(label = SampleID), box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50', max.overlaps=Inf)

save_file <- paste("clr_labelled.pdf", sep="")
ggsave(save_file, path = beta, scale = 1, width = 8, height = 5, units = c("in"), dpi = 300)

# Plot batches CLR PCA
plot_guide <- guides(shape = guide_legend(order=1, override.aes = list(size = 5, color="black", alpha=1)),
                     color = "none")
ordi_data_sub_batch$Sequence_batch <- gsub("Batch1", "Batch 1", ordi_data_sub_batch$Sequence_batch)
ordi_data_sub_batch$Sequence_batch <- gsub("Batch2", "Batch 2", ordi_data_sub_batch$Sequence_batch)
ordi_data_sub_batch$Sequence_batch <- gsub("Batch3", "Batch 3", ordi_data_sub_batch$Sequence_batch)
ordi_data_sub_batch$Sequence_batch <- gsub("Batch4", "Batch 4", ordi_data_sub_batch$Sequence_batch)

clr_plot_batch <- ggplot(ordi_data_sub_batch, aes(x=PC1, y = PC2)) +
    geom_hline(yintercept=0, size=.2, linetype = "dashed", color="black") + geom_vline(xintercept=0, size=.2, linetype = "dashed", color="black") +
    geom_point(color="black", size=4, shape=21, aes(fill=Sequence_batch)) +
    xlab(paste("PC1: 29.2%")) +
    ylab(paste("PC2: 15.3%")) +
    scale_fill_manual(name="Sequence Batch", values = c("#003566", "#ffc300", "#ed105f")) +
    plot_theme + plot_guide

save_file <- paste("clr_plot_batch.pdf", sep="")
ggsave(save_file, path = beta, scale = 1, width = 8, height = 5, units = c("in"), dpi = 300)

# label CLR PCA batch
clr_plot_batch_label <- clr_plot_batch +
    geom_label_repel(aes(label = SampleID), box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50', max.overlaps=Inf)
clr_plot_batch_label

save_file <- paste("clr_plot_batch_labelled.pdf", sep="")
ggsave(save_file, path = beta, scale = 1, width = 8, height = 5, units = c("in"), dpi = 300)

## Sediment Samples
clr_plot_final <- ggplot(ordi_data, aes(x=Axis.1clr, y = Axis2.clr)) +
    geom_hline(yintercept=0, size=.2, linetype = "dashed", color="black") + geom_vline(xintercept=0, size=.2, linetype = "dashed", color="black") +
    geom_point(color="black", size=4, aes(fill=Location_description, shape=Depth_desc)) +
    xlab(paste("PC1: 23.6%")) +
    ylab(paste("PC2: 10.8%")) +
    scale_shape_manual(name = "Depth layer", values=c(25, 21, 22)) +
    scale_fill_manual(name = "Location", values=c("#dcab6b", "#9be564")) +
    plot_theme + plot_guide
clr_plot_final

# save file
save_file <- paste("PCA_clr.pdf", sep="")
ggsave(save_file, path = beta, scale = 1, width = 9, height = 5, units = c("in"), dpi = 300)

# label plot
clr_plot_label <- clr_plot_final +
geom_label_repel(aes(label = SampleID), box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50', max.overlaps=Inf)
pca_clr_label

# save file
save_file <- paste("PCA_clr_labelled.pdf", sep="")
ggsave(save_file, path = beta, scale = 1, width = 9, height = 5, units = c("in"), dpi = 300)

### Plot PhILR PCA
## Water Samples
plot_guide <- guides(fill = guide_colourbar(frame.linewidth = 1, frame.colour = "black", ticks = TRUE, ticks.colour = "black", ticks.linewidth = 1, reverse=T),
    shape = guide_legend(order=1, override.aes = list(size = 5, color="black", alpha=1)),
    color = "none")

philr_plot_final <- ggplot(ordi_data_sub_NoBatch4, aes(x=Axis.1philr, y = Axis.2philr)) +
    geom_hline(yintercept=0, size=.2, linetype = "dashed", color="black") + geom_vline(xintercept=0, size=.2, linetype = "dashed", color="black") +
    geom_point(color="black", size=4, aes(fill=Depth_meters, shape=strat_time)) +
    xlab(paste("PC1: 53.8%")) +
    ylab(paste("PC2: 14.5%")) +
    scale_shape_manual(name = "Stratification?", values=c(21, 23, 24)) +
    scale_fill_gradient2(name = "Depth (m)", breaks = c(0, 2, 6, 10), limits = c(0, 10), low="#fffab1", mid="#ed105f", high="#023e73", midpoint=5) +
     plot_theme + plot_guide
philr_plot_final

save_file <- paste("philr.pdf", sep="")
ggsave(save_file, path = beta, scale = 1, width = 8, height = 5, units = c("in"), dpi = 300)

# label PhILR PCA plot
philr_plot_label <- philr_plot_final +
    geom_label_repel(aes(label = SampleID), box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50', max.overlaps=Inf)

save_file <- paste("philr_labelled.pdf", sep="")
ggsave(save_file, path = beta, scale = 1, width = 8, height = 5, units = c("in"), dpi = 300)

# plot batches PhILR PCA
plot_guide <- guides(shape = guide_legend(order=1, override.aes = list(size = 5, color="black", alpha=1)),
                     color = "none")

philr_plot_batch <- ggplot(ordi_data_sub_batch, aes(x=Axis.1philr, y = Axis.2philr)) +
    geom_hline(yintercept=0, size=.2, linetype = "dashed", color="black") + geom_vline(xintercept=0, size=.2, linetype = "dashed", color="black") +
    geom_point(color="black", size=4, shape=21, aes(fill=Sequence_batch)) +
    xlab(paste("PC1: 53.8%")) +
    ylab(paste("PC2: 14.5%")) +
    scale_fill_manual(name="Sequence Batch", values = c("#003566", "#ffc300", "#ed105f")) +
    plot_theme + plot_guide

save_file <- paste("philr_plot_batch.pdf", sep="")
ggsave(save_file, path = beta, scale = 1, width = 8, height = 5, units = c("in"), dpi = 300)

# label PhILR PCA batch
philr_plot_batch_label <- philr_plot_batch +
    geom_label_repel(aes(label = SampleID), box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50', max.overlaps=Inf)
philr_plot_batch_label

save_file <- paste("philr_plot_batch_labelled.pdf", sep="")
ggsave(save_file, path = beta, scale = 1, width = 8, height = 5, units = c("in"), dpi = 300)

## Sediment Samples
philr_plot_final <- ggplot(ordi_data, aes(x=Axis.1philr, y = Axis.2philr)) +
    geom_hline(yintercept=0, size=.2, linetype = "dashed", color="black") + geom_vline(xintercept=0, size=.2, linetype = "dashed", color="black") +
    geom_point(color="black", size=4, aes(fill=Location_description, shape=Depth_desc)) +
    xlab(paste("PC1: 52.7%")) +
    ylab(paste("PC2: 8.7%")) +
    scale_shape_manual(name = "Depth layer", values=c(25, 21, 22)) +
    scale_fill_manual(name = "Location", values=c("#dcab6b", "#9be564")) +
    plot_theme + plot_guide
philr_plot_final

# save file
save_file <- paste("PCA_philr.pdf", sep="")
ggsave(save_file, path = beta, scale = 1, width = 9, height = 5, units = c("in"), dpi = 300)

# label plot
philr_plot_label <- philr_plot_final +
geom_label_repel(aes(label = SampleID), box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50', max.overlaps=Inf)
pca_philr_label

# save file
save_file <- paste("PCA_philr_labelled.pdf", sep="")
ggsave(save_file, path = beta, scale = 1, width = 9, height = 5, units = c("in"), dpi = 300)

### combine CLR and PhILR PCA plots (Figure 3c, Figure 5c, Figure S6, Figure S10)
both <- plot_grid(clr_plot_final + theme(legend.position="none"),
                 philr_plot_final + theme(legend.position="none"),
                 ncol=2, align = "v", axis="b")

save_file <- paste("Combo_ordination.pdf", sep="")
ggsave(save_file, path = beta, plot = both, scale = 1, width = 10, height = 5, units = c("in"), dpi = 300)

# Water Samples Labelled
both <- plot_grid(clr_plot_batch_label + theme(legend.position="none"),
                 philr_plot_batch_label + theme(legend.position="none"),
                 ncol=3, align = "v", axis="b")

# Sediment Samples Labelled
both <- plot_grid(clr_plot_label + theme(legend.position="none"),
                 philr_plot_label + theme(legend.position="none"),
                 ncol=3, align = "v", axis="b")

save_file <- paste("Combo_ordination_batch_labelled.svg", sep="")
ggsave(save_file, path = beta, plot = both, scale = 1, width = 15, height = 5, units = c("in"), dpi = 300)

#### Stats on beta diversity ####

# obtain distance matrix CLR
ordi_clr_dist <- as.matrix(phyloseq::distance(ps_clr, method="euclidean"))
head(ordi_clr_dist)

# obtain distance matrix PhILR
ordi_philr_dist <- as.matrix(gp.dist)
head(ordi_philr_dist)

# Water ADONIS
adonis2(ordi_clr_dist~DATA_PHYLOSEQ_FIXED$Depth_cat + DATA_PHYLOSEQ_FIXED$Location_description + DATA_PHYLOSEQ_FIXED$Season, permutations=999, by="margin", method="euclidean")
adonis2(ordi_philr_dist~DATA_PHYLOSEQ_FIXED$Depth_cat + DATA_PHYLOSEQ_FIXED$Location_description + DATA_PHYLOSEQ_FIXED$Season, permutations=999, by="margin", method="euclidean")

# Sediment ADONIS
adonis2(ordi_clr_dist~DATA_PHYLOSEQ_FIXED$Location + DATA_PHYLOSEQ_FIXED$Depth_cat, permutations=999, by="margin", method="euclidean")
adonis2(ordi_philr_dist~DATA_PHYLOSEQ_FIXED$Location + DATA_PHYLOSEQ_FIXED$Depth_cat, permutations=999, by="margin", method="euclidean")

# Water Geochem ADONIS
adonis2(ordi_clr_dist~
        DATA_PHYLOSEQ_FIXED$Fe_mg_L+
        DATA_PHYLOSEQ_FIXED$ChlA_RFU+
        DATA_PHYLOSEQ_FIXED$DO_percSat+
        DATA_PHYLOSEQ_FIXED$TDS_ppt+
        DATA_PHYLOSEQ_FIXED$pH+
        DATA_PHYLOSEQ_FIXED$Temp_C+
        DATA_PHYLOSEQ_FIXED$Fe_mg_L, permutations=999, by="margin", method="euclidean")

adonis2(ordi_philr_dist~
    DATA_PHYLOSEQ_FIXED$Fe_mg_L+
    DATA_PHYLOSEQ_FIXED$ChlA_RFU+
    DATA_PHYLOSEQ_FIXED$DO_percSat+
    DATA_PHYLOSEQ_FIXED$TDS_ppt+
    DATA_PHYLOSEQ_FIXED$pH+
    DATA_PHYLOSEQ_FIXED$Temp_C+
    DATA_PHYLOSEQ_FIXED$Fe_mg_L, permutations=999, by="margin", method="euclidean")
