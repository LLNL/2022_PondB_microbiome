# 05. Barplots
# Figure 3a, Figure 5a, Figure S4
# Ref for phyloseq: McMurdie, P. J. & Holmes, S. phyloseq: An R Package for Reproducible Interactive Analysis and Graphics of Microbiome Census Data. PLOS ONE 8, e61217 (2013).
# Water and sediment sequences were analyzed separately

#### Prepare phyloseq object ####
filt_tree = read.tree("<path to tree.newick>")
filt_tree

DATA_PHYLOSEQ <- read_excel("<path to metadata>", na = "NA")
DATA_PHYLOSEQ_FIXED <- data.frame(DATA_PHYLOSEQ, row.names = 1)
DATA_PHYLOSEQ_FIXED$Sample_abbrev <- row.names(DATA_PHYLOSEQ_FIXED)
DATA_PHYLOSEQ_FIXED

ps.noncontam.tree <- phyloseq(
    otu_table(otu_table(ps_filt5), taxa_are_rows = TRUE),
    sample_data(DATA_PHYLOSEQ_FIXED),
    tax_table(tax_table(ps_filt5)),
    phy_tree(filt_tree))
ps.noncontam.tree

#### Set Directory ####
barplots <- file.path(paste(path_phy, "/barplots", sep=""))
dir.create(barplots, showWarnings = FALSE)

#### Get percentages of each taxa ####
glom <- tax_glom(ps.noncontam.tree, taxrank="Family") # Change taxrank
trans <- transform_sample_counts(glom, function(OTU) OTU/sum(OTU) )
OTU1 = as(otu_table(trans), "matrix")
if(taxa_are_rows(ps.noncontam.tree)){OTU1 <- t(OTU1)}

# Coerce to data.frame
OTUdf = as.data.frame(t(OTU1))
write.csv(OTUdf, "Family_percentages.csv")

#### Plot barplots ####
# Set colors
colors <- read_excel("<path to excel file with one column for taxa and another column for hex color>")
head(colors)

# Set up theme
plot_theme <- theme(panel.background = element_rect(fill = "white", colour = "black", size = 1, linetype = "solid"),
    panel.border = element_rect(colour="black", size=1, fill=NA),
    strip.background=element_rect(fill='#d9d9d9', size = 1, colour="black"),
    strip.text = element_text(face="bold", size=30),
    panel.spacing.x=unit(2, "lines"),
    panel.grid.major = element_line(size = 0),
    panel.grid.minor = element_line(size = 0),
    axis.text = element_text(size=30, colour="black"),
    axis.title = element_text(face="bold", size=30),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
    legend.position="bottom",
    legend.key = element_rect(fill = "white"),
    legend.title = element_text(face="bold", size=30),
    legend.text = element_text(size=30),
    plot.margin=unit(c(3,1,1,2),"cm"))
plot_guides <- guides(colour="none", fill = guide_legend(ncol=5))
plot_nomargins_y <- scale_y_continuous(expand = expansion(mult = c(0, 0)), labels = function(x) paste0(x*100, "%"))
plot_nomargins_x <- scale_x_discrete(expand = expansion(mult = c(0, 0)))

# Get each taxonomic level by relative abundance
    glom <- tax_glom(ps.noncontam.tree, taxrank="Class")
    trans <- transform_sample_counts(glom, function(OTU) OTU/sum(OTU) )
    dat <- psmelt(trans) # create dataframe from phyloseq object
    dat$Class <- as.character(dat$Class) # convert taxrank to a character vector from a factor
    dat <- dat[!(dat$Abundance == 0),]
    dat$Class[dat$Abundance <= 0.01] <- 'Less than 1%' # find taxrank whose rel. abund. is less than 1%

# group the data (Water samples)
    columns = c("Depth_cat", "Location", "Season", "Class")
    dat_grouped <- dat %>%
        group_by_at(vars(one_of(columns))) %>%
        summarise(grouped_Abundance = sum(Abundance), avg = mean(Abundance), sd = sd(Abundance), n=n())

# group the data by taxrank, Sample, and Location (Sediment samples)
    columns = c("Depth_cat", "Location_description", "Depth_desc")
    dat_grouped <- dat %>%
        group_by_at(vars(one_of(level, columns))) %>%
        summarise(grouped_Abundance = sum(Abundance), avg = mean(Abundance), sd = sd(Abundance), n=n())

# Take average and stdev of spring samples (Water samples)
    dat_grouped_sub <- dat[ which(dat$Season=='Spring'), ] # substitute the Spring values
    dat_grouped_sub$Depth_cat <- 'Incipient Stratification Avg'
    columns = c("Depth_cat", "Location", "Season", "Class")
    dat_grouped_sub <- dat_grouped_sub %>%
        group_by_at(vars(one_of(columns))) %>%
        summarise(grouped_Abundance = sum(Abundance), avg=mean(Abundance), sd=sd(Abundance), n=n())
    dat_grouped_sub <- dat_grouped_sub[ which(dat_grouped_sub$Location!='Middle-Outlet'), ]
    dat_grouped_summer <- dat_grouped[ which(dat_grouped$Season=='Summer'), ]
    dat_grouped_fall <- dat_grouped[ which(dat_grouped$Season=='Fall'), ]
    dat_grouped_fall$Location <- gsub('Middle', 'Middle (Stratified ~7 months)', dat_grouped_fall$Location)
    dat_grouped_summer$Location <- gsub('Middle', 'Middle (Stratified ~3 months)', dat_grouped_summer$Location)
    dat_grouped_sub$Location <- gsub('Middle', 'Middle (Stratified ~3 months)', dat_grouped_sub$Location) # match labeling of summer samples
    dat_grouped2 <- rbind(dat_grouped_summer, dat_grouped_sub, dat_grouped_fall)
                                     
# Order the groups (Water samples)
    dat_grouped2$Condition = factor(dat_grouped2$Location, levels=c("Inlet","Inlet-Middle","Middle", "Middle-Outlet", "Outlet"), ordered = TRUE)
    dat_grouped2$Depth_cat = factor(dat_grouped2$Depth_cat, levels=c("Incipient Stratification Avg","0-1m","0.75m","1-2m","2-3m","3-4m","4-5m","5-6m","6-7m","7-8m","8-9m","9-10m"), ordered = TRUE)
    dat_grouped2 = dat_grouped2[order(dat_grouped2$Location, dat_grouped2$Depth_cat), ] # order the data frame as desired

# Order the groups (Sediment samples)
    dat_grouped$Depth_cat = factor(dat_grouped$Depth_cat, levels=c("0-1cm","1-2cm","2-3cm","3-4cm","4-5cm","6-7cm", "8-9cm", "10-11cm", "15-16cm", "20-21cm"), ordered = TRUE)
    dat_grouped$Location_description = factor(dat_grouped$Location_description, levels=c("Inlet", "Outlet"), ordered = TRUE)
    dat_grouped$Depth_desc = factor(dat_grouped$Depth_desc, levels=c("Top 5 cm", "5-10 cm", "10-21 cm"), ordered = TRUE)

# manipulation of the taxrank list
    dat_grouped2$Class <- str_replace(dat_grouped2$Class, "NA", "Unclassified") #rename "NA" to "Unclassified"
    dat_grouped2$Class <- str_replace(dat_grouped2$Class, "uncultured", "Unclassified") #rename "uncultured" to "Unclassified"
    dat_grouped2$Class <- fct_relevel(dat_grouped2$Class, "Unclassified", after = Inf) #make "Unclassified" last
    dat_grouped2$Class <- fct_relevel(dat_grouped2$Class, "Less than 1%", after = Inf) #make "Less than 1%" last

# Color code
    col <- as.character(colors$Color)
    names(col) <- as.character(colors$Taxa)

# Plot
                                     
    plot <- ggplot(data=dat_grouped2, aes(x=Depth_cat, y=grouped_Abundance, fill=Class)) +
        scale_fill_manual(values=col) +
        facet_grid(. ~ Location, scales="free", space="free") +
        labs(y = "Relative Abundance\n", x = "Depth (m)", fill = "Class") +
        geom_bar(aes(), stat="identity", position="fill", color="black", width=0.9) +
        plot_nomargins_y + plot_nomargins_x + plot_theme + plot_guides

    # save plot
        save_file_plot <- paste("barplot.basic.Class.great1perc.split.avg_colorsMod_2.pdf", sep="") #change the file name if need to
        ggsave(save_file_plot, path = barplots, scale = 1, width = 40, height = 15, units = c("in"), dpi = 300)
        write.csv(dat_grouped2, paste0("grouped_rel_abundances_great1perc.Class.split_avg.csv"))

#### Function to plot all levels >1% (Here, used for plotting Incipient Stratification Figure S4) ####
# plot each taxonomic level by relative abundance
taxrank <- c("Kingdom","Phylum","Class","Order")
func_plotbar <- llply(as.list(taxrank), function(level, x) {
        glom <- tax_glom(x, taxrank=level)
        trans <- transform_sample_counts(glom, function(OTU) OTU/sum(OTU) )
        dat <- psmelt(trans) # create dataframe from phyloseq object
        dat[[level]] <- as.character(dat[[level]]) # convert taxrank to a character vector from a factor
        dat <- dat[!(dat$Abundance == 0),]
        dat[[level]][dat$Abundance <= 0.01] <- 'Less than 1%' # find taxrank whose rel. abund. is less than 1%
    # group the data by taxrank, Sample, and Location
        columns = c("Depth_cat", "Location", "Season")
        dat_grouped <- dat %>%
            group_by_at(vars(one_of(level, columns))) %>%
            summarise(grouped_Abundance = sum(Abundance), avg=mean(Abundance), sd=sd(Abundance), n=n())
    # manipulation of the taxrank list
        dat_grouped[[level]] <- str_replace(dat_grouped[[level]], "NA", "Unclassified") #rename "NA" to "Unclassified"
        dat_grouped[[level]] <- fct_relevel(dat_grouped[[level]], "Unclassified", after = Inf) #make "Unclassified" last
        dat_grouped[[level]] <- fct_relevel(dat_grouped[[level]], "Less than 1%", after = Inf) #make "Less than 1%" last
    # Extract only spring samples
        dat_grouped_sub <- dat_grouped[ which(dat_grouped$Season=='Spring'), ]
    # Set color palette to accommodate the number of taxa
        colourCount = length(unique(dat_grouped_sub[[level]]))
        getPalette = colorRampPalette(brewer.pal(8, "Accent"))
        print(colourCount)
    # Order the groups
        dat_grouped_sub$Condition = factor(dat_grouped_sub$Location, levels=c("Inlet","Inlet-Middle","Middle", "Middle-Outlet", "Outlet"), ordered = TRUE)
        dat_grouped_sub$Depth_cat = factor(dat_grouped_sub$Depth_cat, levels=c("0-1m","1-2m","2-3m","3-4m","4-5m","5-6m","6-7m","7-8m"), ordered = TRUE)
        dat_grouped_sub = dat_grouped_sub[order(dat_grouped_sub$Location, dat_grouped_sub$Depth_cat), ] # order the data frame as desired
    # plot
        plot <- ggplot(data=dat_grouped_sub, aes_string(x="Depth_cat", y="avg", fill=level)) +
            scale_fill_manual(values=getPalette(colourCount)) +
            facet_grid(. ~ Location, scales="free", space="free") +
            labs(y = "Relative Abundance\n\n", x = "Depth (m)", fill = level) +
            geom_bar(aes(), stat="identity", position="fill", color="black", width=0.9) +
            plot_nomargins_y + plot_nomargins_x + plot_theme + plot_guides
    print(plot)
    # save plot
        save_file_plot <- paste("unstratified.barplot.basic.", level,"great1perc_avg.pdf", sep="") #change the file name if need to
        ggsave(save_file_plot, path = barplots, scale = 1, width = 30, height = 11, units = c("in"), dpi = 300)
        write.csv(dat_grouped_sub, paste0("unstratified.grouped_rel_abundances_great1perc.",level,"_avg.csv"))
}, ps.noncontam.tree)
