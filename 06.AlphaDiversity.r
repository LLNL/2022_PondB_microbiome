# 06. Alpha Diversity
# Figure 3b, Figure 5b, Figure S5, Figure S9
# Ref for phyloseq: McMurdie, P. J. & Holmes, S. phyloseq: An R Package for Reproducible Interactive Analysis and Graphics of Microbiome Census Data. PLOS ONE 8, e61217 (2013).
# Ref for picante: Kembel, S. W. et al. Picante: R tools for integrating phylogenies and ecology. Bioinformatics 26, 1463–1464 (2010).
# Water and sediment sequences were analyzed separately

#### Phyloseq Object ####
ps.noncontam.tree

#### Set Directory ####
alpha <- file.path(paste(path_phy, "/alpha_diversity", sep=""))
dir.create(alpha, showWarnings=FALSE)
setwd(alpha)

#### Get many alpha diversity incides ####
alpha.div <- microbiome::alpha(ps.noncontam.tree, index = "all")
alpha.div.phyloseq <- estimate_richness(ps.noncontam.tree, split = TRUE, measures = NULL)
head(alpha.div)
head(alpha.div.phyloseq)

write.table(alpha.div,paste("alpha.div.pkgmicrobiome.csv", sep=""), sep=",", row.names = T)
write.table(alpha.div.phyloseq,paste("alpha.div.pkgphyloseq.csv", sep=""), sep=",", row.names = T)

#### Calculate Faith's PD ####
obj_ps.pd <- pd(t(otu_table(ps.noncontam.tree)), filt_tree, include.root=T)
head(obj_ps.pd)
write.table(obj_ps.pd,paste("alpha.div.faithsPD.csv", sep=""), sep=",", row.names = T)

#### Preliminary plots ####
# Function for plotting multiple parameters (Water samples)
parameters <- c("Depth_cat", "Location_description", "Season", "Strat_group", "Sequence_batch")

# Function for plotting multiple parameters (Sediment samples)
parameters <- c("Depth_cat", "Location_description", "Depth_desc")

func_richness <- llply(as.list(parameters), function(para, x) {
    plot_alpha <- plot_richness(physeq = x, x = para) +
                    geom_boxplot(notch=FALSE) +
                    geom_jitter(shape=16, position=position_jitter(0.2))
    save_file_plot <- paste("alpha_div.barplot.", para,".pdf", sep="")
    ggsave(save_file_plot, path = alpha, scale = 1, width = 20, height = 6, units = c("in"), dpi = 300)
}, ps.noncontam.tree)

#### Statistical tests (Repeat 1-6 for every variable combinations) ####
# 1. Anova
aov.obs = aov(Observed ~ Depth_cat, data=alpha.div.metadata2)
summary.lm(aov.obs)

# 2. Extract the residuals & Run Shapiro-Wilk test
aov_residuals <- residuals(object = aov.obs)
shapiro.test(x = aov_residuals)

# 3. Check homogeneity
leveneTest(Observed ~ Loc_sec, data = alpha.div.metadata2)

# 4. Kruskal-Wallis test
kruskal.test(Observed ~ Loc_sec, data=alpha.div.metadata2)

# 5. To do all the pairwise comparisons between groups and correct for multiple comparisons
pairwise.wilcox.test(alpha.div.metadata2$Observed, alpha.div.metadata2$Loc_sec, p.adjust.method="fdr", exact=FALSE, paired=FALSE)

# 6. Plot
par(mfrow=c(2,2))
plot(aov.obs)

#### Prepare for plots Figure S5 and S9 ####

# Combine data
combo_plot <- as.data.frame(as.matrix(sample_data(ps.noncontam.tree)))
combo_plot$Observed <- alpha.div$observed[match(rownames(alpha.div), rownames(combo_plot))]
combo_plot$Shannon <- alpha.div$diversity_shannon[match(rownames(alpha.div), rownames(combo_plot))]
combo_plot$diversity_inverse_simpson <- alpha.div$diversity_inverse_simpson[match(rownames(alpha.div), rownames(combo_plot))]
combo_plot$PD <- obj_ps.pd$PD[match(rownames(obj_ps.pd), rownames(combo_plot))]

head(combo_plot)

### Melt the dataframe and fix the variable and value column
# Water samples
combo_plot$Season <- gsub("Spring", "Incipient Stratification", combo_plot$Season)
combo_plot$Season <- gsub("Summer", "Stratification (~3 months)", combo_plot$Season)
combo_plot$Season <- gsub("Fall", "Stratification (~7 months)", combo_plot$Season)
combo_plot$Strat_group <- gsub("Unstratified", "Incipient Stratification2", combo_plot$Strat_group)
combo_plot_melt <- reshape2::melt(combo_plot, id=c("SampleID", "Observed", "Shannon", "diversity_inverse_simpson", "PD"))
combo_plot_melt <- subset(combo_plot_melt, variable=="Location_description" | variable=="Strat_group" | variable=="Season")
combo_plot_melt$value = factor(combo_plot_melt$value, levels=c("Incipient Stratification", "Stratification (~3 months)", "Stratification (~7 months)", "Inlet", "Inlet-Middle", "Middle", "Middle-Outlet", "Outlet", "Incipient Stratification2", "Epilimnion", "Thermocline", "Hypolimnion"), ordered = TRUE)
combo_plot_melt$variable <- gsub("Location_description", "Location", combo_plot_melt$variable)
combo_plot_melt$variable <- gsub("Strat_group", "Stratification Layer", combo_plot_melt$variable)
combo_plot_melt$variable <- gsub("Season", "Stratification Time", combo_plot_melt$variable)
combo_plot_melt$variable <- gsub("Sequence_batch", "Sequence batch", combo_plot_melt$variable)
combo_plot_melt$variable = factor(combo_plot_melt$variable, levels=c("Location", "Stratification Time", "Stratification Layer", "Sequence batch"))
combo_plot_melt$variable2 = combo_plot_melt$variable
combo_plot_melt$variable <- NULL
head(combo_plot_melt)
write.csv(combo_plot_melt, paste0("combo_plot_melt.csv"))

# Sediment samples grouped by Location and Depth
combo_plot_select <- combo_plot[, c("SampleID", "Observed", "Shannon", "diversity_inverse_simpson", "PD", "Location_description", "Depth_desc")]
combo_plot_melt <- reshape2::melt(combo_plot_select, id=c("SampleID", "Observed", "Shannon", "diversity_inverse_simpson", "PD", "Location_description"))
combo_plot_melt$value = factor(combo_plot_melt$value, levels=c("Top 5 cm", "5-10 cm", "10-21 cm"), ordered = TRUE)
head(combo_plot_melt)
write.csv(combo_plot_melt, paste0("combo_plot_melt_LocDepth.csv"))

### Wilcoxon rank sum test
# For Sediment samples grouped by Location and Depth, change 'group_by(variable2)' to 'group_by(Location_description)'
stat.test_Obj <- combo_plot_melt %>%
  group_by(variable2) %>%
  rstatix::wilcox_test(Observed ~ value) %>%
  rstatix::adjust_pvalue(method = "fdr") %>%
  rstatix::add_significance()
stat.test_Obj
write.csv(stat.test_Obj, paste0("wilcox_test_Observed.csv"))

stat.test_PD <- combo_plot_melt %>%
  group_by(variable2) %>%
  rstatix::wilcox_test(PD ~ value) %>%
  rstatix::adjust_pvalue(method = "fdr") %>%
  rstatix::add_significance()
stat.test_PD
write.csv(stat.test_PD, paste0("wilcox_test_FaithsPD.csv"))

stat.test_Shannon <- combo_plot_melt %>%
  group_by(variable2) %>%
  rstatix::wilcox_test(Shannon ~ value) %>%
  rstatix::adjust_pvalue(method = "fdr") %>%
  rstatix::add_significance()
stat.test_Shannon
write.csv(stat.test_Shannon, paste0("wilcox_test_Shannon.csv"))

stat.test_InvSimpson <- combo_plot_melt %>%
  group_by(variable2) %>%
  rstatix::wilcox_test(diversity_inverse_simpson ~ value) %>%
  rstatix::adjust_pvalue(method = "fdr") %>%
  rstatix::add_significance()
stat.test_InvSimpson
write.csv(stat.test_InvSimpson, paste0("wilcox_test_InvSimpson.csv"))

#### Plot Figure S5 and S9 ####
# Set up theme
plot_theme <- theme(panel.background = element_rect(fill = "white", colour = "black", size = 1, linetype = "solid"),
    panel.border = element_rect(colour="black", size=1, fill=NA),
    strip.background=element_rect(fill='white', colour='white', size = 0),
    strip.text = element_text(face="bold", size=15),
    panel.spacing.x=unit(0.5, "lines"),
    panel.grid.major = element_line(size = 0),
    panel.grid.minor = element_line(size = 0),
    axis.text = element_text(size=15, colour="black"),
    axis.title = element_text(face="bold", size=15),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
    legend.position="right",
    legend.key = element_rect(fill = "white"),
    legend.title = element_text(face="bold", size=15),
    legend.text = element_text(size=15))

# Plot all the indices

### Observed
# Water samples
plot_nomargins_Obs <- scale_y_continuous(expand = expansion(mult = c(0, 0)), limits=c(200,3200), breaks=c(500, 1000, 1500, 2000, 2500, 3000))
stat.test_Obj <- stat.test_Obj %>% rstatix::add_y_position(fun = "max", scales = "free", step.increase=0.13, comparisons = list(c("Incipient Stratification2", "Epilimnion"), c("Incipient Stratification2", "Hypolimnion"),c("Epilimnion", "Thermocline"), c("Epilimnion", "Hypolimnion"), c("Thermocline", "Hypolimnion")))

# Sediment samples (pound out 'stat_pvalue_manual' for grouping by Sediment and Depth)
plot_nomargins_Obs <- scale_y_continuous(expand = expansion(mult = c(0, 0)), limits=c(0,4200), breaks=c(0, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000))

# Plot Observed
alpha_plot_Obs <- ggplot(combo_plot_melt, aes(x=value, y=Observed)) +
    stat_boxplot(geom = "errorbar", width = 0.2) +
    geom_boxplot(outlier.size=0, outlier.colour="white", fill="gray88") +
    geom_jitter(size=2, position=position_jitter(0.2), alpha=0.8, shape=21, fill="grey88") +
    ylab("Observed") +
    facet_grid(. ~ variable2, scales="free") +
    stat_pvalue_manual(stat.test_Obj, hide.ns=TRUE, label="p.adj.signif", tip.length=0.01, label.size = 7, bracket.size=0.7, step.group.by="variable2") +
    plot_theme + plot_nomargins_Obs
alpha_plot_Obs

### PD
# Water samples
plot_nomargins_PD <- scale_y_continuous(expand = expansion(mult = c(0, 0)), limits=c(40,250), breaks=c(50, 100, 150, 200, 250))
stat.test_PD <- stat.test_PD %>% rstatix::add_y_position(fun = "max", scales = "free", step.increase=0.13, comparisons = list(c("Incipient Stratification", "Stratification (~3 months)"), c("Incipient Stratification", "Stratification (~7 months)"), c("Incipient Stratification2", "Thermocline"), c("Incipient Stratification2", "Hypolimnion"),c("Epilimnion", "Thermocline"), c("Epilimnion", "Hypolimnion"), c("Thermocline", "Hypolimnion")))

# Sediment samples (pound out 'stat_pvalue_manual' for grouping by Sediment and Depth)
plot_nomargins_PD <- scale_y_continuous(expand = expansion(mult = c(0, 0)), limits=c(0,310), breaks=c(0, 50, 100, 150, 200, 250, 300))

# Plot PD
alpha_plot_PD <- ggplot(combo_plot_melt, aes(x=value, y=PD)) +
    stat_boxplot(geom = "errorbar", width = 0.2) +
    geom_boxplot(outlier.size=0, outlier.colour="white", fill="gray88") +
    geom_jitter(size=2, position=position_jitter(0.2), alpha=0.8, shape=21, fill="grey88") +
    ylab("Faith's PD") +
    facet_grid(. ~ variable2, scales="free") +
    stat_pvalue_manual(stat.test_PD, hide.ns=TRUE, label="p.adj.signif", tip.length=0.01, label.size = 7, bracket.size=0.7, step.group.by="variable2") +
    plot_theme + plot_nomargins_PD
alpha_plot_PD

### Shannon
# Water samples
plot_nomargins_Shannon <- scale_y_continuous(expand = expansion(mult = c(0, 0)), limits=c(3,8.8), breaks=c(4, 5, 6, 7, 8))
stat.test_Shannon <- stat.test_Shannon %>% rstatix::add_y_position(fun = "max", scales = "free", step.increase=0.13, comparisons = list(c("Incipient Stratification", "Stratification (~3 months)"), c("Incipient Stratification", "Stratification (~7 months)"), c("Inlet", "Middle"), c("Incipient Stratification2", "Thermocline"), c("Incipient Stratification2", "Hypolimnion"), c("Epilimnion", "Thermocline"), c("Epilimnion", "Hypolimnion"), c("Thermocline", "Hypolimnion")))

# Sediment samples (pound out 'stat_pvalue_manual' for grouping by Sediment and Depth)
plot_nomargins_Shannon <- scale_y_continuous(expand = expansion(mult = c(0, 0)), limits=c(3.7,7.5), breaks=c(4, 4.5, 5, 5.5, 6, 6.5, 7))

# Plot Shannon
alpha_plot_Shannon <- ggplot(combo_plot_melt, aes(x=value, y=Shannon)) +
    stat_boxplot(geom = "errorbar", width = 0.2) +
    geom_boxplot(outlier.size=0, outlier.colour="white", fill="gray88") +
    geom_jitter(size=2, position=position_jitter(0.2), alpha=0.8, shape=21, fill="grey88") +
    ylab("Shannon") +
    facet_grid(. ~ variable2, scales="free") +
    stat_pvalue_manual(stat.test_Shannon, hide.ns=TRUE, label="p.adj.signif", tip.length=0.01, label.size = 7, bracket.size=0.7, step.group.by="variable2") +
    plot_theme + plot_nomargins_Shannon
alpha_plot_Shannon

### Inverse Simpson
# Water Samples
plot_nomargins_InvSimpson <- scale_y_continuous(expand = expansion(mult = c(0, 0)), limits=c(0,325), breaks=c(0, 100, 200, 300))
stat.test_InvSimpson <- stat.test_InvSimpson %>% rstatix::add_y_position(fun = "max", scales = "free", step.increase=0.13, comparisons = list(c("Incipient Stratification2", "Epilimnion"), c("Incipient Stratification2", "Hypolimnion"), c("Epilimnion", "Thermocline"), c("Epilimnion", "Hypolimnion"), c("Thermocline", "Hypolimnion")))

# Sediment samples (pound out 'stat_pvalue_manual' for grouping by Sediment and Depth)
plot_nomargins_InvSimpson <- scale_y_continuous(expand = expansion(mult = c(0, 0)), limits=c(0,475), breaks=c(0, 100, 200, 300, 400))

# Plot Inverse Simpson
alpha_plot_InvSimpson <- ggplot(combo_plot_melt, aes(x=value, y=diversity_inverse_simpson)) +
    stat_boxplot(geom = "errorbar", width = 0.2) +
    geom_boxplot(outlier.size=0, outlier.colour="white", fill="gray88") +
    geom_jitter(size=2, position=position_jitter(0.2), alpha=0.8, shape=21, fill="grey88") +
    ylab("Inverse Simpson") +
    facet_grid(. ~ variable2, scales="free") +
    stat_pvalue_manual(stat.test_InvSimpson, hide.ns=TRUE, label="p.adj.signif", tip.length=0.01, label.size = 7, bracket.size=0.7, step.group.by="variable2") +
    plot_theme + plot_nomargins_InvSimpson
alpha_plot_InvSimpson

# Combine
both <- plot_grid(alpha_plot_Obs + theme(axis.text.x = element_blank(), legend.position="none"),
                  alpha_plot_Shannon + theme(axis.text.x = element_blank(), legend.position="none", strip.text = element_blank()),
                  alpha_plot_InvSimpson + theme(axis.text.x = element_blank(), legend.position="none", strip.text = element_blank()),
                  alpha_plot_PD + theme(legend.position="none", strip.text = element_blank()),
                  ncol=1, align = "v", axis="b", rel_heights = c(0.7,0.7,0.7,1))
both
save_file <- paste("alpha.div.shannon.simpson.pd.obs.variables.svg", sep="")
ggsave(save_file, path = alpha, plot = both, scale = 1, width = 10, height = 16, units = c("in"), dpi = 300)

#### Prepare for plotting only Stratification Layer (Figure 3b, Figure 5b)  ####
# Water samples
combo_plot_melt_sub <- subset(combo_plot_melt, variable2=="Stratification Layer")
combo_plot_melt_sub$value <- gsub('Incipient Stratification2', 'Incipient Stratification', combo_plot_melt_sub$value)
combo_plot_melt_sub$value = factor(combo_plot_melt_sub$value, levels=c("Incipient Stratification", "Epilimnion", "Thermocline", "Hypolimnion"), ordered = TRUE)

head(combo_plot_melt_sub)

# Sediment samples grouped by Location
combo_plot_melt_sub <- reshape2::melt(combo_plot_select, id=c("SampleID", "Observed", "Shannon", "diversity_inverse_simpson", "PD"))
combo_plot_melt_sub <- subset(combo_plot_melt_sub, variable=="Location_description")
combo_plot_melt_sub$value = factor(combo_plot_melt_sub$value, levels=c("Inlet", "Outlet"), ordered = TRUE)
head(combo_plot_melt_sub)
write.csv(combo_plot_melt, paste0("combo_plot_melt_Loc.csv"))

### Wilcoxon rank sum test
# For Sediment samples grouped by Location, pound out 'group_by(variable2)' and change 'value' to 'Location_description'
stat.test_Obj <- combo_plot_melt_sub %>%
  group_by(variable2) %>%
  rstatix::wilcox_test(Observed ~ value) %>%
  rstatix::adjust_pvalue(method = "fdr") %>%
  rstatix::add_significance()
stat.test_Obj
write.csv(stat.test_Obj, paste0("wilcox_test_Observed.csv"))

stat.test_PD <- combo_plot_melt_sub %>%
  group_by(variable2) %>%
  rstatix::wilcox_test(PD ~ value) %>%
  rstatix::adjust_pvalue(method = "fdr") %>%
  rstatix::add_significance()
stat.test_PD
write.csv(stat.test_PD, paste0("wilcox_test_FaithsPD.csv"))

stat.test_Shannon <- combo_plot_melt_sub %>%
  group_by(variable2) %>%
  rstatix::wilcox_test(Shannon ~ value) %>%
  rstatix::adjust_pvalue(method = "fdr") %>%
  rstatix::add_significance()
stat.test_Shannon
write.csv(stat.test_Shannon, paste0("wilcox_test_Shannon.csv"))

stat.test_InvSimpson <- combo_plot_melt_sub %>%
  group_by(variable2) %>%
  rstatix::wilcox_test(diversity_inverse_simpson ~ value) %>%
  rstatix::adjust_pvalue(method = "fdr") %>%
  rstatix::add_significance()
stat.test_InvSimpson
write.csv(stat.test_InvSimpson, paste0("wilcox_test_InvSimpson.csv"))

#### Plot Figure 3b, Figure 5b ####

# Set up theme
plot_theme <- theme(panel.background = element_rect(fill = "white", colour = "black", size = 1, linetype = "solid"),
    panel.border = element_rect(colour="black", size=1, fill=NA),
    strip.background=element_rect(fill='white', colour='white', size = 0),
    strip.text = element_text(face="bold", size=15),
    panel.spacing.x=unit(0.5, "lines"),
    panel.grid.major = element_line(size = 0),
    panel.grid.minor = element_line(size = 0),
    axis.text = element_text(size=15, colour="black"),
    axis.title = element_text(face="bold", size=15),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
    legend.position="right",
    legend.key = element_rect(fill = "white"),
    legend.title = element_text(face="bold", size=15),
    legend.text = element_text(size=15),
    plot.margin=unit(c(t = 0.5, r = 0.2, b = 0, l = 0.2),"cm"))

### Observed
# Water samples
plot_nomargins_Obs <- scale_y_continuous(expand = expansion(mult = c(0, 0)), limits=c(200,3200), breaks=c(200, 1000, 2000, 3000))
stat.test_Obj <- stat.test_Obj %>% rstatix::add_y_position(fun = "max", scales = "free", step.increase=0.13)

# Sediment samples
plot_nomargins_Obs <- scale_y_continuous(expand = expansion(mult = c(0, 0)), limits=c(0,4200), breaks=c(0, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000))
stat.test_Obj <- stat.test_Obj %>% rstatix::add_y_position(fun = "max", scales = "free")

# Plot Observed
alpha_plot_Obs <- ggplot(combo_plot_melt_sub, aes(x=value, y=Observed)) +
    stat_boxplot(geom = "errorbar", width = 0.2) +
    geom_boxplot(outlier.size=0, outlier.colour="white", fill="gray88") +
    geom_jitter(size=2, position=position_jitter(0.2), alpha=0.8, shape=21, fill="grey88") +
    ylab("Observed") +
    stat_pvalue_manual(stat.test_Obj, hide.ns=TRUE, label="p.adj.signif", tip.length=0.01, label.size = 7, bracket.size=0.7) +
    plot_theme + plot_nomargins_Obs
alpha_plot_Obs

### PD
# Water samples
plot_nomargins_PD <- scale_y_continuous(expand = expansion(mult = c(0, 0)), limits=c(40,250), breaks=c(50, 100, 150, 200, 250))
stat.test_PD <- stat.test_PD %>% rstatix::add_y_position(fun = "max", scales = "free", step.increase=0.13)

# Sediment samples
plot_nomargins_PD <- scale_y_continuous(expand = expansion(mult = c(0, 0)), limits=c(0,310), breaks=c(0, 50, 100, 150, 200, 250, 300))
stat.test_PD <- stat.test_PD %>% rstatix::add_y_position(fun = "max", scales = "free")

# Plot PD
alpha_plot_PD <- ggplot(combo_plot_melt_sub, aes(x=value, y=PD)) +
    stat_boxplot(geom = "errorbar", width = 0.2) +
    geom_boxplot(outlier.size=0, outlier.colour="white", fill="gray88") +
    geom_jitter(size=2, position=position_jitter(0.2), alpha=0.8, shape=21, fill="grey88") +
    ylab("Faith's PD") +
    stat_pvalue_manual(stat.test_PD, hide.ns=TRUE, label="p.adj.signif", tip.length=0.01, label.size = 7, bracket.size=0.7) +
    plot_theme + plot_nomargins_PD
alpha_plot_PD

### Shannon
# Water samples
plot_nomargins_Shannon <- scale_y_continuous(expand = expansion(mult = c(0, 0)), limits=c(3,8.8), breaks=c(4, 5, 6, 7, 8))
stat.test_Shannon <- stat.test_Shannon %>% rstatix::add_y_position(fun = "max", scales = "free", step.increase=0.13)

# Sediment samples
plot_nomargins_Shannon <- scale_y_continuous(expand = expansion(mult = c(0, 0)), limits=c(3.7,7.5), breaks=c(4, 4.5, 5, 5.5, 6, 6.5, 7))
stat.test_Shannon <- stat.test_Shannon %>% rstatix::add_y_position(fun = "max", scales = "free")

# Plot Shannon
alpha_plot_Shannon <- ggplot(combo_plot_melt_sub, aes(x=value, y=Shannon)) +
    stat_boxplot(geom = "errorbar", width = 0.2) +
    geom_boxplot(outlier.size=0, outlier.colour="white", fill="gray88") +
    geom_jitter(size=2, position=position_jitter(0.2), alpha=0.8, shape=21, fill="grey88") +
    ylab("Shannon") +
    stat_pvalue_manual(stat.test_Shannon, hide.ns=TRUE, label="p.adj.signif", tip.length=0.01, label.size = 7, bracket.size=0.7) +
    plot_theme + plot_nomargins_Shannon
alpha_plot_Shannon

### Inverse Simpson
# Water Samples
plot_nomargins_InvSimpson <- scale_y_continuous(expand = expansion(mult = c(0, 0)), limits=c(0,350), breaks=c(0, 100, 200, 300))
stat.test_InvSimpson <- stat.test_InvSimpson %>% rstatix::add_y_position(fun = "max", scales = "free", step.increase=0.13)

# Sediment Samples
plot_nomargins_InvSimpson <- scale_y_continuous(expand = expansion(mult = c(0, 0)), limits=c(0,475), breaks=c(0, 100, 200, 300, 400))
stat.test_InvSimpson <- stat.test_InvSimpson %>% rstatix::add_y_position(fun = "max", scales = "free")

alpha_plot_InvSimpson <- ggplot(combo_plot_melt_sub, aes(x=value, y=diversity_inverse_simpson)) +
    stat_boxplot(geom = "errorbar", width = 0.2) +
    geom_boxplot(outlier.size=0, outlier.colour="white", fill="gray88") +
    geom_jitter(size=2, position=position_jitter(0.2), alpha=0.8, shape=21, fill="grey88") +
    ylab("Inverse Simpson") +
    stat_pvalue_manual(stat.test_InvSimpson, hide.ns=TRUE, label="p.adj.signif", tip.length=0.01, label.size = 7, bracket.size=0.7) +
    plot_theme + plot_nomargins_InvSimpson
alpha_plot_InvSimpson

### Combine plots
both <- plot_grid(alpha_plot_Obs,
                  alpha_plot_Shannon,
                  alpha_plot_InvSimpson,
                  alpha_plot_PD,
                  ncol=4, align = "v", axis="b", rel_widths=c(0.8,0.8,0.8,0.8)) + theme(plot.margin=unit(c(t = 0, r = 0.5, b = 0, l = 1.5),"cm"))
both

save_file <- paste("alpha.div.shannon.simpson.pd.obs.variables.svg", sep="")
ggsave(save_file, path = alpha, plot = both, scale = 1, width = 10, height = 5, units = c("in"), dpi = 300)
