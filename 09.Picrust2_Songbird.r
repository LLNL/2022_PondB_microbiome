#  09.Picrust2_Songbird.r
#
# Figure 6, Table S13, Table S14, Table S15
# Ref for PICRUSt2: Douglas, G. M. et al. PICRUSt2 for prediction of metagenome functions. Nat. Biotechnol. 38, 685–688 (2020).
# Ref for Songbird: Morton, J. T. et al. Establishing microbial composition measurement standards with reference frames. Nat. Commun. 10, 2719 (2019).
# Ref for Qurro: Fedarko, M. W. et al. Visualizing ’omic feature rankings and log-ratios using Qurro. NAR Genomics Bioinforma. 2, (2020).
# Water and sediment sequences were analyzed separately

#### Phyloseq Object ####
ps.noncontam.tree

#### Set Directory ####
picrust <- file.path(paste(path_phy, "PICRUSt2", sep=""))
dir.create(picrust, showWarnings=FALSE)
setwd(picrust)

#### Prepare data for PICRUSt2 ####

# make a biomformat otu table
otu <- as(otu_table(ps.noncontam.tree),"matrix")
otu_biom <- make_biom(data=otu)
write_biom(otu_biom,"otu_biom.biom")
write.table(otu_table(ps.noncontam.tree), file = "otu_table.txt", sep = "\t", row.names = TRUE, col.names = NA)

# add sequences from "tax_table_for_tree.txt" file in folder "tree"

# export metadata table
write.table(sample_data(ps.noncontam.tree), file = "metadata_for_picrust2.txt", sep = "\t", row.names = TRUE, col.names = NA)

#### In Bash ####
wd=<path to PICRUSt2 working dir>
cd $wd
input=${wd}/otu_table.txt
output=${wd}/sequences_for_picrust2.fa

sed 's/"//g' metadata_for_picrust2.txt > metadata_for_picrust2_fixed.txt
# also add #SampleID to header

conda activate picrust2

# change these variables
study_seq=sequences_for_picrust2.fa
abundance_table=otu_biom.biom
date=<date>
rep=<repeat number>

# do not change these variables
output=${date}_picrust2_out_rep${rep}
db_path_wd=<path to PICRUSt2 default files>
ref_db=${db_path_wd}/prokaryotic/pro_ref
traits_table=${db_path_wd}/prokaryotic/ec-modified.txt.gz
copy_number_16S=${db_path_wd}/prokaryotic/16S.txt.gz
map_file=${db_path_wd}/pathway_mapfiles/metacyc_path2rxn_struc_filt_pro-modified.txt
regroup_map_file=${db_path_wd}/pathway_mapfiles/ec_level4_to_metacyc_rxn-modified.tsv

cd $wd

# run commands
picrust2_pipeline.py -s ${study_seq} -i ${abundance_table} -o ${output} -p 1 \
    --ref_dir ${ref_db} \
    --custom_trait_tables ${traits_table} \
    --marker_gene_table ${copy_number_16S} \
    --pathway_map ${map_file} \
    --reaction_func ${traits_table} \
    --regroup_map ${regroup_map_file} \
    --verbose \
    --stratified --per_sequence_contrib --coverage

# after done
add_descriptions.py -i ${output}/pathways_out/path_abun_unstrat.tsv.gz \
                    --custom_map_table ${db_path_wd}/description_mapfiles/metacyc_pathways_info-modified.txt.gz \
                    -o ${output}/pathways_out/path_abun_unstrat_descrip.tsv.gz

add_descriptions.py -i ${output}/ec-modified.txt_metagenome_out/pred_metagenome_unstrat.tsv.gz \
                    --custom_map_table ${db_path_wd}/description_mapfiles/ec_level4_info-modified.tsv.gz \
                    -o ${output}/ec-modified.txt_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz

#### the file OTU contribution for EC file is large; use bash to obtain lines of interest ####

# first create a file listing the EC of interest

EC_list=pathways_out_per_seq_contrib/path_list_all_selected.tsv
data=pathways_out_per_seq_contrib/path_abun_contrib.tsv

# move the pathway column (second column) to the first column
awk 'BEGIN {OFS="\t"} { print $2,$1,$3,$4,$5,$6,$7,$8,$9 }' ${data} > ${data}2

# search for pathways from EC_list file and print out matches from the datafile
awk 'FNR==NR {a[$1]; next}; $1 in a' ${EC_list} ${data}2 2>&1 | tee -a pathways_out_per_seq_contrib/path_abun_contrib_all_selected.tsv

#### Back to R: Prepare data ####

# Import the predicted Metacyc pathways
pred_pathway <- read.delim(<"path to file path_abun_unstrat_descrip.tsv">), row.names=1, sep = "\t", stringsAsFactors=FALSE, fileEncoding="latin1")
use_table <- pred_pathway

# Create "otu table" and "tax table"
tax_table_picrust2 <- use_table[names(use_table) %in% c("description") ]
tax_table_picrust2$ID <- rownames(tax_table_picrust2)
otu_table_picrust2 <- use_table[,-1]
otu_table_picrust2 <- t(otu_table_picrust2)

# Convert to matrix
tax_table_picrust2 <- as.matrix(tax_table_picrust2)
otu_table_picrust2 <- as.matrix(otu_table_picrust2)
rownames(tax_table_picrust2) <- gsub(':', '.', rownames(tax_table_picrust2))
colnames(otu_table_picrust2) <- gsub(':', '.', colnames(otu_table_picrust2))

# Retrieve metadata
DATA_PHYLOSEQ_FIXED_picrust2 <- DATA_PHYLOSEQ_FIXED
rownames(DATA_PHYLOSEQ_FIXED_picrust2) <- DATA_PHYLOSEQ_FIXED_picrust2$SampleID

# Make sure it looks good
head(tax_table_picrust2)
head(otu_table_picrust2)

# Make Phyloseq object
ps_picrust <- phyloseq(
    otu_table(otu_table_picrust2, taxa_are_rows = FALSE),
    sample_data(DATA_PHYLOSEQ_FIXED_picrust2),
    tax_table(tax_table_picrust2)
)

ps_picrust

#### Songbird for picrust2: prepare data ####
songbird_picrust <- file.path(paste(picrust, "/Songbird", sep=""))
dir.create(songbird_picrust, showWarnings=FALSE)
setwd(songbird_picrust)

# format taxonomy table
tax <- as(tax_table(ps_picrust),"matrix")
tax <- as.data.frame(tax)
tax$picrust <- paste("p", rownames(tax), sep="__")
tax$description <- paste("d", tax$description, sep="__")
tax_cols <- c("picrust", "description")
tax$concat <- do.call(paste, c(tax[tax_cols], sep=";"))
for(co in tax_cols) tax[co]<-NULL
head(tax)
write.table(tax, "tax_for_qiime2.txt", quote=FALSE, col.names=FALSE, sep="\t")

# make a biomformat otu table
otu <- as(otu_table(ps_picrust),"matrix")
otu <- t(otu)
otu_biom <- make_biom(data=otu)
write_biom(otu_biom,"otu_biom.biom")
write.table(otu_table(ps_picrust), file = "picrust_otu_table.txt", sep = "\t", row.names = TRUE, col.names = NA)

# export metadata table
write.table(sample_data(ps_picrust), file = "metadata_for_qiime2.txt", sep = "\t", row.names = TRUE, col.names = NA)

#### Songbird for picrust2: to Bash ####
conda activate qiime2-2020.6
wd=<path to PICRUSt2 folder>
cd $wd

# prepare data for QIIME2
sed 's/"//g' metadata_for_qiime2.txt > metadata_for_qiime2_fixed.txt
# also add #SampleID to header

biom convert -i otu_biom.biom -o otu_biom_HDF5.biom --to-hdf5

biom add-metadata -i otu_biom_HDF5.biom -o otu_wTax_metadata.biom --observation-metadata-fp tax_for_qiime2.txt --sc-separated taxonomy --observation-header OTUID,taxonomy --sample-metadata-fp metadata_for_qiime2_fixed.txt

# Import to QIIME2
qiime tools import \
    --input-path otu_biom_HDF5.biom \
    --type 'FeatureTable[Frequency]' \
    --input-format BIOMV210Format \
    --output-path feature-table.qza

qiime tools import \
    --type 'FeatureData[Taxonomy]' \
    --input-format HeaderlessTSVTaxonomyFormat \
    --input-path tax_for_qiime2.txt \
    --output-path taxonomy.qza

# Check import
qiime feature-table summarize \
    --i-table feature-table.qza \
    --m-sample-metadata-file metadata_for_qiime2_fixed.txt \
    --o-visualization summary_vis.qzv

qiime tools view summary_vis.qzv

# Make model
dir=Strat_group_ref # For Sediment: Location_description
mkdir ${dir}

qiime songbird multinomial \
    --i-table feature-table.qza \
    --m-metadata-file metadata_for_qiime2_fixed.txt \
    --p-formula "C(Strat_group, Treatment('Unstratified'))" \ # For Sediment: --p-formula "Location_description"
    --p-epochs 10000 \
    --p-differential-prior 0.5 \
    --p-summary-interval 1 \
    --p-num-random-test-examples 5 \
    --o-differentials ${dir}/differentials.qza \
    --o-regression-stats ${dir}/regression-stats.qza \
    --o-regression-biplot ${dir}/regression-biplot.qza \
    --verbose

# Make null model
null_dir=null_model
mkdir ${null_dir}

qiime songbird multinomial \
    --i-table feature-table.qza \
    --m-metadata-file metadata_for_qiime2_fixed.txt \
    --p-formula "1" \
    --p-epochs 10000 \
    --p-differential-prior 0.5 \
    --p-summary-interval 1 \
    --o-differentials ${null_dir}/null-diff.qza \
    --o-regression-stats ${null_dir}/null-stats.qza \
    --o-regression-biplot ${null_dir}/null-biplot.qza \
    --p-num-random-test-examples 5 \
    --verbose

# Visualize the first model's regression stats and the null model's
qiime songbird summarize-paired \
    --i-regression-stats ${dir}/regression-stats.qza \
    --i-baseline-stats ${null_dir}/null-stats.qza \
    --o-visualization ${dir}/paired-summary.qzv

qiime tools view ${dir}/paired-summary.qzv

# Qurro
qiime qurro differential-plot \
    --i-ranks ${dir}/differentials.qza \
    --i-table feature-table.qza \
    --m-sample-metadata-file metadata_for_qiime2_fixed.txt \
    --m-feature-metadata-file tax_for_qiime2.txt \
    --verbose \
    --o-visualization ${dir}/qurro_plot_q2.qzv

qiime tools view ${dir}/qurro_plot_q2.qzv

# Export Songbird data
qiime metadata tabulate \
    --m-input-file ${dir}/differentials.qza \
    --o-visualization ${dir}/differentials-viz.qzv

qiime tools export \
  --input-path ${dir}/differentials-viz.qzv \
  --output-path ${dir}/differentials

#### Songbird-picrust2 top categories ####

# all pathways
use_pathway <- pred_pathway[,-1]
head(use_pathway)

# load pathways to focus on (csv file has columns for description, Order, broad_categorization1, broad_categorization2, broad_categorization3)
pathway_list <- read.delim(<path to csv file with select pathways>, sep = ",", stringsAsFactors=FALSE, fileEncoding="latin1")
rownames(pathway_list) <- pathway_list[,1]
pathway_list[,1] <- NULL
pathway_list

# subset pathways
focus_pathways <- use_pathway %>% filter(rownames(use_pathway) %in% rownames(pathway_list))
nrow(focus_pathways)
nrow(pathway_list)

focus_pathways <- as.data.frame(t(focus_pathways))
focus_pathways <- tibble::rownames_to_column(focus_pathways, "SampleID")
head(focus_pathways)

# setup denominator
pred_pathway_denominator <- use_pathway %>% filter(!rownames(use_pathway) %in% rownames(pathway_list))
head(pred_pathway_denominator)

# sum denominator
sum_denominator <- as.data.frame(colSums(pred_pathway_denominator))
colnames(sum_denominator) <- column_names
sum_denominator <- tibble::rownames_to_column(sum_denominator, "SampleID")
head(sum_denominator)

# calculate natural log ratios (Note: log is natural log in R)
df1 <- focus_pathways
df2 <- sum_denominator
df3 <- cbind(df1[1], log(df1[, -1] / df2[match(df1$SampleID, df2$SampleID), -1]))
head(df3)

# add descriptions
rownames(df3) <- df3$SampleID
df3$SampleID <- NULL
df3 <- as.data.frame(t(df3))
df3$description <- pathway_list$description[match(rownames(df3), rownames(pathway_list))]
df3$broad_categorization1 <- pathway_list$broad_categorization1[match(rownames(df3), rownames(pathway_list))]
df3$broad_categorization2 <- pathway_list$broad_categorization2[match(rownames(df3), rownames(pathway_list))]
df3$broad_categorization3 <- pathway_list$broad_categorization3[match(rownames(df3), rownames(pathway_list))]
df3$Order <- pathway_list$Order[match(rownames(df3), rownames(pathway_list))]
df3 <- tibble::rownames_to_column(df3, "pathway")
head(df3)

# add stratificaiton group layer
df3_melt <- reshape2::melt(df3, id=c("pathway", "description", "broad_categorization1", "broad_categorization2", "broad_categorization3", "Order"))
colnames(df3_melt) <- c("pathway", "description", "broad_categorization1", "broad_categorization2", "broad_categorization3", "Order", "SampleID", "Natural_Log_Ratio")
df3_melt$SampleID <- gsub("\\.", "-", df3_melt$SampleID)
df3_melt$Strat_group <- DATA_PHYLOSEQ_FIXED$Strat_group[match(df3_melt$SampleID, DATA_PHYLOSEQ_FIXED$SampleID)]
df3_melt$strat_time <- ordi_data$strat_time[match(df3_melt$SampleID, ordi_data$SampleID)]

# order the strat group (Water Samples)
df3_melt$Strat_group <- gsub('Unstratified', 'Incipient Stratification', df3_melt$Strat_group)
df3_melt$Strat_group <- factor(df3_melt$Strat_group, levels=c("Incipient Stratification","Epilimnion","Thermocline", "Hypolimnion"), ordered = TRUE)

# order the strat group (Sediment Samples)
df3_melt$Location_description <- factor(df3_melt$Location_description, levels=c("Inlet","Outlet"), ordered = TRUE)
df3_melt$Depth_desc <- factor(df3_melt$Depth_desc, levels=c("Top 5 cm","5-10 cm","10-21 cm"), ordered = TRUE)

# remove Inf values
df3_melt <- df3_melt[is.finite(df3_melt$Natural_Log_Ratio),]
head(df3_melt)
write.csv(df3_melt, "top_contributing_pathways_Limited2.csv")

#### Plot Figure 6 ####
## Prepare for plotting heatmap (Water Samples)
df3_melt_avg <- df3_melt %>%
  group_by(pathway, description, broad_categorization1, broad_categorization2, broad_categorization3, Order, Strat_group, strat_time) %>%
  summarise_at(vars(Natural_Log_Ratio), list(mean = mean)) %>%
  ungroup()
df3_melt_avg <- as.data.frame(df3_melt_avg)
df3_melt_avg$broad_categorization2 <- factor(df3_melt_avg$broad_categorization2, levels=c('Carbon', 'Iron'), ordered = TRUE)
df3_melt_avg$Order <- factor(df3_melt_avg$Order, levels=c('1', '2', '3', '4', '5', '6' ,'7', '8', '9', '10', '11', '12', '13', '14'))
df3_melt_avg$broad_categorization3 <- factor(df3_melt_avg$broad_categorization3, levels=c('Methanogenesis', 'CO2 Fixation', 'C1 compound utilization', 'Fermentation', 'Aromatic Compounds', 'Iron microbes', 'Siderophore biosynthesis', 'Sulfur'), ordered = TRUE)
df3_melt_avg$broad_categorization1 <- gsub('Iron reducers', 'Iron reduction: Fe(III) to Fe(II)', df3_melt_avg$broad_categorization1)
df3_melt_avg$broad_categorization1 <- gsub('Magnetotactic bacteria', 'Magnetosome formation', df3_melt_avg$broad_categorization1)
df3_melt_avg$broad_categorization1 <- factor(df3_melt_avg$broad_categorization1, levels = unique(df3_melt_avg$broad_categorization1[order(df3_melt_avg$Order)]))
head(df3_melt_avg)

## Prepare for plotting heatmap (Sediment samples)
df3_melt_avg <- df3_melt %>%
  group_by(pathway, desc, broad_categorization1, broad_categorization2, broad_categorization3, Order, Location_description) %>%
  summarise_at(vars(Natural_Log_Ratio), list(mean = mean)) %>%
  ungroup()
df3_melt_avg <- as.data.frame(df3_melt_avg)
df3_melt_avg$broad_categorization2 <- factor(df3_melt_avg$broad_categorization2, levels=c('Carbon', 'Iron'), ordered = TRUE)
df3_melt_avg$Order <- factor(df3_melt_avg$Order, levels=c('1', '2', '3', '4', '5', '6' ,'7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24'))
df3_melt_avg$broad_categorization3 <- factor(df3_melt_avg$broad_categorization3, levels=c('Methanogens', 'CO2 Fixation', 'Carbohydrate Degradation', 'Fermentation', 'Alcohol Degradation', 'Aldehyde Degradation', 'Amine and Polyamine Degradation', 'Aromatic Compound Degradation', 'Siderophore biosynthesis'), ordered = TRUE)
df3_melt_avg$broad_categorization1 <- factor(df3_melt_avg$broad_categorization1, levels = unique(df3_melt_avg$broad_categorization1[order(df3_melt_avg$Order)]))
head(df3_melt_avg)

# Set up theme
plot_theme <- theme(panel.background = element_rect(fill = "white", colour = "white", size = 0, linetype = "solid"),
    panel.border = element_rect(colour="white", size=0, fill=NA),
    strip.background=element_rect(fill='white', colour='white', size = 0),
    strip.text = element_text(face="bold", size=10),
    strip.text.y = element_text(angle = 360, hjust=0),
    panel.spacing.x=unit(0.5, "lines"),
    panel.grid.major = element_line(size = 0),
    panel.grid.minor = element_line(size = 0),
    axis.text = element_text(size=10, colour="black"),
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
    legend.position="right",
    legend.key = element_rect(fill = "white"),
    legend.title = element_text(face="bold", size=10),
    legend.text = element_text(size=10))
plot_guide <- guides(fill = guide_colourbar(frame.linewidth = 1, frame.colour = "black", ticks = TRUE, ticks.colour = "black", ticks.linewidth = 1, reverse=F))

# Plot water samples
ggplot(df3_melt_avg, aes(Strat_group, broad_categorization1)) +
    geom_tile(aes(fill = mean, width=0.9), colour = "black", size=0.5) +
    scale_fill_gradientn(colors=c( "#7a024e", "#fdece9", "#b0d66d", "#34a0a4", "#02303f"), breaks=c(-8, -10, -12, -14, -16)) +
    facet_grid(broad_categorization3 ~ ., scales = "free", space = "free") +
    labs(fill="Natural Log Ratio\n(Average)") +
    scale_y_discrete(limits=rev) + plot_theme + plot_guide

save_file_plot <- paste("top_contributing_pathways_heatmap_water.svg", sep="") #change the file name if need to
ggsave(save_file_plot, path = songbird_picrust, scale = 1, width = 7, height =6, units = c("in"), dpi = 300)

# Plot sediment samples
ggplot(df3_melt_avg, aes(Location_description, broad_categorization1)) +
    geom_tile(aes(fill = mean, width=0.9), colour = "black", size=0.5) +
    scale_fill_gradientn(colors=c("#213e1b", "#9be564", "#fee231", "#e78e23", "#341209"), breaks=c(-6, -8, -10, -12, -14)) +
    facet_grid(broad_categorization3 ~ ., scales = "free", space = "free") +
    labs(fill="Natural Log Ratio\n(Average)") +
    scale_y_discrete(limits=rev) + plot_theme + plot_guide

save_file_plot <- paste("top_contributing_pathways_heatmap_sediment.svg", sep="") #change the file name if need to
ggsave(save_file_plot, path = songbird_picrust, scale = 1, width = 7, height =6, units = c("in"), dpi = 300)
