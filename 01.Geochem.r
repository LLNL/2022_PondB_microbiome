# 01.Geochem.r
# Figure 2, Figure S2, Figure S3, Table S1, Table S2, Table S3

#### Load Packages ####
packages <- c("plyr", "dplyr", "ggplot2", "reshape2", "tidyr", "fuzzyjoin", "mltools", "cutr", "plotly", "utils", "stringr", "data.table", "ggpubr", "cowplot", "readxl")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages], dependencies=TRUE)
}

# Load Packages
invisible(lapply(packages, library, character.only = TRUE))

print(paste("R", getRversion()))
print("-------------")
for (package_name in sort(loadedNamespaces())) {
    print(paste(package_name, packageVersion(package_name)))
}

#### Set Paths ####
mypath = <path to main dir>
setwd(mypath)

folder_ff <- c("Figures_PondBMicrobePaper")
dir.create(folder_ff, showWarning=FALSE)

#### Get Data ####
Fe_data <- read.delim(paste0(mypath, folder_ff, "/Fe_data.csv"), row.names=1, sep = ",", stringsAsFactors=FALSE, fileEncoding="latin1")
SO4_data <- read.delim(paste0(mypath, folder_ff, "/SO4_data.csv"), row.names=1, sep = ",", stringsAsFactors=FALSE, fileEncoding="latin1")
TOC_data <- read.delim(paste0(mypath, folder_ff, "/TOC_data.csv"), row.names=1, sep = ",", stringsAsFactors=FALSE, fileEncoding="latin1")
CTD_data <- read.delim(paste0(mypath, folder_ff, "/CTD_data.csv"), row.names=1, sep = ",", stringsAsFactors=FALSE, fileEncoding="latin1")

#### Organize Data ####
# Fix locations
Fe_data$Location <- factor(Fe_data$Location, ordered = TRUE, levels = c("Inlet", "Inlet-Middle", "Middle", "Middle-Outlet", "Outlet"))
CTD_data$Loc <- factor(CTD_data$Loc, ordered = TRUE, levels = c("Inlet", "Inlet-Middle", "Middle", "Middle-Outlet", "Outlet"))
SO4_data$Loc <- factor(SO4_data$Loc, ordered = TRUE, levels = c("Inlet", "Inlet-Middle", "Middle", "Middle-Outlet", "Outlet"))
TOC_data$Loc <- factor(TOC_data$Loc, ordered = TRUE, levels = c("Inlet", "Inlet-Middle", "Middle", "Middle-Outlet", "Outlet"))

# Remove weird characters
CTD_data$Measurement <- gsub('Â', '', CTD_data$Measurement)

# Convert to mM
Fe_data$Average_mM <- Fe_data$Average_mg_L / 55.845
Fe_data$SD_mM <- Fe_data$SD / 55.845
SO4_data$Average_uM <- (SO4_data$Average / 96.06)*1000
SO4_data$SD_uM <- (SO4_data$SD / 96.06)*1000

# Modify Date to stratification time
Fe_data$Strat_time <- Fe_data$Date
Fe_data$Strat_time <- gsub('Oct-20', 'Stratified (~7 months)', Fe_data$Strat_time)
Fe_data$Strat_time <- gsub('Jun-19', 'Stratified (~3 months)', Fe_data$Strat_time)
Fe_data$Strat_time <- gsub('Mar-20', 'Incipient Stratification', Fe_data$Strat_time)
Fe_data$Strat_time <- gsub('Feb-21', 'Unstratified (~3 months)', Fe_data$Strat_time)
Fe_data$Strat_time <- gsub('Dec-20', 'Unstratified (~1 month)', Fe_data$Strat_time)
Fe_data$Strat_time <- factor(Fe_data$Strat_time, ordered = TRUE, levels = c("Unstratified (~1 month)", "Unstratified (~3 months)", "Incipient Stratification", "Stratified (~3 months)", "Stratified (~7 months)"))
SO4_data$Strat_time <- SO4_data$Month
SO4_data$Strat_time <- gsub('June', 'Stratified (~3 months)', SO4_data$Strat_time)
SO4_data$Strat_time <- gsub('March', 'Incipient Stratification', SO4_data$Strat_time)
SO4_data$Strat_time <- factor(SO4_data$Strat_time, ordered = TRUE, levels = c("Incipient Stratification", "Stratified (~3 months)"))
TOC_data$Strat_time <- TOC_data$Month
TOC_data$Strat_time <- gsub('June', 'Stratified (~3 months)', TOC_data$Strat_time)
TOC_data$Strat_time <- gsub('March', 'Incipient Stratification', TOC_data$Strat_time)
TOC_data$Strat_time <- factor(TOC_data$Strat_time, ordered = TRUE, levels = c("Incipient Stratification", "Stratified (~3 months)"))
CTD_data$Strat_time <- gsub('7 months est.', '~7 months', CTD_data$Strat_time)
CTD_data$Strat_time <- gsub('3 months est.', '~3 months', CTD_data$Strat_time)
CTD_data$Strat_time <- gsub('Spring turnover', 'Incipient Stratification', CTD_data$Strat_time)
CTD_data$Strat_time <- gsub('3 months est.', '~3 months', CTD_data$Strat_time)
CTD_data$Strat_time <- gsub('1 month est.', '~1 month', CTD_data$Strat_time)
CTD_data$Strat_time <- factor(CTD_data$Strat_time, ordered = TRUE, levels = c("Unstratified (~1 month)", "Unstratified (~3 months)", "Incipient Stratification", "Stratified (~3 months)", "Stratified (~7 months)"))

# Remove rows with NA in Strat_time
Fe_data <- Fe_data[complete.cases(Fe_data[ , 7]),]
CTD_data <- CTD_data[complete.cases(CTD_data[ , 9]),]
SO4_data <- SO4_data[complete.cases(SO4_data[ , 11]),]
TOC_data <- TOC_data[complete.cases(TOC_data[ , 9]),]

# Order the data
Fe_data <- Fe_data[order(Fe_data$Location, Fe_data$Strat_time),]
CTD_data <- CTD_data[order(CTD_data$Loc, CTD_data$Strat_time),]
SO4_data <- SO4_data[order(SO4_data$Loc, SO4_data$Strat_time),]
TOC_data <- TOC_data[order(TOC_data$Loc, TOC_data$Strat_time),]

#### Plot Figure 2 ####
# Plot configurations
plot_color <- scale_colour_manual(values=c("Unstratified (~1 month)" = "#2a9d8f", "Unstratified (~3 months)" = "#264653", "Incipient Stratification" = "#ce4257", "Stratified (~3 months)" = "#582f0e", "Stratified (~7 months)" = "#fb8b24"))
plot_shape <- scale_shape_manual(values=c("Unstratified (~1 month)" = 17, "Unstratified (~3 months)" = 17, "Incipient Stratification" = 8, "Stratified (~3 months)" = 16, "Stratified (~7 months)" = 16))
plot_theme <- theme(panel.background = element_rect(fill = "white", colour = "black", size = 1, linetype = "solid"),
                    panel.border = element_rect(colour="black", size=1, fill=NA),
                    strip.background=element_rect(fill='white', colour='white'),
                    strip.text = element_text(color="black", size=20, face="bold"),
                    panel.grid.major = element_line(size = 0),
                    panel.grid.minor = element_line(size = 0),
                    axis.text = element_text(size=20, colour="black"),
                    axis.title.x = element_text(color="black", size=20, face="bold"),
                    axis.title.y = element_text(color="black", size=20, face="bold"),
                    legend.position="right",
                    legend.key = element_rect(fill = "white"),
                    legend.title = element_blank(),
                    legend.text = element_text(size=20),
                    panel.spacing = unit(2, "lines"),
                    plot.title = element_blank())
plot_guides <- guides(colour=guide_legend(ncol=1))

# Plot Total Fe
cat = c("Total")
cat2 = c("Middle")
Fe_data_subTFe <- subset(Fe_data, Fe_data$Total_Filtered == cat)
Fe_data_subTFe <- subset(Fe_data_subTFe, Fe_data_subTFe$Location == cat2)
fig_totalFe <- ggplot(Fe_data_subTFe, aes(x = Depth_m, y = Average_mM, colour=Strat_time, shape=Strat_time)) +
        geom_errorbar(aes(ymin = Average_mM - SD_mM, ymax = Average_mM + SD_mM), width=0.1) +
        geom_line() +
        geom_point(size = 3) +
        geom_vline(xintercept=2, linetype="dashed", color = "black", size=0.5) +
        geom_vline(xintercept=6.5, linetype="dashed", color = "black", size=0.5) +
        geom_vline(xintercept=0, linetype="solid", color = "black", size=0.5) +
        coord_flip() +
        labs(x = "Depth (m)", y="Total Fe (mM)") +
        scale_x_reverse(lim=c(11,-0.5), breaks = seq(10, 0, by = -2)) +
        scale_y_continuous(lim=c(0,0.7), breaks = seq(0,0.6, by = 0.2)) + #, labels = scales::scientific
        plot_theme + plot_guides + plot_color + plot_shape

# Plot Filtered vs Total Fe
cat = c("Middle")
cat2 = c("Stratified (~7 months)")
Fe_data_subFFe <- subset(Fe_data, Fe_data$Location == cat)
Fe_data_subFFe <- subset(Fe_data_subFFe, Fe_data_subFFe$Strat_time == cat2)
fig_filteredFe <- ggplot(Fe_data_subFFe, aes(x = Depth_m, y = Average_mM, colour=Strat_time, shape=Total_Filtered, linetype=Total_Filtered)) +
        geom_errorbar(aes(ymin = Average_mM - SD_mM, ymax = Average_mM + SD_mM), width=0.1) +
        geom_line() +
        geom_point(size = 3) +
        geom_vline(xintercept=2, linetype="dashed", color = "black", size=0.5) +
        geom_vline(xintercept=6.5, linetype="dashed", color = "black", size=0.5) +
        geom_vline(xintercept=0, linetype="solid", color = "black", size=0.5) +
        coord_flip() +
        labs(x = "Depth (m)", y="Total Fe (mM)") +
        scale_x_reverse(lim=c(11,-0.5), breaks = seq(10, 0, by = -2)) +
        scale_y_continuous(lim=c(0,0.7), breaks = seq(0,0.6, by = 0.2)) + #, labels = scales::scientific
        plot_theme + plot_color +
        scale_shape_manual(values=c("Total" = 16, "Filtered" = 1)) +
        scale_linetype_manual(values=c("Total" = "solid", "Filtered" = "dashed")) +
        guides(shape = guide_legend(override.aes = list(colour = "#fb8b24") ) )

# Plot Temp
cat = c("Temperature (ºC)")
cat2 = c("Middle")
CTD_data_sub <- subset(CTD_data, CTD_data$Measurement == cat)
CTD_data_sub <- subset(CTD_data_sub, CTD_data_sub$Loc == cat2)
figure_temp <- ggplot(CTD_data_sub, aes(x = Depth_m, y = Average, colour=Strat_time, shape=Strat_time)) +
        geom_errorbar(aes(ymin = Average - SD, ymax = Average + SD), width=0.1) +
        geom_line() +
        geom_point(size = 3) +
        geom_vline(xintercept=2, linetype="dashed", color = "black", size=0.5) +
        geom_vline(xintercept=6.5, linetype="dashed", color = "black", size=0.5) +
        geom_vline(xintercept=0, linetype="solid", color = "black", size=0.5) +
        coord_flip() +
        labs(x = "Water Column Depth (m)", y=cat) +
        scale_x_reverse(lim=c(11,-0.5), breaks = seq(10, 0, by = -2)) +
        scale_y_continuous(lim=c(8,32), breaks = seq(10,30, by = 5)) + #, labels = scales::scientific
        plot_theme + plot_guides + plot_color + plot_shape

# Plot DO
cat = c("Dissolved Oxygen (%)")
cat2 = c("Middle")
CTD_data_sub <- subset(CTD_data, CTD_data$Measurement == cat)
CTD_data_sub <- subset(CTD_data_sub, CTD_data_sub$Loc == cat2)
figure_DO <- ggplot(CTD_data_sub, aes(x = Depth_m, y = Average, colour=Strat_time, shape=Strat_time)) +
        geom_errorbar(aes(ymin = Average - SD, ymax = Average + SD), width=0.1) +
        geom_line() +
        geom_point(size = 3) +
        geom_vline(xintercept=2, linetype="dashed", color = "black", size=0.5) +
        geom_vline(xintercept=6.5, linetype="dashed", color = "black", size=0.5) +
        geom_vline(xintercept=0, linetype="solid", color = "black", size=0.5) +
        coord_flip() +
        labs(x = "Depth (m)", y=cat) +
        scale_x_reverse(lim=c(11,-0.5), breaks = seq(10, 0, by = -2)) +
        scale_y_continuous(lim=c(-0.5,105), breaks = seq(0,100, by = 20)) + #, labels = scales::scientific
        plot_theme + plot_guides + plot_color + plot_shape

# Combine
all <- plot_grid(figure_temp + theme(legend.position="none"),
                 figure_DO + theme(legend.position="none", axis.title.y = element_blank()),
                 fig_totalFe + theme(legend.position="none", axis.title.y = element_blank()),
                 fig_filteredFe + theme(legend.position="none", axis.title.y = element_blank()),
                 ncol=4, align = "v", axis="b")

save_file <- paste("Figure_02.pdf", sep="")
ggsave(save_file, path = folder_ff, plot = all, scale = 1, width = 15, height = 5, units = c("in"), dpi = 300)

#### Plot Figure S2 ####

# Plot Total Fe
cat = c("Total")
Fe_data_subTFe <- subset(Fe_data, Fe_data$Total_Filtered == cat)
fig_totalFe <- ggplot(Fe_data_subTFe, aes(x = Depth_m, y = Average_mM, colour=Strat_time, shape=Strat_time)) +
        geom_errorbar(aes(ymin = Average_mM - SD_mM, ymax = Average_mM + SD_mM), width=0.1) +
        geom_line() +
        geom_point(size = 3) +
        geom_vline(xintercept=2, linetype="dashed", color = "black", size=0.5) +
        geom_vline(xintercept=6.5, linetype="dashed", color = "black", size=0.5) +
        geom_vline(xintercept=0, linetype="solid", color = "black", size=0.5) +
        coord_flip() +
        labs(x = "Depth (m)", y="Total Fe (mM)") +
        facet_wrap( ~ Location, nrow=1) +
        scale_x_reverse(lim=c(11,-0.5), breaks = seq(10, 0, by = -2)) +
        scale_y_continuous(lim=c(0,0.7), breaks = seq(0,0.6, by = 0.2)) +
        plot_theme + plot_guides + plot_color + plot_shape

# Plot TOC
fig_TOC <- ggplot(TOC_data, aes(x = Depth_m, y = Average, colour=Strat_time, shape=Strat_time)) +
        geom_errorbar(aes(ymin = Average - SD, ymax = Average + SD), width=0.1) +
        geom_line() +
        geom_point(size = 3) +
        geom_vline(xintercept=2, linetype="dashed", color = "black", size=0.5) +
        geom_vline(xintercept=6.5, linetype="dashed", color = "black", size=0.5) +
        geom_vline(xintercept=0, linetype="solid", color = "black", size=0.5) +
        coord_flip() +
        labs(x = "Depth (m)", y="TOC (mg/L)") +
        facet_wrap( ~ Loc, nrow=1) +
        scale_x_reverse(lim=c(11,-0.5), breaks = seq(10, 0, by = -2)) +
        scale_y_continuous(lim=c(0,13), breaks = seq(0,12, by = 2)) +
        plot_theme + plot_guides + plot_color + plot_shape

# Plot SO4
fig_SO4 <- ggplot(SO4_data, aes(x = Depth_m, y = Average_uM, colour=Strat_time, shape=Strat_time)) +
        geom_errorbar(aes(ymin = Average_uM - SD_uM, ymax = Average_uM + SD_uM), width=0.1) +
        geom_line() +
        geom_point(size = 3) +
        geom_vline(xintercept=2, linetype="dashed", color = "black", size=0.5) +
        geom_vline(xintercept=6.5, linetype="dashed", color = "black", size=0.5) +
        geom_vline(xintercept=0, linetype="solid", color = "black", size=0.5) +
        coord_flip() +
        labs(x = "Depth (m)", y="Sulfate (µM)") +
        facet_wrap( ~ Loc, nrow=1) +
        scale_x_reverse(lim=c(11,-0.5), breaks = seq(10, 0, by = -2)) +
        scale_y_continuous(lim=c(0,7), breaks = seq(0,7, by = 2)) +
        plot_theme + plot_guides + plot_color + plot_shape

# Plot Temp
cat = c("Temperature (ºC)")
CTD_data_sub <- subset(CTD_data, CTD_data$Measurement == cat)
figure_temp <- ggplot(CTD_data_sub, aes(x = Depth_m, y = Average, colour=Strat_time, shape=Strat_time)) +
        geom_errorbar(aes(ymin = Average - SD, ymax = Average + SD), width=0.1) +
        geom_line() +
        geom_point(size = 3) +
        geom_vline(xintercept=2, linetype="dashed", color = "black", size=0.5) +
        geom_vline(xintercept=6.5, linetype="dashed", color = "black", size=0.5) +
        geom_vline(xintercept=0, linetype="solid", color = "black", size=0.5) +
        coord_flip() +
        labs(x = "Depth (m)", y=cat) +
        facet_wrap( ~ Loc, nrow=1) +
        scale_x_reverse(lim=c(11,-0.5), breaks = seq(10, 0, by = -2)) +
        scale_y_continuous(lim=c(8,32), breaks = seq(10,30, by = 5)) + #, labels = scales::scientific
        plot_theme + plot_guides + plot_color + plot_shape

# Plot DO
cat = c("Dissolved Oxygen (%)")
CTD_data_sub <- subset(CTD_data, CTD_data$Measurement == cat)
figure_DO <- ggplot(CTD_data_sub, aes(x = Depth_m, y = Average, colour=Strat_time, shape=Strat_time)) +
        geom_errorbar(aes(ymin = Average - SD, ymax = Average + SD), width=0.1) +
        geom_line() +
        geom_point(size = 3) +
        geom_vline(xintercept=2, linetype="dashed", color = "black", size=0.5) +
        geom_vline(xintercept=6.5, linetype="dashed", color = "black", size=0.5) +
        geom_vline(xintercept=0, linetype="solid", color = "black", size=0.5) +
        coord_flip() +
        labs(x = "Depth (m)", y=cat) +
        facet_wrap( ~ Loc, nrow=1) +
        scale_x_reverse(lim=c(11,-0.5), breaks = seq(10, 0, by = -2)) +
        scale_y_continuous(lim=c(-0.5,110), breaks = seq(0,100, by = 20)) + #, labels = scales::scientific
        plot_theme + plot_guides + plot_color + plot_shape

# Combine
all <- plot_grid(figure_temp + theme(legend.position="none"),
                 figure_DO + theme(legend.position="none", strip.text = element_blank(), strip.background = element_blank()),
                 fig_TOC + theme(legend.position="none", strip.text = element_blank(), strip.background = element_blank()),
                 fig_totalFe + theme(legend.position="none", strip.text = element_blank(), strip.background = element_blank()),
                 fig_SO4 + theme(legend.position="none", strip.text = element_blank(), strip.background = element_blank()),
                 ncol=1, align = "v", axis="b")

save_file <- paste("Figure_S2.pdf", sep="")
ggsave(save_file, path = folder_ff, plot = all, scale = 1, width = 14, height = 17, units = c("in"), dpi = 300)

#### Prepare to Plot Figure S3 ####
# Set path
mypath = <path to CTD data for contour plots>
setwd(mypath)

folder_f <- c("summary_figures")
folder_c <- c("contour_plots")
dir.create(folder_f, showWarning=FALSE)
dir.create(folder_c, showWarning=FALSE)

# Retrieve all the summary table files
temp = list.files(pattern="*.txt")

print("List of all files:")
temp

# Create list of data frame names without the ".txt" part
names <- sapply(strsplit(temp, ".txt"), `[`, 1)

for (i in names) {
    filepath <- file.path(paste(i, ".txt",sep=""))
    var_name <- gsub("prelim_table_v3_","contour_",i)
    assign(var_name, read.delim(filepath, sep = "\t", stringsAsFactors=FALSE, fileEncoding="latin1"))
}

# Set location to plot
Location <- c("Loc3")
df_list <- as.list(ls(pattern=glob2rx(paste("contour*",Location,sep=""))))

print("List of all the data that will be plotted in the contour:")
df_list

print("The number of columns in these dataframes need to be the same")
for(i in df_list){
    data <- as.data.frame(get(i))
    data$Measurement <- gsub('Â', '', data$Measurement)
    unmelt <- reshape(data, idvar = c("Depth_m", "Date"), timevar = "Measurement", direction = "wide")
    print(paste("Num. of columns in ", i, ": ", ncol(unmelt), sep=""))
    assign(paste(i,"_unmelt", sep=""), unmelt)
}

list_unmelt <- mget(ls(pattern=glob2rx(paste("contour*",Location,"_unmelt",sep=""))))

combo <- do.call(rbind, list_unmelt) # Combine all the dataframes
row.names(combo) <- NULL # Remove the rownames

# Change the Date to a Date class
combo$Date <- factor(combo$Date)
combo$Date <- as.Date(combo$Date, format = "%m/%d/%y")

# Check the data
head(combo)

# get the max and min of depth
max_y <- ceiling(max(combo$Depth_m))
min_y <- floor(min(combo$Depth_m))

# get the earliest and latest date
max_x <- max(combo$Date)
min_x <- min(combo$Date)

# Configuration of the x-axis and y-axis
x_config <- list(
    tickcolor = toRGB("black"), tickfont = list(size = 20), tickwidth = 1,
    title="", titlefont = list(size = 20),
    color = "black",
    showline = TRUE, linecolor = toRGB("black"), linewidth = 1,
    mirror="allticks",
    range = c(min_x, max_x))

y_config <- list(
    tickcolor = toRGB("black"), tickwidth = 1, tickfont = list(size = 20),
    title = "Depth (m)", titlefont = list(size = 20),
    color = "black",
    showline = TRUE,
    linecolor = toRGB("black"), linewidth = 1,
    mirror="allticks",
    rangemode = 'tozero', autorange="reversed")

# Configuration of the legend
leg_config <- list(x = 100, y = 0.5)

parameter <- c("Average.Chlorophyll_A_Fluorescence_RFU", "Average.Conductivity_µS_cm", "Average.Dissolved_Oxygen_Sat", "Average.ORP_mV", "Average.pH", "Average.Specific_Conductivity..µS.cm", "Average.Temperature_C", "Average.Total_Dissolved_Solids_ppt")

max <- as.data.frame(apply(combo[,parameter], 2, max))
min <- as.data.frame(apply(combo[,parameter], 2, min))
max_min <- cbind(min, max)
colnames(max_min) <- c("Min", "Max")
format(max_min, digits = 2, scientific=F)

# Check the measurements that can be plotted. These are the column names
colnames(combo)

# Save the combined file
filename <- paste(folder_c, "/contour_plots_table_", Location, ".txt", sep="")
write.table(combo, filename, sep="\t")

# Load new combined file
filepath <- paste(folder_c, "/contour_plots_table_", Location, ".txt", sep="")
combo <- read.delim(filepath, sep = "\t", stringsAsFactors=FALSE, fileEncoding="latin1")

# Correct column names
names(combo) <- gsub(x = names(combo), pattern = "Â", replacement = "")

# Change the Date to a Date class
combo$Date <- factor(combo$Date)
combo$Date <- as.Date(combo$Date, format = "%m/%d/%y")

# Check new file loaded correctly
head(combo)

#### Plot the contours: Figure S3 ####

### TEMPERATURE ###

## Change these values ##
min <- 10
max <- 30
interval <- 0.5
numcontours <- 15
legend_interval <- 5
#########################

# Configuration of the contours
contour_fig <- list(
    coloring = 'heatmap',
    start = min, end = max, size = interval,
    showlabels = TRUE, labelfont = list(size = 15, color = 'white'),
    showlines = FALSE))

palette <- colorRampPalette(c("#071183", "#114c79", "#0d7c9b", "#bfab6b", "#f9aa0f", "#a60808"))

# Configuration of the contour lines
line_fig <- list(smoothing = 1.3, width=0)

# Contour plot
Temp <- plot_ly(combo,
        x = ~Date,
        y = ~Depth_m,
        z = ~Average.Temperature_C,
        zmin = min,
        zmax = max,
        type = "contour",
        autocontour = FALSE,
        ncontours = numcontours,
        #contours = contour_fig,
        contours = contour_fig,
        #colorscale = 'Portland',
        colorscale = cbind(seq(0, 1, by=1/(length(combo) - 1)), palette(length(combo))),
        line = line_fig,
        colorbar = list(title = list(text = "Temperature\n(ºC)", font = list(size = 15, color = 'black')),
                      tickfont = list(size=15, color='black'), tick0 = min, dtick = legend_interval,
                      ticks = "outside", ticklen = 5, tickcolor = 'black'),
        width = 1000, height = 500) %>%
        layout(yaxis = y_config, xaxis = x_config, legend = leg_config)

# Save the file
filename <- paste(folder_c, "/contour_plots_temp_", Location, ".pdf", sep="")
orca(Temp, filename, width = 4 * 300, height = 2 * 300)

### ORP ###

## Change these values ##
min <- -500
max <- 350
interval <- 1
numcontours <- 15
legend_interval <- 100
#########################

# Configuration of the contours
contour_fig <- list(
    coloring = 'heatmap',
    start = min, end = max, size = interval,
    showlabels = TRUE, labelfont = list(size = 15, color = 'white'),
    showlines = FALSE))
    
# Configuration of the contour lines
line_fig <- list(smoothing = 1.3, width=0)

# Contour plot
ORP <- plot_ly(combo,
        x = ~Date,
        y = ~Depth_m,
        z = ~Average.ORP_mV,
        zmin = min,
        zmax = max,
        type = "contour",
        autocontour = TRUE,
        ncontours = numcontours,
        contours = contour_fig,
        colorscale = 'Viridis',
        line = line_fig,
        colorbar = list(title = list(text = "ORP\n(mV)", font = list(size = 15, color = 'black')),
                      tickfont = list(size=15, color='black'), tick0 = min, dtick = legend_interval,
                      ticks = "outside", ticklen = 5, tickcolor = 'black'),
        autocontour = TRUE,
        width = 1000, height = 500) %>%
        layout(yaxis = y_config, xaxis = x_config, legend = leg_config)

# Save the file
filename <- paste(folder_c, "/contour_plots_ORP_", Location, ".pdf", sep="")
orca(ORP, filename, width = 4 * 300, height = 2 * 300)

### CHLOROPHYLL A ###

## Change these values ##
min <- 0
max <- 50
interval <- 1
numcontours <- 15
legend_interval <- 10
#########################

# Configuration of the contours
contour_fig <- list(
    coloring = 'heatmap',
    start = min, end = max, size = interval,
    showlabels = TRUE, labelfont = list(size = 15, color = 'white'),
    showlines = FALSE))
    
# Configuration of the contour lines
line_fig <- list(smoothing = 1.3, width=0)

# Contour plot
ChlA <- plot_ly(combo,
        x = ~Date,
        y = ~Depth_m,
        z = ~Average.Chlorophyll_A_Fluorescence_RFU,
        zmin = min,
        zmax = max,
        type = "contour",
        autocontour = TRUE,
        ncontours = numcontours,
        contours = contour_fig,
        #colorscale = 'Portland',
        line = line_fig,
        colorbar = list(title = list(text = "Chlorophyll A\n(RFU)", font = list(size = 15, color = 'black')),
                      tickfont = list(size=15, color='black'), tick0 = min, dtick = legend_interval,
                      ticks = "outside", ticklen = 5, tickcolor = 'black'),
        autocontour = TRUE,
        width = 1000, height = 500) %>%
        layout(yaxis = y_config, xaxis = x_config, legend = leg_config)

# Save the file
filename <- paste(folder_c, "/contour_plots_ChlA_", Location, ".pdf", sep="")
orca(ChlA, filename, width = 4 * 300, height = 2 * 300)

### DISSOLVED OXYGEN ###

## Change these values ##
min <- 0
max <- 100
interval <- 1
numcontours <- 15
legend_interval <- 20
#########################

# Configuration of the contours
contour_fig <- list(
    coloring = 'heatmap',
    start = min, end = max, size = interval,
    showlabels = TRUE, labelfont = list(size = 15, color = 'white'),
    showlines = FALSE))
    
# Configuration of the contour lines
line_fig <- list(smoothing = 1.3, width=0)

# Contour plot
DO <- plot_ly(combo,
        x = ~Date,
        y = ~Depth_m,
        z = ~Average.Dissolved_Oxygen_Sat,
        zmin = min,
        zmax = max,
        type = "contour",
        autocontour = TRUE,
        ncontours = numcontours,
        contours = contour_fig,
        colorscale = 'Blues',
        line = line_fig,
        colorbar = list(title = list(text = "DO\n(% Sat)", font = list(size = 15, color = 'black')),
                      tickfont = list(size=15, color='black'), tick0 = min, dtick = legend_interval,
                      ticks = "outside", ticklen = 5, tickcolor = 'black'),
        autocontour = TRUE,
        width = 1000, height = 500) %>%
        layout(yaxis = y_config, xaxis = x_config, legend = leg_config)

# Save the file
filename <- paste(folder_c, "/contour_plots_DO_", Location, ".pdf", sep="")
orca(DO, filename, width = 4 * 300, height = 2 * 300)

### pH ###

## Change these values ##
min <- 4
max <- 8
interval <- 0.1
numcontours <- 15
legend_interval <- 1
#########################

# Configuration of the contours
contour_fig <- list(
    coloring = 'heatmap',
    start = min, end = max, size = interval,
    showlabels = TRUE, labelfont = list(size = 15, color = 'white'),
    showlines = FALSE))
    
# Configuration of the contour lines
line_fig <- list(smoothing = 1.3, width=0)

# Contour plot
PH <- plot_ly(combo,
        x = ~Date,
        y = ~Depth_m,
        z = ~Average.pH,
        zmin = min,
        zmax = max,
        type = "contour",
        autocontour = TRUE,
        ncontours = numcontours,
        contours = contour_fig,
        colorscale = 'YlGnBu',
        line = line_fig,
        colorbar = list(title = list(text = "pH", font = list(size = 15, color = 'black')),
                      tickfont = list(size=15, color='black'), tick0 = min, dtick = legend_interval,
                      ticks = "outside", ticklen = 5, tickcolor = 'black'),
        autocontour = TRUE,
        width = 1000, height = 500) %>%
        layout(yaxis = y_config, xaxis = x_config, legend = leg_config)

# Save the file
filename <- paste(folder_c, "/contour_plots_pH_", Location, ".pdf", sep="")
orca(PH, filename, width = 4 * 300, height = 2 * 300)

### CONDUCTIVITY ###

## Change these values ##
min <- 0
max <- 100
interval <- 1
numcontours <- 15
legend_interval <- 20
#########################

# Configuration of the contours
contour_fig <- list(
    coloring = 'heatmap',
    start = min, end = max, size = interval,
    showlabels = TRUE, labelfont = list(size = 15, color = 'white'),
    showlines = FALSE))
    
# Configuration of the contour lines
line_fig <- list(smoothing = 1.3, width=0)

# Configuration of the colorscale
palette <- colorRampPalette(c("#ec9d44", "#462a54", "#07000a"))

# Contour plot
COND <- plot_ly(combo,
        x = ~Date,
        y = ~Depth_m,
        z = ~Average.Conductivity_µS_cm,
        zmin = min,
        zmax = max,
        type = "contour",
        autocontour = TRUE,
        ncontours = numcontours,
        contours = contour_fig,
        colorscale = 'Greens',
        #colorscale = cbind(seq(0, 1, by=1/(length(combo) - 1)), palette(length(combo))),
        line = line_fig,
        colorbar = list(title = list(text = "Conductivity\n(µS/cm)", font = list(size = 15, color = 'black')),
                      tickfont = list(size=15, color='black'), tick0 = min, dtick = legend_interval,
                      ticks = "outside", ticklen = 5, tickcolor = 'black'),
        autocontour = TRUE,
        width = 1000, height = 500) %>%
        layout(yaxis = y_config, xaxis = x_config, legend = leg_config)

# Save the file
filename <- paste(folder_c, "/contour_plots_Cond_", Location, ".pdf", sep="")
orca(COND, filename, width = 4 * 300, height = 2 * 300)
