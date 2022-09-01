#
# DiffBind.R
#
#  running comparison analyses & plots for pairs of sample peaks sets
#
# ####################################################################
#
#      CopyLeft (C) 2022 - Computational Genomics Lab @ UB
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
# ####################################################################
#

library(DiffBind);

print(paste0("DiffBind version: ",
             packageVersion("DiffBind")));

### Read Arguments from command-line
options(width=150,
        digits=3);

Args <- commandArgs(TRUE);

inputDir        <- Args[1]; # "$ATCDR/DiffBind_NF/anterior_atac_chip"
#setwd(inputDir)
outputDir       <- Args[2]; # "$BDIR/images/DiffBind_NF"
sampleTable     <- paste0(inputDir,"/",Args[3]); # "ATAC12h_anterior_atac_chip_vs_notum_Table_macs2_newG_NF.csv"
outputPrefix    <- Args[4]; # "anterior_atac_chip_notum"
maskA           <- Args[5]; # "anterior"
maskB           <- Args[6]; # "notum"
thy.summits     <- as.numeric(Args[7]); # 50
                   
       
outputBedFile   <- paste0(outputPrefix,"_consensus_peaks");
outDatabaseFile <- paste0(outputPrefix,"_dataTable");

    
### Variables
    
overLapCount <- 2;
fdrThreshold <- 0.050;
usePvalue    <- FALSE;
diffMethod   <- DBA_EDGER;
useTagwise   <- FALSE;

####################

# print out the working directory.
print(getwd());

# print the sample table.
print(sampleTable);

## Load the sample sheet and print out the sample sheet names.
samples <- read.csv(sampleTable);
names(samples);
samples;

## Read in the peakset using diffbind function and print the DBA object.
# This shows how many peaks are in each peakset, as well as
# (in the first line) the total number of unique peaks
# after merging and overlapping.
experiment <- dba(sampleSheet=sampleTable);
experiment;

## Occupancy Heat Map
# a correlation heatmap can be generated which gives
# an initial clustering of the samples using the cross-correlations
# of each row of the binding matrix
# Correlation heatmap, using occupancy (peak caller score) data.
png(paste0(outputDir,"/",outputPrefix,".00_occupancy_correlation_heatmap.png"),
    width = 10, height = 8, units = "in", res = 300);
plot(experiment, density = "none", keysize=1);
dev.off();

png(paste0(outputDir,"/",outputPrefix,".01_occupancy_correlation_heatmap.png"),
    width = 10, height = 8, units = "in", res = 300);
testing <- dba.plotHeatmap(experiment, density = "none", keysize=1);
dev.off();

# One reason to do an occupancy-based analysis is
# to determine what candidate sites should be used
# in a subsequent affinity-based analysis 
olap.rate <- dba.overlap(experiment,mode=DBA_OLAP_RATE);
olap.rate;
# The returned data in olap.rate is a vector containing the number of peaks
# that appear in at least one, two, three, and so on up to all eleven peaksets.
# These values can be plotted to show the overlap rate drop-off curve:

png(paste0(outputDir,"/",outputPrefix,".02_overlap_peaksets.png"),
    width = 10, height = 8, units = "in", res = 300);
plot(olap.rate,type='b',ylab='# peaks', xlab='Overlap at least this many peaksets');
dev.off();

# 
png(paste0(outputDir,"/",outputPrefix,".03_venn_sets_",maskA,".png"),
    width = 8, height = 8, units = "in", res = 300);
dba.plotVenn(experiment, experiment$masks[[maskA]],
             main=paste0("Overlaps between replicates: ",maskA));
dev.off();

png(paste0(outputDir,"/",outputPrefix,".04_venn_sets_",maskB,".png"),
    width = 8, height = 8, units = "in", res = 300);
dba.plotVenn(experiment, experiment$masks[[maskB]],
             main=paste0("Overlaps between replicates: ",maskB));
dev.off();

png(paste0(outputDir,"/",outputPrefix,".05_venn_sets_",maskA,"-x-",maskB,".png"),
    width = 8, height = 8, units = "in", res = 300);
dba.plotVenn(experiment, experiment$masks[[maskA]] | experiment$masks[[maskB]],
             main=paste0("Overlaps between replicates: ",maskA," x ",maskB));
dev.off();

    
# Next step is to calculate a binding matrix with scores based on read counts
# for every sample (affinity scores), rather than confidence scores for only
# those peaks called in a specific sample (occupancy scores).
# Re-center each peak around the point of greatest enrichment with summits=250,
# the peaks will be 500bp, extending 250bp up- and down- stream of the summit

experiment <- dba.count(experiment, summits=thy.summits,
                        score=DBA_SCORE_RPKM, minOverlap=overLapCount);

png(paste0(outputDir,"/",outputPrefix,".06_occupancy_correlation_heatmap.png"),
    width = 10, height = 8, units = "in", res = 300);
dba.plotHeatmap(experiment, density = "none", keysize=1);
dev.off();

#print out the dba object
experiment;

## Establishing a contrast
# Tell diffbind what to compare
# Then uses the Condition metadata (Responsive vs. Resistant) to set up a a contrast
experiment <- dba.contrast(experiment, categories=DBA_CONDITION, minMembers=2);

## Performing the differential analysis
# The main differential analysis function is invoked as follows:
experiment <- dba.analyze(experiment, method=diffMethod, bTagwise=useTagwise);
experiment;

# Correlation heatmap, using only significantly differentially bound sites.
png(paste0(outputDir,"/",outputPrefix,".07_occupancy_correlation_heatmap.significant_diff_bind_sites.png"),
    width = 10, height = 8, units = "in", res = 300);
plot(experiment, contrast=1, th=fdrThreshold, bUsePval=usePvalue, method=diffMethod);
dev.off()

# Retrieving the differentially bound sites
cur.wdir <- getwd();
setwd(outputDir);
experiment.DB <- dba.report(experiment, bCounts = TRUE,
                            file = outDatabaseFile, bCalled = TRUE, th=1,
                            bCalledDetail = TRUE, method=diffMethod,
                            bNormalized = FALSE)
setwd(cur.wdir);
experiment.DB;

# A PCA plot,  using affinity data for all sites,
# which includes normalized read counts
# for all the binding sites, can be obtained as follows:
png(paste0(outputDir,"/",outputPrefix,".08_PCA_all_bind_sites.png"),
    width = 10, height = 8, units = "in", res = 300);
dba.plotPCA(experiment, DBA_TISSUE, label=DBA_CONDITION);
dev.off();

# A PCA plot using only the differentially bound sites (corresponding to Figure 3),
# using an FDR threshold of 0.05, can be drawn as follows:
png(paste0(outputDir,"/",outputPrefix,".09_PCA_diff_bind_sites.png"),
    width = 10, height = 8, units = "in", res = 300);
dba.plotPCA(experiment, contrast=1, label=DBA_TISSUE,
            th=fdrThreshold, bUsePval=usePvalue, method=diffMethod);
dev.off();

# If you want to see where the replicates
# for each of the unique cell lines lies then
png(paste0(outputDir,"/",outputPrefix,".10_PCA_all_bind_sites_2factor.png"),
    width = 10, height = 8, units = "in", res = 300);
dba.plotPCA(experiment, attributes=c(DBA_TISSUE,DBA_CONDITION),
            label=DBA_REPLICATE, method=diffMethod);
dev.off();

# Seeing the first three principal components
# can be a useful exploratory exercise
png(paste0(outputDir,"/",outputPrefix,".11_PCA_all_bs_3pca.png"),
    width = 10, height = 8, units = "in", res = 300);
dba.plotPCA(experiment, b3D=T);
dev.off();

# On next plot each point represents a binding site,
# with points in red representing sites identified as differentially bound. 
# MA plot of Resistant-Responsive contrast.
# Sites identified as significantly differentially bound shown in red.
png(paste0(outputDir,"/",outputPrefix,".12_MAplot_binding_affinity.png"),
    width = 8, height = 8, units = "in", res = 300);
dba.plotMA(experiment, th=fdrThreshold,
           bUsePval=usePvalue, method=diffMethod);
dev.off();

# This same data can also be shown with the concentrations
# of each sample groups plotted against each other
png(paste0(outputDir,"/",outputPrefix,".13_MAplot_concentration.png"),
    width = 8, height = 8, units = "in", res = 300);
dba.plotMA(experiment, bXY=TRUE, th=fdrThreshold,
           bUsePval=usePvalue, method=diffMethod);
dev.off();

###
### NO dba.plotVolcano for DiffBind version < 2
###
# # Volcano plot of Resistant-Responsive contrast.
# # Sites identified as significantly differentially bound shown in red. 
# png(paste0(outputDir,"/",outputPrefix,".14_volcanoplot.png"),
#     width = 10, height = 8, units = "in", res = 300);
# dba.plotVolcano(experiment, th=fdrThreshold,
#                 bUsePval=usePvalue, method=diffMethod);
# dev.off();

# Box plots of read distributions for significantly
# differentially bound (DB) sites. Tamoxifen resistant
# samples are shown in red, and responsive samples are shown in blue.
# Left two boxes show distribution of reads over all
# DB sites in the Resistant and Responsive groups;
# middle two boxes show distributions of reads in DB sites that increase
# in affinity in the Responsive group; last two boxes show distributions
# of reads in DB sites that increase in affinity in the Resistant group. 
png(paste0(outputDir,"/",outputPrefix,".15_boxplot.png"),
    width = 10, height = 8, units = "in", res = 300);
pvals <- dba.plotBox(experiment, th=fdrThreshold,
                     bUsePval=usePvalue, method=diffMethod,
                     notch=FALSE);
dev.off();

# dba.plotBox returns a matrix of p-values (computed using a two-sided Wilcoxon 'Mann-Whitney' test,
# paired where appropriate, indicating which of these distributions are significantly
# different from another distribution.
pvals;

## Heatmaps
# Binding affinity heatmap showing affinities for differentially bound sites.
# Samples cluster first by whether they are responsive to tamoxifen treatment, then by cell line.
# Clusters of binding sites show distinct patterns of affinity levels.
png(paste0(outputDir,"/",outputPrefix,".16_heatmap_diff_bound_sites.png"),
    width = 10, height = 8, units = "in", res = 300);
corvals <- dba.plotHeatmap(experiment, contrast=1,
                           correlations=FALSE, maxSites=100000,
                           density = "none", keysize=1,
                           th=fdrThreshold, bUsePval=usePvalue,
                           method=diffMethod);
dev.off();

proc.time();

cat("# DONE...\n");
