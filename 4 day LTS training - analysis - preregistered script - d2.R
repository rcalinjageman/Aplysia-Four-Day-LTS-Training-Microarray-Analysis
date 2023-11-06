# Planned analysis of 4-day LTS training microarray 
# v1.1  - Bob Calin-Jageman, 7/1/2022 - Preregistered Script
# v1.2  - Bob Calin-Jageman, 11/5/2023 - Final params for real analysis

# Things to be updated after data collections:
# * We may update the microarry platform based on long-read sequencing.
#   if so, the data for comparing to the 1-day training protocol (d1rep)
#   will need to be re-keyed to the updated array, or if that is not possible
#   that comparison condition will be dropped.  It is only serving as a reference
#   to help anchor evaluation of similarity between the 1 day and 5 day timepoints
#
#   - No updates were made to the microarray platform
#
#
# * Funding may not permit running the 11-day microarray.  If so, 
#   all analyses with the 11-day data will be deleted
#
#   - 11-day was completed; adjustments made below for those analyses
#
# * Normexp offset has to be adjusted manually after seeing the MA plots
# * Some graph parameters will need to be adjusted (e.g. x/y limits) based on
#   data obtained
#
#   - No updates needed

#Notes
# For final data analysis --------------------------------------------------
# This script contains minor (documented) updates to our pre-registered
#  microarray analysis plan to enable it to run with the real data.
#  To document all changes, the script and analysis files were published to
#  github, and all changes committed/explained:
#  https://github.com/rcalinjageman/Aplysia-Four-Day-LTS-Training-Microarray-Analysis
#
# * Updated notes
# * Updated params for real analysis, Targets.txt for slide names, and
#    sample size for 1-day due to one slide being omitted for a processing error.
# * Updated call to falseNegativeCount (wrong arguments were being passed) and
#   fixed flaw in that function that allowed negative proportions of false
#   negatives to be returned.
# * Fixed error in build of table 1
# * Added planned interaction test from d1 to d11 (planned but had not been
#   included in original script due to uncertainty if funding would permit
#   d11 to be run).
# * Added some additional scripting to pull together outputs into one file.  
#   And fixed output of MA2avg to csv and ordered targets by condition.
#
#
# Notes from pre-registered analysis script --------------------------------
# This is a pre-planned analysis script for a project (https://osf.io/wvx6z/) 
#  examining the transcriptional changes that accompany maintenance of a very
#  long lasting form of sensitization in Aplysia.
# We expect wide-spread transcriptional regulation 1-day after training
#  The key question is if this transcriptional response is largely maintained
#  as maintenance progresses, or if the memory trace can become transcriptionally
#  neutral.
# This r script analyses differential expression for trained vs. untrained 
#   samples of pleural ganglia from Aplysia from 3 different conditions:
#    1) 1 day after extended 4-day LTS training
#    2) 5 days after extended 4-day LTS training
#    3) (If funding permits) 11 days after extended 4-day LTS training
#  Extended 4-day LTS training will consist of 4 consecutive days of training.
#   Each day will have 4 rounds of shocks, 30 minutes between each round.
#   Each round will be 90mA of 0.5s on 0.5s off for 10s applied to one side of
#   the body
#  The goal is 8 biological replicates were harvested for each condition.  
#  A two-color approach (trained vs. control) will be used on a 
#   custom 8x60k Agilent microarray.  
#  Dye color will be counterbalanced across conditions.  
#  Samples from each condition were also stratefied across slides.

# This script first completes basic mcirarraoy processing:
# It reads a design file, reads the agilent file, adjusts for background signal 
#   using normexp algorithm, normalizes for dye intensity within each array, 
#   filters out control genes, then averages by probeset (for transcripts measured
#   multiple probes.

# Next, this script tests a number of research questions about regulation of 
#   expression in this experiment; see comments below

# Development of script--was based primarily on the Limma user guide and 
#  this guide for two-color Agilent analysis:
#  http://matticklab.com/index.php?title=Two_channel_analysis_of_Agilent_microarray_data_with_Limma
#  with some additional help from the references cited below.

# Decisions:
#  Median or mean signal?  Zahurak et al. (2007) report the choice makes no 
#    practical difference in their analysis of Agilent wo-color arrays.  
#    Therefore, default of median was selected.
#  Background subtraction?  Zahurak et al. (2007) report best results with 
#    no correction.  However, in a more extensive comparison.
#	   Ritchie et al. (2007) report best results with normexp+offset for background 
#    correction.  
#	Normexp offset?  Offset needs to be empirically determined by inspection of 
#     MA plots, but analysis is not super-sensitive to choice of offset value (k).  
#     Ritchie et al. suggest an offset of 50 based on their previous work, 
#     but also suggest inspection of MA plots to check other values.  
#     This will be done after results are obtained.  
#	Normalize Within Arrays?  Zahurak et al. (2007) found best results for 
#     Agilent two-color with loess normalization.  This also 
#		  seems to have become the standard within-array normalization strategy.
# How to assess practical significance?  Most statistical testing is done 
#   against an LFC of 0, which leads to some results which are significant 
#   but not meaningful in terms of the magnitude of regulation.  
#   This can be fixed by adding an LFC threshold after significance testing.  
#   Alternatively, McCarthy & Smith, 2009 recommend conducting statistical testing
#		against a more reasonable null hypothesis, thus testing for practical and 
#   statistical significance at the same time.  The approach of checking 
#   against a more rigorous null was used.  Specifically, an FC of 1.1 was selected.

# -----------------
# Setup - load needed libraries.  Be sure these are installed before running the script
#   Note that to install limma and OrdererList one must install bioconductor and install via bioconductor:
#   install.packages("BiocManager")
#   BiocManager::install("limma")
#   BiocManager::install("OrderedList")
# -----------------
  library(limma)	
  library(ellipse)
  library(ggplot2)
  library(viridis)
  library(OrderedList)
  library(psych)
  source("estimateProportions.R")


# -----------------
# Helper functions
# -----------------
  falseNegativeCount <- function(resultsTable) {
    # This function returns the estimated proportion and count of false positives by analyzing the distribution of p values in a resultstable.
    # First, we get the total number of rows in the results table
    # Then we use propTrueNull to obtain an estimate of the true null effects, given the distribution of p values.  
    #   This uses the convest approach described by Langaas et al., 2005 which is the most sophisticated approach available
    #   Note that uncrorrected p values are used--corrected p values are not suitable for this type of analysis
    # We also obtain the positives (true and false) by tallying all transcripts significantly regulated (with adjustment) and dividing by total
    # Finally, we estimate false negative rate as 1 - proportion true nulls - proportion positives
    # The total true positives is then calculated as proportion true positive * total 
    # Get total number of transcripts
    total <- nrow(resultsTable)
    # Get estimated pTrueNull
    pTrueNull <- propTrueNull(resultsTable$P.Value, method="convest")
    # Estimate total positives
    phit <- nrow(resultsTable[resultsTable$adj.P.Val < .05, ]) / total
    # Now estimate proportion falsenegatives
    fn <- 1 - pTrueNull - phit
    # Updated to avoid negative estimates of false negatives
    if (fn < 0) fn <- 0
    # Finally, calculate total number of possible misses
    fnTotal <- fn * total
    # Make a list to return results
    res <- list(
      fn = fn,
      fnTotal = fnTotal
    )
    return(res)
  }


# -----------------
# Set parameters for this analysis
# -----------------
  # The runtype variable enables this script to switch between sample/test and real runs
  # Set runtpe to be Test or Real 
  runtype <- "Real"
  #runtype <- "Test"

  # runparams stores the different parameters for the different types of runs
  # Parameters for running the script
  runparams <- list(
    Test= list (
      prefix = "Test",
      targets = "Targets - test.txt",
      arraydir = "sample_arraydata",
      nCondition <- c(8,8,8,8)
    ),
    Real = list (
      prefix = "Real",
      targets = "Targets.txt",
      arraydir = "arraydata",
      nCondition <- c(7,8,8,8)
    )
  )
  nCondition <- runparams[[runtype]][[4]]
  
  # TO BE UPDATED - we may not run the 11-day, in which case the
  #  params for the 'real' run will need to be adjusted
  
# -----------------
# Read targets and data
# -----------------
  # Read targets, load sample sizes, and read raw data
  targets <- readTargets(runparams[[runtype]]$targets, row.names = "Name")
  targets <- targets[order(targets$Name), ]
  RG <- read.maimages(targets, source="agilent.median", path=runparams[[runtype]]$arraydir)

  # Now load and read the targets from 1 day after a shorter LTS protocol.
  #  This will be used as a benchmark that should represent fairly high similariy
  # TO BE UPDATED - Targets - replication will be updated to contain the d1 files and the prior d1 files 
  # rep_targets <- readTargets("Targets - replication.txt", row.names = "Name")
  # repRG <- read.maimages(rep_targets, source="agilent.median", path=runparams[[runtype]]$arraydir)
  
# -----------------
# Pre-processing and quality check outputs
# -----------------
  # Plot MAs of non-adjusted data
  plotMA3by2(RG, path="output", main=targets$Name, status=RG$genes$ControlType, prefix=paste(runparams[[runtype]]$prefix, "prenormMA", sep="-"))
  # Background correct
  RG <- backgroundCorrect(RG, method="normexp", offset=50)
    # TO BE UPDATED - using default offset of 50 but MA plots will be inspected to check other values
  # Normalize
  MA <- normalizeWithinArrays(RG, method="loess")
  # Output the table for GEO - dye assignment is not corrected, so M values will be inverted on dye-swap samples
  write.table(MA, file=paste("./output/", runparams[[runtype]]$prefix, "- GeoTable.txt", sep=""))
  # Plot MAs after normalization
  plotMA3by2(MA, path="output", main=targets$Name, status=MA$genes$ControlType, prefix=paste(runparams[[runtype]]$prefix, "postnormMA", sep="-"))
  # Remove control genes
  MA2 <- MA[MA$genes$ControlType==0,]
  # Average across probes designed to the same transcript
      if(is.null(MA2$genes$GeneName)) {
        MA2$genes$GeneName <- MA2$genes$SystematicName
      }
  MA2avg <- avereps(MA2, ID=MA2$genes$GeneName)
  # Save MA2 table - again, dye assignment is not corrected
  write.table(MA2avg, file=paste("./output/", runparams[[runtype]]$prefix, " - MA2avg.csv", sep = ""), sep = ",")
  # Plot MDS of all samples
    # Note we use dye swaps so samples will cluster not only by condition but also by dye assignment
  plotMDS(MA2avg, labels = targets$Name)

  #Same as above but for the replication benchmark.  No need to output MA plots, though
  # repRG <- backgroundCorrect(repRG, method = "normexp", offset = 50)
  # repMA <- normalizeWithinArrays(repRG, method = "loess")
  # repMA2 <- repMA[repMA$genes$ControlType==0,]
  #     if(is.null(repMA2$genes$GeneName)) {
  #       repMA2$genes$GeneName <- repMA2$genes$SystematicName
  #     }
  # repMA2avg <- avereps(repMA2, ID=repMA2$genes$GeneName)
  
    
# -----------------
# Analysis - Design
# -----------------
  # Create a basic design matrix comparing each condition to control
  design <- modelMatrix(targets, ref="Control")
  # Add a dye factor to estimate dye effects
  design <- cbind(Dye = 1, design)
  # This command adds a plate factor as well. 
  design <- cbind(Plate = targets$SlideNumber, design)
  
  # Create design matrix for comparison with replication set
  # repdesign <- modelMatrix(rep_targets, ref = "Control")
  # repdesign <- cbind(Dye = 1, repdesign)
  # repdesign <- cbind(Plate = rep_targets$SampleNumber, repdesign)
  

# -----------------
# Basic Gene Lists
# -----------------
  # Fit the basic model and test against an interval null LFC: -1.1 to 1.1
  fit <- lmFit(MA2avg, design)
  fit <- treat(fit, lfc=log2(1.1), trend=TRUE)

  # Same for comparing d1 to d1 replication  
  # repfit <- lmFit(repMA2avg, repdesign)
  # repfit <- treat(repfit, lfc=log2(1.1), trend=TRUE)
  
  # Create a results table for each condition
  resultsTable <- list(
    d1 = topTreat(fit, coef="d1", number=nrow(fit), adjust.method="BH"),
    d5 = topTreat(fit, coef="d5", number=nrow(fit), adjust.method="BH"),
    d11 = topTreat(fit, coef="d11", number=nrow(fit), adjust.method="BH"),
    d1_rep = topTreat(fit, coef="d1rep", number=nrow(fit), adjust.method="BH")
  )
  # Write results tables
  for (i in 1:4) {
    write.csv(resultsTable[[i]], file = paste("./output/", runparams[[runtype]]$prefix, " - resultsTable - ", names(resultsTable[i]), ".csv", sep=""))
  }
  
  
  # Same as above, but for each condition filter down to transcripts that are significantly > LFC 1.1
  geneList <- list(
    d1 = resultsTable$d1[resultsTable$d1$adj.P.Val < .05, ],
    d5 = resultsTable$w1[resultsTable$w1$adj.P.Val < .05, ],
    d11 = resultsTable$sav[resultsTable$sav$adj.P.Val < .05, ],
    d1_rep = resultsTable$d1_rep[resultsTable$d1_rep$adj.P.Val < .05, ]
  )

  # Create a fit object that is restricted just to what is regulated at 1d - we'll use this in later analyses
  fit_1donly <- fit[fit$genes$GeneName %in% geneList$d1$GeneName, ]
  
  # Diagnostic - How complete are the gene lists?
  # For each condition, estimate the false negative (misses) proportion and counts
    falseNegatives <- list(
      d1 = falseNegativeCount(resultsTable$d1),
      d5 = falseNegativeCount(resultsTable$d5),
      d11 = falseNegativeCount(resultsTable$d11),
      d1_rep = falseNegativeCount(resultsTable$d1_rep)
    )

  
  # Table 1- summarize gene list and estimated false negatives
    # Summarize regulation by condition
    # To do this, filter down to just yes/no, BH adjustment on each condition - same as the resultstable above, 
    #    but just yes/no calls and all 3 conditions together
    dtests <- decideTests(fit, adjust.method = "BH", method = "separate")
    # Print a summary, showing up, down, and no call for each condition, including dye and plate effects
    rtable <- summary(dtests)
    # Create data frame
    table1 <- data.frame(Group = c("1day", "5day", "11day", "1day_weaker_training"))
    # Store results from the decide test... not very elegant, but rtable comes out as some type of table and as.data.frame doesn't work well on it
    table1$up_regulated <- c(rtable["Up", "d1"], rtable["Up", "d5"], rtable["Up", "d11"], rtable["Up", "d1rep"])
    table1$down_regulated <- c(rtable["Down", "d1"], rtable["Down", "d5"], rtable["Down","d11"], rtable["Down", "d1rep"])
    table1$false_negative_rate <- c(falseNegatives$d1$fn, falseNegatives$d5$fn, falseNegatives$d11$fn, falseNegatives$d1_rep$fn)
    table1$false_negative_count <- c(falseNegatives$d1$fnTotal, falseNegatives$d5$fnTotal, falseNegatives$d11$fnTotal, falseNegatives$d1_rep$fnTotal)
    table1
    write.csv(table1, file = paste("./output/", runparams[[runtype]]$prefix, " - Table 1 - gene list counts.csv", sep = ""))  
  

# -----------------
# Research question: How does overlap from d5 compare to overlap from d1?
# -----------------
  # Start with a venn diagram showing overlap of regulation
    # Save the decide tests again; helpful to have this here for debugging
    dtests <- decideTests(fit, adjust.method = "BH", method = "separate")
    # Drop the plate and dye effects
    dtests <- dtests[ , colnames(dtests) != "Plate"]
    dtests <- dtests[ , colnames(dtests) != "Dye"]
    # Also drop d1_rep to make venn diagram more readable
    dtests <- dtests[ , colnames(dtests) != "d1rep"]
    # Make a simplified vennDiagram showing just the up/down in the d1, d5, and d11
    jpeg(paste("./output/", runparams[[runtype]]$prefix, " - figure_venn.jpg", sep = ""))
    vennDiagram(dtests, counts.col = c("black", "black", "black", "black", "black", "black", "red"))
    dev.off()

  # Now make a direct contrast between d5 and d11 among transcripts that had been regulated at d1
    # This is a direct test for an interaction, checking for transcripts regulated significantly differently at d5 and d11
    # Make a contrast matrix for comparing d1 and d5
    contrast.matrix <- makeContrasts(d5vd1=d5-d1, levels=c("Plate", "Dye", "d1", "d11", "d1rep","d5"))
    # Now analyze using contrasts.fit and the fit_1donly object which has only the transcripts regulated at 1d
    int_fit <- contrasts.fit(fit_1donly, contrast.matrix)
    int_fit <- treat(int_fit, lfc=log2(1.1), trend=TRUE)
    int_list <- topTreat(int_fit, coef="d5vd1", number=nrow(int_fit), adjust.method="BH")
    write.csv(int_list, file = paste("./output/", runparams[[runtype]]$prefix, " - Table - d5_d1 interaction gene list.csv", sep = ""))
    int_regulated <- nrow(int_list[int_list$adj.P.Val < .05, ])
    print(paste("Number of transcripts at d5 regulated differently than at d1: ", int_regulated))
    

    # Make a contrast matrix for comparing d1 and d5
    contrast.matrix.d11 <- makeContrasts(d11vd1=d11-d1, levels=c("Plate", "Dye", "d1", "d11", "d1rep","d5"))
    # Now analyze using contrasts.fit and the fit_1donly object which has only the transcripts regulated at 1d
    int_fit_d11 <- contrasts.fit(fit_1donly, contrast.matrix.d11)
    int_fit_d11 <- treat(int_fit_d11, lfc=log2(1.1), trend=TRUE)
    int_list_d11 <- topTreat(int_fit_d11, coef="d11vd1", number=nrow(int_fit), adjust.method="BH")
    write.csv(int_list_d11, file = paste("./output/", runparams[[runtype]]$prefix, " - Table - d11_d1 interaction gene list.csv", sep = ""))
    int_regulated_d11 <- nrow(int_list_d11[int_list_d11$adj.P.Val < .05, ])
    print(paste("Number of transcripts at d11 regulated differently than at d1: ", int_regulated_d11))
    
        
    
    #Now, let's get the proportion of d1 and d5 overlap and compare these proportions
    # Store a temporary copies of the d1, d5, and d11 results
    qdata <- resultsTable$d1
    temp_d5 <- resultsTable$d5
    temp_d11 <- resultsTable$d11
    # For each result set, define a column "time" that represents days since training
    qdata$time = 1
    temp_d5$time = 5
    temp_d11$time = 11
    
    qdata$condition <- "Newly Maintained"
    temp_d5$condition <- "Moderate Maintenance"
    temp_d11$condition <- "Long Maintenance"
    
    # Ensure each data set is in the same order
    qdata <- qdata[order(qdata$GeneName), ]
    temp_d5 <- temp_d5[order(temp_d5$GeneName), ]
    temp_d11 <- temp_d11[order(temp_d11$GeneName), ]
    # Make a 'type' column that contcatenates significance status at each time point
    qdata$type <- factor(paste(qdata$adj.P.Val < .05, temp_d5$adj.P.Val < .05, temp_d11$adj.P.Val < .05))
    # Also, make an intype column that represents if the gene shows a significant interaction
    qdata$intype <- "Sig Different"
    qdata[qdata$GeneName %in% int_list[int_list$adj.P.Val > .05, ]$GeneName, ]$intype <- "Not significant"
    # Relabel the type column into logical groupings based on time-course of regulation
    #  e.g. encoding_only is regulated at d1 but not d5 or d11
    #  Of critical interest will be savings_reactivated--transcripts regulated at d1, not at w1, but again at savings
    levels(qdata$type)[levels(qdata$type) == "FALSE FALSE FALSE"] <- "unregulated"
    levels(qdata$type)[levels(qdata$type) == "TRUE FALSE FALSE"] <- "Early only"
    levels(qdata$type)[levels(qdata$type) == "FALSE TRUE FALSE"] <- "Mid-term only"
    levels(qdata$type)[levels(qdata$type) == "FALSE FALSE TRUE"] <- "Late only"
    levels(qdata$type)[levels(qdata$type) == "TRUE TRUE TRUE"] <- "Always_regulated"
    levels(qdata$type)[levels(qdata$type) == "TRUE FALSE TRUE"] <- "Up-Down-Up"
    levels(qdata$type)[levels(qdata$type) == "FALSE TRUE TRUE"] <- "Late-developing"
    levels(qdata$type)[levels(qdata$type) == "TRUE TRUE FALSE"] <- "Early and Mid"
    # Store type variable in w1 and sav, requires same order, which was set above
    temp_d5$type <- qdata$type
    temp_d11$type <- qdata$type
    temp_d5$intype <- qdata$intype
    temp_d11$intype <- qdata$intype
    
    #get overlaps counts from d1 to each group
    d5_overlap_count <- nrow(qdata[qdata$type == "Always_regulated", ] ) + nrow(qdata[qdata$type == "Early and Mid", ] )
    d11_overlap_count <- nrow(qdata[qdata$type == "Always_regulated", ] ) + nrow(qdata[qdata$type == "Up-Down-Up", ] )
    d1rep_overlap_count <- nrow(geneList$d1_rep[geneList$d1_rep$GeneName %in% geneList$d1$GeneName, ])
    d1_total <- nrow (geneList$d1)
    #Now estimate and print the difference in proportion of overlap
    overlap_diff <- estimateProportionDifference.numeric(d5_overlap_count, d1_total, d11_overlap_count, d1_total, 
                                                         caselabels = c("Shared", "d1_only"), 
                                                         grouplabels = c("day5", "day11")
                                                         )
    d1rep_overall <- estimateProportion(d1rep_overlap_count, d1_total, caselabels = c("Shared", "d1_only"))
    d1rep_overall$summary_data <- cbind(data.frame("Condition" = "D1_replication"), d1rep_overall$summary_data)
    overlap_diff$summary_data <- rbind(overlap_diff$summary_data, d1rep_overall$summary_data)
    overlap_diff$summary_data
    write.csv(overlap_diff$summary_data, file = paste("./output/", runparams[[runtype]]$prefix, " - Table - overlap change.csv", sep = ""))
  
        
  # Now make a figure showing the time-course of regulation for the genes regulated at encoding, colored by what happens to them
    # Combine all three time points together, long format for ggplot
    qdata <- rbind(qdata, temp_d5, temp_d11)
    intdata <- qdata
    # Select just the transcripts regulated at some time point
    qdata <- qdata[qdata$type != 'unregulated', ]

    # Now build the plot
    qplot <- ggplot() 
    qplot <- ggplot(data = qdata, aes(x = time, y = logFC, group = GeneName, color = type)) + geom_line(aes(alpha = type))
    # Adjust these based on number of regulated transcripts of each type
    # qplot <- qplot + scale_alpha_manual(values = c(0.05, 0.1, 1, 1)) 
    # qplot <- qplot + scale_color_manual(values = c("black", "red", "blue", "gold"))
    # qplot <- qplot + guides(colour = guide_legend(override.aes = list(alpha = 1)))
    qplot <- qplot + theme_classic()
    qplot <- qplot + geom_hline(yintercept = 0, linetype = "dashed")
    if(runtype != "Real") {
      qplot <- qplot + annotate("text", label="sample data", x = 3, y = 0, angle = 45, alpha = .1, size = 14)
    }
    qplot
    ggsave(plot = qplot, filename = paste("./output/", runparams[[runtype]]$prefix, " - Figure encoding fate.jpg", sep = ""))    


  # Create a figure showing quantity of up- and down-reguated transcripts at each time point
    # Initialize a data frame to store the summary data
    ldata <- data.frame(time = NA, type = NA, m = NA, n = NA)
    # Initialize a vector to represent the d1, w1, and savings conditions in terms of days from training
    timepoint <- c(1, 5, 11)
    # Cycle through the conditions
    for (i in 1:3) {
      # Count number of transcripts up-regulated
      nup = nrow(geneList[[i]][geneList[[i]]$logFC > 0, ])
      if (is.null(nup)) nup <- 0
      # Calculate mean of up-regulated transcripts, scoring a 0 if there are none
      if (nup == 0) {
        mup = 0
      } else {
        mup <- mean(geneList[[i]][geneList[[i]]$logFC > 0, ]$logFC)
      }
      # Now count and get mean for down-regulated transcripts
      ndown = nrow(geneList[[i]][geneList[[i]]$logFC < 0, ])
      if (is.null(ndown)) ndown <- 0
      if (ndown == 0) {
        mdown = 0
      } else {
        mdown <- mean(geneList[[i]][geneList[[i]]$logFC < 0, ]$logFC)
      }
      # store in the data frame
      ldata <- rbind(ldata, data.frame(time = timepoint[i], type = "Up", m = mup, n = nup))
      ldata <- rbind(ldata, data.frame(time = timepoint[i], type = "Down", m = -mdown, n = -ndown))
    }
    # delete the initial NA row
    ldata <- ldata[2:nrow(ldata), ]
    # Create a plot of up- and down-regulated hits in d1, w1, and savings
    lplot <- ggplot(data = ldata, aes(x = time, y = n, group = type, color = type)) + geom_line(size = 2)
    lplot <- lplot + scale_color_manual(values = c("blue", "red"))
    lplot <- lplot + geom_hline(yintercept = 0, linetype="dashed")
    lplot <- lplot + theme_classic()
    if(runtype != "Real") {
      lplot <- lplot + annotate("text", label="sample data", x = 3, y = 0, angle = 45, alpha = .1, size = 14)
    }
    
    lplot
    ggsave(plot = lplot, filename = paste("./output/", runparams[[runtype]]$prefix, " - Figure count by timeoint.jpg", sep = ""))    
    

# -----------------
# Research question: How similar are the gene lists, qualitatively?
# -----------------
    
  # Calculate rank similarity across gene lists
    # First, make a column that ranks genes from most upregulated to least regulated to most down-regulated
    # We do this by taking the log of the P.value (best evidence -> lowest p value -> large negative number)
    # Then we multiple this by the sign of the fold change 
    # This gives strongly up-regulated transcripts a large negative number, note regulated a 0 (p = 1, so log10(1) = 0), and down-regulated large positive numbers
    # And this sorts so that up-regulated at top, unregulated in the middle, and down-regulated at the bottom
    for (i in 1:4) {
      resultsTable[[i]]$sortit <- log10(resultsTable[[i]]$P.Value) * sign(resultsTable[[i]]$logFC)
    }

    # Now we do a rank comparison of the gene lists, first comparing d1 to w1, then comparing d1 to sav, and finally d1 to d1_rep
    # Note there is no point in comparing sav to w1 because there is not very much regulation at w1, so no much information it the sorted list
    # The question is if savings is dissimilar to d1, as predicted by re-covery, or highly similar to d1, as predicted by re-learning
    
    # Use the compare lists function, feeding in d1 and w1 in sorted order
    #   Rather than let the function tune alpha, we set it to .01 which corresponds to gene lists of ~1000 in length
    #   This holds the analysis constant for both comparisons
      # First, the comparLists function calculates the rank comparisons
      cl_d1d5 <- compareLists(resultsTable$d1[order(resultsTable$d1$sortit), ]$GeneName, resultsTable$d5[order(resultsTable$d5$sortit), ]$GeneName, alphas = c(0.01))
      # Then, the getOverlap function gives an overlap score and p value
      ol_d1d5 <- getOverlap(cl_d1d5)
      # Print both objects
      cl_d1d5
      ol_d1d5
      # Make a graph showing observed vs. expected rank overlaps by rank length
      plot(ol_d1d5)
      jpeg(paste("./output/", runparams[[runtype]]$prefix, " - overlap_d1d5.jpg", sep = ""))
      plot(ol_d1d5)
      dev.off()
      # Plot observed similarity scores against histogram of similarity scores from random shuffles 
      plot(ol_d1d5, "scores")
      jpeg(paste("./output/", runparams[[runtype]]$prefix, " - similarity_d1d5.jpg", sep = ""))
      plot(ol_d1d5, "scores")
      dev.off()
      # Plot a scatter plot of the expression of the overlapping genes in the gene list
      plot(x = fit[fit$genes$GeneName %in% ol_d1d5$intersect, ]$coefficients[, "d1"], fit[fit$genes$GeneName %in% ol_d1d5$intersect, ]$coefficients[, "d5"])
      
      
      # Repeat the same steps as above, but comparing d1 to sav
      cl_d1d11 <- compareLists(resultsTable$d1[order(resultsTable$d1$sortit), ]$GeneName, resultsTable$d11[order(resultsTable$d11$sortit), ]$GeneName, alphas = c(0.01))
      ol_d1d11 <- getOverlap(cl_d1d11)
      cl_d1d11
      ol_d1d11
      plot(ol_d1d11)
      jpeg(paste("./output/", runparams[[runtype]]$prefix, " - overlap_d1d11.jpg", sep = ""))
      plot(ol_d1d11)
      dev.off()
      plot(ol_d1d11, "scores")
      jpeg(paste("./output/", runparams[[runtype]]$prefix, " - similarity_d1d11.jpg", sep = ""))
      plot(ol_d1d11, "scores")
      dev.off()
      plot(x = fit[fit$genes$GeneName %in% ol_d1d11$intersect, ]$coefficients[, "d1"], fit[fit$genes$GeneName %in% ol_d1d11$intersect, ]$coefficients[, "d11"])
    
      # And again comparing d1 to d1rep
      cl_d1d1rep <- compareLists(resultsTable$d1[order(resultsTable$d1$sortit), ]$GeneName, resultsTable$d1_rep[order(resultsTable$d1_rep$sortit), ]$GeneName, alphas = c(0.01))
      ol_d1d1rep <- getOverlap(cl_d1d1rep)
      cl_d1d1rep
      ol_d1d1rep
      plot(ol_d1d1rep)
      jpeg(paste("./output/", runparams[[runtype]]$prefix, " - overlap_d1d1rep.jpg", sep = ""))
      plot(ol_d1d1rep)
      dev.off()
      plot(ol_d1d1rep, "scores")
      jpeg(paste("./output/", runparams[[runtype]]$prefix, " - similarity_d1d1rep.jpg", sep = ""))
      plot(ol_d1d1rep, "scores")
      dev.off()
      # This seems to work, but could require making sure both fit objects are sorted int the same order first
      plot(x = fit[fit$genes$GeneName %in% ol_d1d1rep$intersect, ]$coefficients[, "d1"], y = fit[fit$genes$GeneName %in% ol_d1d1rep$intersect, ]$coefficients[, "d1rep"])
      
      # Now assemble and save output
      table_ro <- data.frame("Measure" = c("Weighted overlap score to d1 Group", "p value"))
      table_ro$to_d5 <- c(ol_d1d5$score, ol_d1d5$pvalue)
      table_ro$to_new_replication <- c(ol_d1d1rep$score, ol_d1d1rep$pvalue)
      table_ro$to_d11 <- c(ol_d1d11$score, ol_d1d11$pvalue)
      write.csv(table_ro, file = paste("./output/", runparams[[runtype]]$prefix, " - Table - overlap scores.csv", sep = ""))
      
      
# -----------------
# Research question: How similar is expression, quantitatively?
# -----------------
  # Rather than compare gene lists, which are based on categorical calls, let's look quantitatively at the degree of similarity in expression
  # To do this, we will look only at the transcripts regulated at d1, because otherwise our comparison will be swamped by noise among unregulated transcripts
  #    As before, no point in doing the same with 1-week transcipts as there will probably be only a handful of them, and even among those regulated most are at 1d too
  
    # 2a - similarity using genas
    # First, we'll use the genuine association (genas) function to calculate a correlation adjusted for sample reliability.
    # This approach is ok, but could be a bit optomistic
    # Calculate expression association between d1 and d5 among transcripts regulated at d1
    assoc_d1d5 <- genas(fit_1donly, coef=c("d1", "d5"), subset="all", plot=TRUE, alpha=0.4)
    # Now expression association between d1 and d11
    assoc_d1d11 <- genas(fit_1donly, coef=c("d1", "d11"), subset="all", plot=TRUE, alpha=0.4)
    # Then association between d5 and d11 for calculating difference in correlation
    assoc_d5d11 <- genas(fit_1donly, coef=c("d5", "d11"), subset="all", plot=TRUE, alpha=0.4)
    # Finally, expression association between d1 and d1rep
    assoc_d1d1rep <- genas(fit_1donly, coef=c("d1", "d1rep"), subset = "all", plot = TRUE, alpha = 0.4)
    # Print results
    assoc_d1d5
    assoc_d1d11
    assoc_d1d1rep
    
    #2b - similarity using standard Pearson's r - preferred approach
      # Just the basic correlation between expression at d1 and w1
          cor_d1d5 <- cor.test(x=fit_1donly$coefficients[, "d1"], y=fit_1donly$coefficients[, "d5"])
          cor_d1d5
      cor_d1d5 <- cor(x=fit_1donly$coefficients[, "d1"], y=fit_1donly$coefficients[, "d5"])
      # Same but for d1 and savings
          cor_d1d11 <- cor.test(x=fit_1donly$coefficients[, "d1"], y=fit_1donly$coefficients[, "d11"])
          cor_d1d11
      
      cor_d1d11 <- cor(x=fit_1donly$coefficients[, "d1"], y=fit_1donly$coefficients[, "d11"])
      # Also, calculate d5 to d11, because we need this to correctly compar d1d5 to d1d11
      cor_d5d11 <- cor(x=fit_1donly$coefficients[, "d5"], y=fit_1donly$coefficients[, "d11"])
      # Finally, the correlation between d1 and d1_rep
          cor_d1d1rep <- cor.test(x = fit_1donly$coefficients[, "d1"], y = fit_1donly$coefficients[, "d1rep"])
          cor_d1d1rep
      cor_d1d1rep <- cor(x = fit_1donly$coefficients[, "d1"], y = fit_1donly$coefficients[, "d1rep"])
      print(paste("Quantitative Comprarison: d1_d5 = ", cor_d1d5, " and d1_d11 = ", cor_d1d11, " giving a rdiff = ", cor_d1d5 - cor_d1d11))
    
    #2c - now calclulate difference in correlations, using both standard pearson's r and geneas corrected correlations
      r_result_uncorrected <- paired.r(xy= cor_d1d5, xz = cor_d1d11, yz = cor_d5d11, n = nrow(fit_1donly$coefficients))
      r_result_corrected <-paired.r(xy= assoc_d1d5$biological.correlation, xz = assoc_d1d11$biological.correlation, yz = assoc_d5d11$biological.correlation, n = nrow(fit_1donly$coefficients))
      
    
    # Make table of output
      table_r <- data.frame("To d1" = c("r", "r_corrected"))
      table_r$d5 <- c(cor_d1d5, assoc_d1d5$biological.correlation)
      table_r$d1_rep <- c(cor_d1d1rep, assoc_d1d1rep$biological.correlation)
      table_r$d11 <- c(cor_d1d11, assoc_d1d11$biological.correlation)
      table_r$rdiff <- c(
        cor_d1d5 - cor_d1d11,
        assoc_d1d5$biological.correlation - assoc_d1d11$biological.correlation
      )
      table_r$rdiff_p <- c(
          r_result_uncorrected$p,
          r_result_corrected$p
      )
      write.csv(table_r, file = paste("./output/", runparams[[runtype]]$prefix, " - Table - correlations.csv", sep = ""))
      
      
      # Make a scatter plot of d1 to d5
      # Basic plot using coefficients from fit_1donly, with x = d1 and y = d5
      d5plot <- ggplot(data = as.data.frame(fit_1donly$coefficients), aes(x = d1, y = d5))
      # Set limits and labels
      d5plot <- d5plot + xlim(-2, 2) + ylim(-2,2) + xlab("LFC, 1-Day") + ylab("LFC, 5 day")
        # TO BE CHANGED - might need to adjust the xlim and ylim
      # Add points
      d5plot <- d5plot + geom_point(alpha = 0.1)
      # Add regression line
      d5plot <- d5plot + geom_smooth(method=lm)
      # Add lines through 0 representing no change
      d5plot <- d5plot + geom_hline(yintercept = 0, linetype="dashed")
      d5plot <- d5plot + geom_vline(xintercept = 0, linetype="dashed")
      # Add line with slope = 1 indicating what a perfect expression correlation would look like
      d5plot <- d5plot + geom_abline(intercept = 0, slope = 1, linetype = "dotted")
      # Classic theme
      d5plot <- d5plot + theme_classic()
      # Print the r2 value on the plot
      d5plot <- d5plot + annotate(geom="text", x=-1.8,y=1.9,label = paste("r^2 = ", format(cor_d1d5^2, digits=2)))
      if(runtype != "Real") {
        d5plot <- d5plot + annotate("text", label="sample data", x = 0, y = 0, angle = 45, alpha = .1, size = 14)
      }
      
      d5plot
      ggsave(plot = d5plot, filename = paste("./output/", runparams[[runtype]]$prefix, " - d1 vs d5.jpg", sep = ""))    
      
  
      # Same plot for d1 to d11
      d11plot <- ggplot(data = as.data.frame(fit_1donly$coefficients), aes(x = d1, y = d11))
      d11plot <- d11plot + xlim(-2, 2) + ylim(-2,2) + xlab("LFC, 1-Day") + ylab("LFC, Day 11")
        # TO BE CHANGED - might need to adjust the xlim and ylim
      d11plot <- d11plot + geom_point(alpha = 0.1)
      d11plot <- d11plot + geom_smooth(method=lm)
      d11plot <- d11plot + geom_hline(yintercept = 0, linetype="dashed")
      d11plot <- d11plot + geom_vline(xintercept = 0, linetype="dashed")
      d11plot <- d11plot + geom_abline(intercept = 0, slope = 1, linetype = "dotted")
      d11plot <- d11plot + theme_classic()
      d11plot <- d11plot + annotate(geom="text", x=-1.8,y=1.9,label = paste("r^2 = ", format(cor_d1d11^2, digits=2)))
      if(runtype != "Real") {
        d11plot <- d11plot + annotate("text", label="sample data", x = 0, y = 0, angle = 45, alpha = .1, size = 14)
      }
      d11plot
      ggsave(plot = d11plot, filename = paste("./output/", runparams[[runtype]]$prefix, " - d1 vs d11.jpg", sep = ""))    
      
      
      # Same plot for d1 to d1rep
      repplot <- ggplot(data = as.data.frame(fit_1donly$coefficients), aes(x = d1, y = d1rep))
      repplot <- repplot + xlim(-2, 2) + ylim(-2,2) + xlab("LFC, 1-Day") + ylab("LFC, 1-Day Replication")
      # TO BE CHANGED - might need to adjust the xlim and ylim
      repplot <- repplot + geom_point(alpha = 0.1)
      repplot <- repplot + geom_smooth(method=lm)
      repplot <- repplot + geom_hline(yintercept = 0, linetype="dashed")
      repplot <- repplot + geom_vline(xintercept = 0, linetype="dashed")
      repplot <- repplot + geom_abline(intercept = 0, slope = 1, linetype = "dotted")
      repplot <- repplot + theme_classic()
      repplot <- repplot + annotate(geom="text", x=-1.8,y=1.9,label = paste("r^2 = ", format(cor_d1d1rep^2, digits=2)))
      if(runtype != "Real") {
        repplot <- repplot + annotate("text", label="sample data", x = 0, y = 0, angle = 45, alpha = .1, size = 14)
      }
      repplot
      ggsave(plot = repplot, filename = paste("./output/", runparams[[runtype]]$prefix, " - d1 vs d1replication.jpg", sep = ""))    
  

# -----------------------------------------------------------------------------
# The code below was added to the pre-registered script to organize a master
#   output file.  It does not change/alter/add to pre-registered analyses.
#
      
      # Grab the data from the MAavg object
      MAread<- as.data.frame(MA2avg$M)

      # Update column names to match sample names rather than slide file names
      colnames(MAread) <- targets$Name

      # Update MA table for dye swaps
      targets$swapweight <- 1
      targets[targets$Cy3 != "Control", ]$swapweight <- -1
      MAread <- as.data.frame(t(t(MAread) * targets$swapweight))
      
      # Now combine MA table with statistical analyses
      # Remake results table, just in case
      resultsTable <- list(
        d1 = topTreat(fit, coef="d1", number=nrow(fit), adjust.method="BH"),
        d5 = topTreat(fit, coef="d5", number=nrow(fit), adjust.method="BH"),
        d11 = topTreat(fit, coef="d11", number=nrow(fit), adjust.method="BH"),
        d1_rep = topTreat(fit, coef="d1rep", number=nrow(fit), adjust.method="BH")
      )
      int_list <- topTreat(int_fit, coef="d5vd1", number=nrow(int_fit), adjust.method="BH")
      int_list_d11 <- topTreat(int_fit_d11, coef="d11vd1", number=nrow(int_fit), adjust.method="BH")
      
      # Rename analysis column names
      colnames(resultsTable$d1) <- paste("d1_", colnames(resultsTable$d1), sep = "")
      colnames(resultsTable$d5) <- paste("d5_", colnames(resultsTable$d5), sep = "")
      colnames(resultsTable$d11) <- paste("d11_", colnames(resultsTable$d11), sep = "")
      colnames(resultsTable$d1_rep) <- paste("d1_rep_", colnames(resultsTable$d1_rep), sep = "")
      colnames(int_list) <- paste("int_d1_d5_", colnames(int_list), sep = "")
      colnames(int_list_d11) <- paste("int_d1_d11_", colnames(int_list_d11), sep = "")
      
      MAread$SystematicName <- row.names(MAread)
      
      
      MAread <- merge(MAread, resultsTable$d1[ , c(5, 7:11)], by.x = "SystematicName", by.y = "d1_SystematicName", all.y = TRUE, all.x = TRUE)
      MAread <- merge(MAread, resultsTable$d5[ , c(5, 7:11)], by.x = "SystematicName", by.y = "d5_SystematicName", all.y = TRUE, all.x = TRUE)
      MAread <- merge(MAread, resultsTable$d11[ , c(5, 7:11)], by.x = "SystematicName", by.y = "d11_SystematicName", all.y = TRUE, all.x = TRUE)
      MAread <- merge(MAread, resultsTable$d1_rep[ , c(5, 7:11)], by.x = "SystematicName", by.y = "d1_rep_SystematicName", all.y = TRUE, all.x = TRUE)
      MAread <- merge(MAread, int_list[ , c(5, 7:11)], by.x = "SystematicName", by.y = "int_d1_d5_SystematicName", all.y = TRUE, all.x = TRUE)
      MAread <- merge(MAread, int_list_d11[ , c(5, 7:11)], by.x = "SystematicName", by.y = "int_d1_d11_SystematicName", all.y = TRUE, all.x = TRUE)
      
      # Define which transcripts were previously identified as upregulated 1 
      #   day after a 1-day protocol
      restricted = read.csv("refined_24h.txt")
      MAread$previous <- MAread$SystematicName %in% restricted$Transcripts
      
      
      # Calculate own means for each time point
      MAread$d1_mean <- rowMeans(MAread[ , 2:8])
      MAread$d11_mean <- rowMeans(MAread[ , 9:16])
      MAread$d1_rep_mean <- rowMeans(MAread[ , 17:24])
      MAread$d5_mean <- rowMeans(MAread[ , 25:32])
      
      # Write the full output
      write.csv(MAread, file=paste("./output/", runparams[[runtype]]$prefix, " - all_output.csv", sep = ""))

      