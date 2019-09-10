# This sample code returns the query data in tab-separated values format, which LabKey then
# renders as HTML. Replace this code with your R script. See the Help tab for more details.
#write.table(labkey.data, file = "${tsvout:tsvfile}", sep = "\t", qmethod = "double", col.names=NA)

# UPDATE: return top 3 plots plus sum of ions plot (Sept 2014)

# UPDATE: new sort to get top 3 plots (March 2016)

# UPDATE: in cases where there are medium labeled peptides instead of light labeled peptides, 
#   'medium' Isotope Label is converted to 'light' Isotope Label to calculate variable/constant ratios (April 2016)
# For a reverse curve: H/L = H/M. For a forward curve: L/H = M/H.

# UPDATE: convert all sample_group values to Lo, Med, Hi based on first character (L, M, H) of sample_group (April 2016)
#
# UPDATE: For experiment 4, add other processing info:
# Sample Group -> temperature (4C, FT = freeze/thaw, -70 = freezer)
# Concentration -> time point: for 4C: 0,6,24 hours (may vary)
# for FT: count: 1 or 2 FT cycles (may vary)
# for -70: in weeks: 4 or 5 weeks (may vary)
# Peptide concentration: Should vary by peptide (but have a single concentration for each peptide) -> can check just in case
# Also want to look at difference between groups and at time 0 and maybe test it
#
# UPDATE: For experiment 4 3 in the new template, information of some columns has been changed.
# SampleGroup: remain the same (Control, 4C , FT and Frozen)
# Time:  the number of hours, days, weeks, months prior to analysis of the sample.
# TimeUnits: the units recorded in the Time parameter.
# FreezeThawCycles : the number of freeze thaw cycles.
# TimeAdjusted: for FT: 1 or 2 cycle 
# For Control, 4C and Frozen: (Time TimeUnits)


library(Cairo) # need for producing PNG image using Panoramax
library(Rlabkey)
library(dplyr)
library(ggplot2)

##variable that specifies whether the script takes a standalone csv file vs going through Panorama
standalone <- FALSE
inferReplicateNumber <- FALSE

# Multiple plot function - from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# Sort the Time levels in the order of Control, 4C, Frozen and FT
sortLevel <- function(timeFactors) {
    time_points <- as.character(unique(timeFactors))
    df_tmp <- data.frame(time_point<-as.character(), sample_group<-as.character(), sample_group_type<-as.numeric(), time_value<-as.numeric())
    for (time_point in time_points) {
      sample_group <- strsplit(time_point, " ")[[1]][1]
      if (sample_group == "Control") {
        sample_group_type <- 1
      } else if (sample_group == "Frozen") {
        sample_group_type <- 3
      } else if (sample_group == "FTx1") {
        sample_group_type <- 4
      } else if (sample_group == "FTx2") {
        sample_group_type <- 5
      }
      else {
        sample_group_type <- 2
      }
      time_value <- strsplit(time_point, " ")[[1]][3]
      if (suppressWarnings(!is.na(as.numeric(time_value)))) {
        time_value <- as.numeric(time_value)
        df_tmp1 <- data.frame(time_point=time_point, sample_group=sample_group, sample_group_type=sample_group_type, time_value=time_value)
      } else {
        df_tmp1 <- data.frame(time_point=time_point, sample_group=sample_group, sample_group_type=sample_group_type, time_value=0)
      }
      df_tmp <- rbind(df_tmp, df_tmp1)
    }
    # Sort by sample_group_type, then sort by time_value
    df_tmp <- df_tmp[with(df_tmp, order(sample_group_type, time_value)), ]
    as.character(df_tmp$time_point)
}

 

# ***** plot_QC function *****
plot_QC <- function(plot_fragment_ion_results, input_peptide_sequence, current_ion, times) {
  
  
  if (current_ion == 'all'){
    
    current_ion <- 'sum of ions'
  }
  
  
  plot_title <- paste(input_peptide_sequence, current_ion, sep='\n')
  
  temp <- as.character(plot_fragment_ion_results$sample_group)
  
  plot_fragment_ion_results$Time <- paste(temp, 
                                          " (", 
                                          plot_fragment_ion_results$time, ")",
                                          #plot_fragment_ion_results$Time_scale, ")",
                                          sep="")
  #plot_fragment_ion_results$Time <- gsub("1 cycles", "1 cycle", plot_fragment_ion_results$Time)
  
  
  ##factor Time in such a way that the 4C temps are first
  # Sort the Time based on the order of Control, 4C, Frozen and FT.
  sortedLevels <- sortLevel(factor(plot_fragment_ion_results$Time))
  
  #plot_fragment_ion_results$Time <- factor(plot_fragment_ion_results$Time,
  #                                         levels = c("4C (0 hours)",
  #                                                    "4C (4 hours)",
  #                                                    "4C (6 hours)",
  #                                                    "4C (24 hours)",
  #                                                    "-70 (4 weeks)",
  #                                                    "-70 (5 weeks)",
  #                                                    "-70 (6 weeks)",
  #                                                    "FT (1 cycle)",
  #                                                    "FT (2 cycles)"))
  plot_fragment_ion_results$Time <- factor(plot_fragment_ion_results$Time,
                                           levels = sortedLevels)
  
  colnames(plot_fragment_ion_results)[colnames(plot_fragment_ion_results) == "sample_group"] <- 
    "Condition"
  
  # Change the replicate into replicateNumber
  # In plot_fragment_ion_results, the column Replicate is used as replicate Number now.
  #for (sortedLevel in sortedLevels) {
  #  len1 <- length(plot_fragment_ion_results$replicate[plot_fragment_ion_results$Time==sortedLevel])
  #  plot_fragment_ion_results$Replicate[plot_fragment_ion_results$Time==sortedLevel] <- as.character(1:len1)
  #}

  ##set values that are zero to the smallest non-zero value divided by 2
  ##flag this as well
  plot_fragment_ion_results$zero_values <- FALSE
  plot_fragment_ion_results$zero_values[plot_fragment_ion_results$calculated_area_ratio == 0] <- TRUE
  plot_fragment_ion_results$calculated_area_ratio[plot_fragment_ion_results$calculated_area_ratio == 0] <-
    min(plot_fragment_ion_results$calculated_area_ratio[plot_fragment_ion_results$calculated_area_ratio > 0])/2
  
  ##get minimum and maximum calculated area ratios
  min_area_ratio <- min(plot_fragment_ion_results$calculated_area_ratio)
  max_area_ratio <- max(plot_fragment_ion_results$calculated_area_ratio)
  
  plot_fragment_ion_results$horiz_line <- min_area_ratio
  
  ##perform ANOVA and get p-value for hypothesis of no difference between groups
  lm_null <- lm(calculated_area_ratio ~ 1, data = plot_fragment_ion_results)
  lm_groups <- lm(calculated_area_ratio ~ Time, data = plot_fragment_ion_results)
  p_val <- signif(anova(lm_null, lm_groups)$"Pr(>F)"[2],2)
  
  plot_title <- paste(plot_title, "\n", "P-value from ANOVA: ", p_val, sep="")
  
  g <- ggplot(plot_fragment_ion_results, aes(x=Time, y=calculated_area_ratio, color=Condition,
                                             shape=Replicate)) +
    geom_point() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle(plot_title) +
    coord_trans(y="log10") +
    scale_y_continuous(limits = c(min_area_ratio*0.8,max_area_ratio*1.2)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab("Condition (time)") + ylab("Measured (area ratio) [log-scale]")

  ##if there were any zeros, add in a horizontal line at the smallest non-zero value divided by 2
  if(sum(plot_fragment_ion_results$zero_values) > 0)
  {
    g <- g + geom_hline(aes(yintercept=horiz_line),color="grey",linetype="dashed") 
    warning("Dashed line indicates values of 0 - They are added at the smallest non-zero value divided by 2.")
  }
  
  g
  
} 

if(standalone)
{
  all_data_files <- list.files("../test_data_new_template")[2:116]
  current_number <- 7
  current_data_file <- all_data_files[current_number-1]
  QC_set <- read.csv(paste("../test_data_new_template/",
                           current_data_file,
                           sep=""))
  
  print(current_data_file)
                           
  curve_type <- "forward"
  
  # Transform the column names to match those from embedded panorama query
  
  thenames = tolower(names(QC_set))
  thenames = gsub("\\.","", thenames) 
  names(QC_set) <- thenames
  
  # rename columns in QC_set dataframe (replace Panorama names with new names used by R script)
  QC_set <- dplyr::rename(QC_set,
                          peptide = peptidemodifiedsequence,
                          protein_name = protein,
                          precursor_charge = precursorcharge,
                          product_charge = productcharge,
                          fragment_ion_only = fragmention,
                          time_only = time, 
                          replicate = replicatename,
                          Replicate = replicatenumber,
                          sample_group = exp4samplegroup,
                          isotope_label_type = isotopelabel,
                          temperature = temperature,
                          time_units = timeunits,
                          freeze_thaw_cycles = freezethawcycles,
                          area = area)
  
  QC_set$fragment_ion <- paste(QC_set[ ,'fragment_ion_only'], " (", QC_set[ ,'product_charge'], "+)", sep='' )
  QC_set$time_units <- tolower(QC_set$time_units)
  time <-  sapply(1:nrow(QC_set), function(x) ifelse(QC_set[x,'sample_group'] == "FT", paste(QC_set[x,'freeze_thaw_cycles'], 'cycles', sep=' '), paste(QC_set[x,'time_only'], QC_set[x,'time_units'], sep=' ')))
  QC_set <- cbind(QC_set, time)
  QC_set$time <- gsub("^1 cycles", "1 cycle", QC_set$time)
  QC_set$time <- gsub("^0 hours", "0 hour", QC_set$time)
  QC_set$time <- gsub("^1 hours", "1 hour", QC_set$time)
  QC_set$time <- gsub("^1 days", "1 day", QC_set$time)
  QC_set$time <- gsub("^1 weeks", "1 week", QC_set$time)
  
  input_protein_name <- as.character(unique(QC_set$protein_name))
  input_peptide_sequence <- as.character(unique(QC_set$peptide))
  input_precursor_charge <- as.character(unique(QC_set$precursor_charge))
  
} else {
    #args <- commandArgs(trailingOnly = TRUE)
    #dataset_path <- args[1]
    #input_protein_name <- args[2]
    #input_peptide_sequence <- args[3]
    #input_precursor_charge <- args[4]
    #curve_type <- args[5]

    ##test it out with these parameters
    #dataset_path <- "CPTAC_TEST/Experiment_4_TEST"
    #input_protein_name <- "AURKB_HUMAN"
    #input_peptide_sequence <- "SNVQPTAAPGQK"
    #input_precursor_charge <- 2
    #curve_type <- "reverse"

    #QC_set <- labkey.selectRows(
    #    baseUrl="http://cptac-proliant-linux.esacinc.com/labkey",
    #    folderPath=paste("/CPTAC Assay Portal/",dataset_path,"/Stability",sep=""),
    #    schemaName="targetedms",
    #    queryName="QCAnalysisQuery",
    #    viewName="",
    #    colFilter=makeFilter(c("Protein", "EQUAL", input_protein_name),
    #                         c("PeptideModifiedSequence", "EQUAL", input_peptide_sequence),
    #                         c("PrecursorCharge","EQUAL",input_precursor_charge)),
    #    containerFilter=NULL
  #)

  #current_number = 1234
  
  input_peptide_sequence <- labkey.url.params$query.PeptideModifiedSequence
  input_protein_name <- labkey.url.params$query.Protein
  input_precursor_charge <- labkey.url.params$query.PrecursorCharge
  curve_type <- labkey.url.params$curve_type

  # ***** get data from Panorama *****

  QC_set <- labkey.data[
            labkey.data$peptidemodifiedsequence == input_peptide_sequence &
            labkey.data$protein == input_protein_name &
            labkey.data$precursorcharge == input_precursor_charge,
          ]


  #thenames = tolower(names(QC_set));
  #thenames = gsub(" ","", thenames) 
  #names(QC_set) <- thenames;
  
  #
  # Now back to your regularly scheduled program
  #
  
  #column names in Panorama query:
  #protein
  #peptidemodifiedsequence
  #isotopelabel
  #precursorcharge
  #productcharge
  #fragmention
  #area
  #replicatename
  #time
  #samplegroup
  #timeunits
  #freezethawcycles
  #id
  #runid_folder
  
  # sample row
  #1 YARS.IPI00007074 VDAQFGGIDQR heavy 2 1 y6 3573011.0 GO_QCorig_Broad_1000ng_Interlab_092412_031 3 2 Med 1.0 2 22199 #66fba526-16af-1031-a003-
  
  
  # rename columns in QC_set dataframe (replace Panorama names with new names used by R script)
  colnames(QC_set)[colnames(QC_set)=="peptidemodifiedsequence"] <- "peptide"
  colnames(QC_set)[colnames(QC_set)=="protein"] <- "protein_name"
  colnames(QC_set)[colnames(QC_set)=="replicatename"] <- "replicate"
  colnames(QC_set)[colnames(QC_set)=="replicatenumber"] <- "Replicate"
  colnames(QC_set)[colnames(QC_set)=="precursorcharge"] <- "precursor_charge"
  colnames(QC_set)[colnames(QC_set)=="productcharge"] <- "product_charge"
  colnames(QC_set)[colnames(QC_set)=="fragmention"] <- "fragment_ion_only"
  colnames(QC_set)[colnames(QC_set)=="time"] <- "time_only"      # day
  colnames(QC_set)[colnames(QC_set)=="exp4samplegroup"] <- "sample_group"     # temperature: possible values are 4C, -70, or FT
  colnames(QC_set)[colnames(QC_set)=="isotopelabel"] <- "isotope_label_type"     # light, heavy
  colnames(QC_set)[colnames(QC_set)=="area"] <- "area"
  colnames(QC_set)[colnames(QC_set)=="timeunits"] <- "time_units"
  colnames(QC_set)[colnames(QC_set)=="freezethawcycles"] <- "freeze_thaw_cycles"
  
  QC_set$fragment_ion <- paste(QC_set[ ,'fragment_ion_only'], " (", QC_set[ ,'product_charge'], "+)", sep='' )
  QC_set$time_units <- tolower(QC_set$time_units)
  time <-  sapply(1:nrow(QC_set), function(x) ifelse(QC_set[x,'sample_group'] == "FT", paste(QC_set[x,'freeze_thaw_cycles'], 'cycles', sep=' '), paste(QC_set[x,'time_only'], QC_set[x,'time_units'], sep=' ')))
  QC_set <- cbind(QC_set, time)
  QC_set$time <- gsub("^1 cycles", "1 cycle", QC_set$time)
  QC_set$time <- gsub("^0 hours", "0 hour", QC_set$time)
  QC_set$time <- gsub("^1 hours", "1 hour", QC_set$time)
  QC_set$time <- gsub("^1 days", "1 day", QC_set$time)
  QC_set$time <- gsub("^1 weeks", "1 week", QC_set$time)
  
}


# **********************************


now <- Sys.time()

#plot_filename <- paste(input_peptide_sequence, "_", trunc(as.numeric(now)), ".png", sep='' )

ion_category <- 'error'

# convert columns from character to numeric
QC_set[,'area'] <- as.numeric(as.character(QC_set[,'area']))

# remove factor version
QC_set[,'replicate'] <- as.character(QC_set[,'replicate'])
QC_set[,'Replicate'] <- as.character(QC_set[,'Replicate'])
QC_set[,'time'] <- as.character(QC_set[,'time'])
QC_set[,'fragment_ion'] <- as.character(QC_set[,'fragment_ion'])
QC_set[,'sample_group'] <- as.character(QC_set[,'sample_group'])
QC_set[,'isotope_label_type'] <- as.character(QC_set[,'isotope_label_type'])

# get a list of all unique fragment ions associated with current peptide
#fragment_ion_list <- sort(unique(QC_set[ , 'fragment_ion']))
fragment_ion_list <- unique(QC_set[ , 'fragment_ion'])


# get a list of all unique times associated with current peptide
times <- sort(unique(QC_set[ , 'time']))

# get a list of all unique sample groups (4C, FT, -70) associated with current peptide
sample_groups <- sort(unique(QC_set[ , 'sample_group']))

# Since the current Skyline templates don't have replicateNumer, so replicateNumbers will be assigned based on the column of replicate.
# A new column named Replicate will be used as replicate number.
if (inferReplicateNumber) {
    for (current_ion in fragment_ion_list) {
        for (current_sample in sample_groups) {
            for (current_time in times) {
                replicateName_list <- unique(QC_set$replicate[QC_set$fragment_ion==current_ion & QC_set$sample_group==current_sample & QC_set$time==current_time])
                replicateNumber <- 1
                for (replicateName in replicateName_list) {
                    lenTmp <- length(QC_set$replicate[QC_set$fragment_ion==current_ion & QC_set$sample_group==current_sample & QC_set$time==current_time & QC_set$replicate==replicateName])
                    QC_set$Replicate[QC_set$fragment_ion==current_ion & QC_set$sample_group==current_sample & QC_set$time==current_time & QC_set$replicate==replicateName] <- as.character(rep(replicateNumber, times=lenTmp))
                    replicateNumber <- replicateNumber + 1
                }
                #QC_set$Replicate[QC_set$fragment_ion==current_ion & QC_set$sample_group==current_sample & QC_set$time==current_time] <- as.character(1:lenTmp)
            }
        }
    }
    Replicates <- sort(unique(QC_set[ , 'Replicate']))
} else {
    Replicates <- sort(unique(QC_set[ , 'Replicate']))
}

# get a list of all unique replicates associated with current peptide
replicates <- sort(unique(QC_set[ , 'replicate']))


# new: April 2016
# *** for medium labeled peptides ***
isotope_label_types <- unique(QC_set[ , 'isotope_label_type'])


if(('light' %in% isotope_label_types) & ('medium' %in% isotope_label_types)) {
  
    stop("both light and medium isotope label found in dataset")
  
} else {
  
  QC_set$isotope_label_type[QC_set$isotope_label_type == "medium"] <- "light"
  
}


sum_light_area <- 0
sum_heavy_area <- 0

theoretical_area <- 0
measured_area <- 0

fragment_ion_results <- data.frame()



# # ***** prepare file to print PNG images *****

# image_frame_count <- length(fragment_ion_list)+1   # one for each fragment ion plus one more for sum

# only plot top 3 ions plus one more for sum
image_frame_count <- 4


# PNG using Panorama
# CairoPNG(filename="${imgout:QC_plot_png}", width=image_frame_count*400, height=400, bg="white", units="px", res=600)
CairoPNG(filename="${imgout:QC_plot_png}", width=image_frame_count*480, height=400, bg="white", units="px")
#CairoPNG(filename="${imgout:QC_plot_png}", width=800, height=600, bg="white")


# par(mfrow= c(1, image_frame_count))

# ********************************************





for (current_ion in fragment_ion_list) {


  for (current_time in times) {

    for (current_sample in sample_groups) {
    
      for (current_Rep in Replicates) {


        current_set_count <- 0
        light_area <- 0
        heavy_area <- 0
        theoretical_area <- 0
        measured_area <- 0
        calculated_area_ratio <- 0

        current_set <- QC_set[QC_set$fragment_ion==current_ion & QC_set$time==current_time & QC_set$sample_group==current_sample & QC_set$Replicate==current_Rep, ]
        
        #if (nrow(current_set) > 0) {
            if(nrow(current_set[current_set$isotope_label_type=='light', ]) > 1){

              stop("more than one light isotope")
            }

            if(nrow(current_set[current_set$isotope_label_type=='heavy', ]) > 1){

              stop("more than one heavy isotope")
            }


            current_set_count <- nrow(current_set)

            if (current_set_count == 2) {

              light_area <- current_set[current_set$isotope_label_type=='light', 'area' ]
              heavy_area <- current_set[current_set$isotope_label_type=='heavy', 'area' ]


              if(curve_type=='forward'){

                theoretical_area <- heavy_area
                measured_area <- light_area

              }
              else if (curve_type=='reverse'){

                theoretical_area <- light_area
                measured_area <- heavy_area

              }
              else{

                stop("invalid curve type")
              }



              if(theoretical_area==0 | is.na(theoretical_area) | is.na(measured_area)){
                calculated_area_ratio = NA
                
              }
              else {
                calculated_area_ratio <- measured_area/theoretical_area


              }

              fragment_ion_results <-  rbind(fragment_ion_results, data.frame(Peptide = input_peptide_sequence, Protein_Name = input_protein_name, Precursor_Charge = input_precursor_charge,
                                                                              fragment_ion = current_ion, time = current_time, sample_group = current_sample, Replicate = current_Rep, light_area=light_area, heavy_area=heavy_area,
                                                                              theoretical_area=theoretical_area, measured_area=measured_area,
                                                                              calculated_area_ratio=calculated_area_ratio, ion_category='individual') )
              
            }
            else {
            }
        #}
      } # end current_rep
    
    } # end current_sample
    
  } # end current_day

} # end current_ion


# ***** repeat calculations for sum of ions *****




for (current_time in times) {
  
  for (current_sample in sample_groups) {
    
    # Use replicate number here.
    #for (current_rep in replicates) {
    for (current_Rep in Replicates) {

      sum_light_area <- 0
      sum_heavy_area <- 0
      skipped_count <- 0
      skip_current_sample <- 'true'

      sum_theoretical_area <- 0
      sum_measured_area <- 0


      for (current_ion in fragment_ion_list)  {

        current_set_count <- 0
        calculated_area_ratio <- 0
        light_area <- 0
        heavy_area <- 0
        theoretical_area <- 0
        measured_area <- 0

        current_set <- QC_set[QC_set$fragment_ion==current_ion & QC_set$time==current_time & QC_set$sample_group==current_sample & QC_set$Replicate==current_Rep, ]
        
        #if (row(current_set) > 0) {

            if(nrow(current_set[current_set$isotope_label_type=='light', ]) > 1){

              stop("more than one light isotope")
            }

            if(nrow(current_set[current_set$isotope_label_type=='heavy', ]) > 1){

              stop("more than one heavy isotope")
            }


            current_set_count <- nrow(current_set)

            if (current_set_count == 2) {

              light_area <- current_set[current_set$isotope_label_type=='light', 'area' ]
              heavy_area <- current_set[current_set$isotope_label_type=='heavy', 'area' ]

              if(curve_type=='forward'){

                theoretical_area <- heavy_area
                measured_area <- light_area

              }
              else if (curve_type=='reverse'){

                theoretical_area <- light_area
                measured_area <- heavy_area

              }
              else{

                stop("invalid curve type")
              }

              add_to_sum <- 'false'

              if(theoretical_area==0 | is.na(theoretical_area) | is.na(measured_area)){
                skipped_count <- skipped_count + 1
              }
              else {
                sum_theoretical_area <- sum_theoretical_area + theoretical_area
                sum_measured_area <- sum_measured_area + measured_area

                sum_light_area <- sum_light_area + light_area
                sum_heavy_area <- sum_heavy_area + heavy_area

                add_to_sum <- 'true'
              }



              skip_current_sample <- 'false'


            }
            else {
            }
        #}

      } # end current_ion



      if (skip_current_sample=='false'){

        if(sum_theoretical_area==0){
          calculated_area_ratio <- NA

        }
        else {
          calculated_area_ratio <- sum_measured_area/sum_theoretical_area


        }



        fragment_ion_results <-  rbind(fragment_ion_results, data.frame(Peptide = input_peptide_sequence, Protein_Name = input_protein_name, Precursor_Charge = input_precursor_charge,
                                                                        fragment_ion = 'all', time = current_time, sample_group = current_sample, Replicate = current_Rep, light_area=sum_light_area, heavy_area=sum_heavy_area,
                                                                        theoretical_area=sum_theoretical_area, measured_area=sum_measured_area,
                                                                        calculated_area_ratio=calculated_area_ratio, ion_category='all') )
        

      }
      
    } # end current_rep
    
  } # end current_sample
  
} # end current_day




# ***** calculate CV (Coefficient of Variation) *****

conditions <- c("Control", "Autosampler", "Frozen", "FTx1", "FTx2")
##if condition happens to be something else, then replace 4C by that
for(t in 5:10)
{
  temp <- paste(t,"C",sep="")
  if(length(grep(temp, as.character(fragment_ion_results$sample_group))) > 0)
  {
    conditions[2] <- temp
  }
}
ions <- c(fragment_ion_list, 'all')
# days
# replicates


# make CV summary data frame
CV_results <- data.frame(fragment_ion = ions,
                         control_intra_CV = NA,
                         actual_temp_intra_CV = NA,
                         frozen_intra_CV = NA,
                         FTx1_intra_CV = NA,
                         FTx2_intra_CV = NA,
                         all_intra_CV = NA,
                         all_inter_CV = NA,
                         control_count = NA,
                         control_count_light_area_0 = NA,
                         control_count_heavy_area_0 = NA,
                         actual_temp_count = NA,
                         actual_temp_count_light_area_0 = NA,
                         actual_temp_count_heavy_area_0 = NA,
                         frozen_count = NA,
                         frozen_count_light_area_0 = NA,
                         frozen_count_heavy_area_0 = NA,
                         FTx1_count = NA,
                         FTx1_count_light_area_0 = NA,
                         FTx1_count_heavy_area_0 = NA,
                         FTx2_count = NA,
                         FTx2_count_light_area_0 = NA,
                         FTx2_count_heavy_area_0 = NA)
  
# *** intra-assay CV ***
for (current_ion in ions)  {

  individual_intra_assay_CVs_all <- c()
  
  for(current_cond in conditions) {
  
    avg_intra_assay_CV <-0
    individual_intra_assay_CVs <- c()

    for (current_time in times) {

      #for (current_rep in replicates) {

      current_set <- fragment_ion_results[fragment_ion_results$fragment_ion==current_ion & 
                                            fragment_ion_results$time==current_time &
                                            fragment_ion_results$sample_group==current_cond, ]

      # remove rows with a value of NA for calculated_area_ratio
      current_set <- current_set[complete.cases(current_set[ , 'calculated_area_ratio' ]),]

      if (nrow(current_set) <= 1 ){
        percent_CV <- NA
      }
      else{
        percent_CV <- (sd(current_set$calculated_area_ratio))/(mean(current_set$calculated_area_ratio)) * 100
        individual_intra_assay_CVs <- c(individual_intra_assay_CVs, percent_CV)
        individual_intra_assay_CVs_all <- c(individual_intra_assay_CVs_all, percent_CV)
      }

      #  } # end current_rep


    } # end current_day


    if (length(individual_intra_assay_CVs)==0){
      avg_CV <- NA
      count <- 0

    }
    else{
      avg_CV <- mean(individual_intra_assay_CVs, na.rm = TRUE)
    }
    
    if (current_cond=="Control") {
      CV_results[CV_results$fragment_ion==current_ion, "control_intra_CV"] <- avg_CV
    }
    else if(current_cond==conditions[2]){
      CV_results[CV_results$fragment_ion==current_ion, "actual_temp_intra_CV"] <- avg_CV
    }
    else if(current_cond=="Frozen"){
      CV_results[CV_results$fragment_ion==current_ion, "frozen_intra_CV"] <- avg_CV
    }
    else if(current_cond=="FTx1"){
      CV_results[CV_results$fragment_ion==current_ion, "FTx1_intra_CV"] <- avg_CV
    } else {
      CV_results[CV_results$fragment_ion==current_ion, "FTx2_intra_CV"] <- avg_CV
    }

  } #end current_cond
  
  ##now also get this over all conditions
  CV_results[CV_results$fragment_ion==current_ion, "all_intra_CV"] <- mean(individual_intra_assay_CVs_all, 
                                                                           na.rm = TRUE)
} # current_ion

# round for display purposes:
CV_results[,c("control_intra_CV", "actual_temp_intra_CV", "frozen_intra_CV", "FTx1_intra_CV", "FTx2_intra_CV", "all_intra_CV")] <- 
  sapply(CV_results[,c("control_intra_CV", "actual_temp_intra_CV", "frozen_intra_CV", "FTx1_intra_CV", "FTx2_intra_CV", "all_intra_CV")],
         round, digits=1)
  
# *** END: intra-assay CV ***

# *** get the number of samples with areas == 0 in each condition ***
for (current_ion in ions)  {

  for(current_cond in conditions) {
    
    current_set <- fragment_ion_results[fragment_ion_results$fragment_ion==current_ion & 
                                          fragment_ion_results$sample_group==current_cond, ]
    
    if (current_cond=="Control") {
      CV_results[CV_results$fragment_ion==current_ion,  
                 c("control_count_light_area_0",
                   "control_count_heavy_area_0")] <- 
        c(sum(current_set$light_area == 0),
          sum(current_set$heavy_area == 0))
    }
    if(current_cond==conditions[2]){
      CV_results[CV_results$fragment_ion==current_ion,  
                 c("actual_temp_count_light_area_0",
                   "actual_temp_count_heavy_area_0")] <- 
        c(sum(current_set$light_area == 0),
          sum(current_set$heavy_area == 0))
    }
    else if(current_cond=="Frozen"){
      CV_results[CV_results$fragment_ion==current_ion,  
                 c("frozen_count_light_area_0",
                   "frozen_count_heavy_area_0")] <- 
        c(sum(current_set$light_area == 0),
          sum(current_set$heavy_area == 0))
    }
    else if(current_cond=="FTx1"){
      CV_results[CV_results$fragment_ion==current_ion,  
                 c("FTx1_count_light_area_0",
                   "FTx1_count_heavy_area_0")] <- 
        c(sum(current_set$light_area == 0),
          sum(current_set$heavy_area == 0))
    } else {
      CV_results[CV_results$fragment_ion==current_ion,  
                 c("FTx2_count_light_area_0",
                   "FTx2_count_heavy_area_0")] <- 
        c(sum(current_set$light_area == 0),
          sum(current_set$heavy_area == 0))
    }
    
  } # end current_cond
  
} # end current_ion

# *** END: get the number of samples with areas == 0 in each condition ***

# *** inter-assay CV across all conditions ***
for (current_ion in ions)  {
  
  current_set <- fragment_ion_results[fragment_ion_results$fragment_ion==current_ion, ]
  
  ##for each replicate, get the mean and the standard deviation across all conditions
  mean_reps <- tapply(current_set$calculated_area_ratio, INDEX=current_set$Replicate,
                      FUN = mean)
  sd_reps <- tapply(current_set$calculated_area_ratio, INDEX=current_set$Replicate,
                    FUN = sd)
  
  percent_CV <- sd_reps/mean_reps * 100
  
  ##now calculate overall inter-assay CV!
  ##need to remove NAs in case any of them are left
  CV_results[CV_results$fragment_ion==current_ion, "all_inter_CV"] <- 
    mean(percent_CV, 
         na.rm = TRUE)
} # current_ion

# *** END: inter-assay CV across all conditions ***

# round for display purposes:
CV_results[,"all_inter_CV"] <- 
  sapply(CV_results[,"all_inter_CV"],
         round, digits=1)

# determine counts
for (current_ion in ions){
  
  CV_results[CV_results$fragment_ion==current_ion, "control_count"] <- 
    nrow(fragment_ion_results[fragment_ion_results$sample_group=="Control"&fragment_ion_results$fragment_ion==current_ion & !is.na(fragment_ion_results$calculated_area_ratio), ] )
  CV_results[CV_results$fragment_ion==current_ion, "actual_temp_count"] <- 
    nrow(fragment_ion_results[fragment_ion_results$sample_group==conditions[2]&fragment_ion_results$fragment_ion==current_ion & !is.na(fragment_ion_results$calculated_area_ratio), ] )
  CV_results[CV_results$fragment_ion==current_ion, "frozen_count"] <- 
    nrow(fragment_ion_results[fragment_ion_results$sample_group=="Frozen"&fragment_ion_results$fragment_ion==current_ion & !is.na(fragment_ion_results$calculated_area_ratio), ] )
  CV_results[CV_results$fragment_ion==current_ion, "FTx1_count"] <- 
    nrow(fragment_ion_results[fragment_ion_results$sample_group=="FTx1"&fragment_ion_results$fragment_ion==current_ion & !is.na(fragment_ion_results$calculated_area_ratio), ] )
  CV_results[CV_results$fragment_ion==current_ion, "FTx2_count"] <- 
    nrow(fragment_ion_results[fragment_ion_results$sample_group=="FTx2"&fragment_ion_results$fragment_ion==current_ion & !is.na(fragment_ion_results$calculated_area_ratio), ] )
}


ions_to_plot <- c()

# determine fragment ions to plot
if (length(ions) <= 4 ) {

  ions_to_plot <- ions


} else {
  
  #results_to_plot <- CV_results[CV_results$fragment_ion != 'all' & 
                                  #!is.na(CV_results$control_total_CV) & 
  #                                !is.na(CV_results$actual_temp_total_CV),
                                  #!is.na(CV_results$frozen_total_CV) & 
                                  #!is.na(CV_results$FT_total_CV) , 
  #                                ]
  results_to_plot <- CV_results[CV_results$fragment_ion != 'all',]
  
  # new sort to get Top 3 plots
  results_to_plot <- results_to_plot[order(results_to_plot$all_inter_CV), ]
  
  three_lowest_inter_CV <- head(results_to_plot, 3)
  
  three_lowest_ions <- as.character(three_lowest_inter_CV[ , 'fragment_ion'])
  
  ions_to_plot  <- c(three_lowest_ions, 'all')

}


#graphics.off()
#png(paste("../test_output/plots/Test_plot",current_number-1,".png",sep=""), width=image_frame_count*480, height=400)
# par(mfrow = c(1,length(ions_to_plot)))
all_plots <- list()
for (current_plot_ion in ions_to_plot ){

  ##if there are any NAs, throw a warning
  if(is.na(sum(fragment_ion_results$calculated_area_ratio)))
  {
    warning("Some are ratios could not be calculated. They are excluded from the plot and the tables.")
  }

  plot_fragment_ion_results <- fragment_ion_results[!is.na(fragment_ion_results$calculated_area_ratio), ]
  plot_fragment_ion_results <- plot_fragment_ion_results[plot_fragment_ion_results$fragment_ion == current_plot_ion, ]


  plot_times <- sort(unique(plot_fragment_ion_results[ , 'time']))


  # do not make plots for fragment ions with no data
  if (nrow(plot_fragment_ion_results) != 0) {

  # make QC plot for current ion
  current_plot <- plot_QC(plot_fragment_ion_results, input_peptide_sequence, current_plot_ion, plot_times)


  }

  all_plots[[current_plot_ion]] <- current_plot

}

multiplot(plotlist=all_plots, cols=image_frame_count)    

dev.off()

##add in comparison between interassay CVs for Experiment 4 and median concentration interassay CVs for Experiment 2
#if(current_number == 2)
#{
#  results_exp_2 <- read.table("../test_data_exp2/Test_CV_results1_experiment_2.tsv", header=TRUE)
#  CV_results$perc_diff_CV_exp2_4 <- try(round((results_exp_2$med_inter_CV - CV_results$all_inter_CV)/results_exp_2$med_inter_CV*100,2))
#} else {
#CV_results$perc_diff_CV_exp2_4 <- rep(NA, nrow(CV_results))
#}

# ***** output using Panorama *****

# output to web
write.table(CV_results, file = "${tsvout:CV_results}", sep = "\t", qmethod = "double", col.names=NA)
#write.table(fragment_ion_results, file = "${tsvout:fragment_ion_results}", sep = "\t", qmethod = "double", col.names=NA)
#write.table(QC_set, file = "${tsvout:QC_set}", sep = "\t", qmethod = "double", col.names=NA)

# output to files
write.csv(CV_results, file="${fileout:CV_results.csv}")
#write.csv(CV_results, file=paste("../test_output/tables_summary/Test_CV_results",current_number-1,".csv",sep=""))
#write.csv(QC_set, file="${fileout:cpatc_assay_portal_QC_input_set.csv}")
#write.csv(fragment_ion_results, file="${fileout:cpatc_assay_portal_QC_area_ratios.csv}")


#write.table(QC_set, file = "${tsvout:tsvfile}", sep = "\t", qmethod = "double", col.names=NA)

# *********************************




##dev.off()
