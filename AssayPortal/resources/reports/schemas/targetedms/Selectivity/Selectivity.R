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
# UPDATE: For experiment 3, add other processing info:
# SampleGroup = Spike level (0,5,10)
# Concentration = Cell line
# Replicate = 1 or 2
#
# UPDATE: For experiment 3 in new template. Information of some columns has been changed.
# AnalyteConcentration: three concentration levels (no spike and two analyte spikes)
# SampleGroup:  individual samples of the matrix of interest
# ReplicateNumber: replicate number


library(Cairo) # need for producing PNG image using Panoramax
library(Rlabkey)
library(dplyr)
library(ggplot2)


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

# ***** plot_QC function *****
plot_QC <- function(plot_fragment_ion_results, input_peptide_sequence, current_ion) {


  if (current_ion == 'all'){

    current_ion <- 'sum of ions'
  }


  plot_title <- paste(input_peptide_sequence, current_ion, sep='\n')

  # Expand right side of clipping rect to make room for the legend
  #par(xpd=TRUE, mar=par()$mar+c(0,0,0,4))

  # bty="L",
  #
  # plot(plot_fragment_ion_results$day, plot_fragment_ion_results$calculated_area_ratio, log="y", yaxt="n",
  #      col=ifelse(plot_fragment_ion_results$sample_group=='4C', 'red',
  #                 ifelse(plot_fragment_ion_results$sample_group=='FT', 'blue', 'green')),
  #      pch=ifelse(plot_fragment_ion_results$replicate==1 , 1,
  #                 ifelse(plot_fragment_ion_results$replicate==2 , 0,
  #                        ifelse(plot_fragment_ion_results$replicate==3 , 2,
  #                               ifelse(plot_fragment_ion_results$replicate==4 , 5,
  #                                      ifelse(plot_fragment_ion_results$replicate==5 , 6,
  #                               4))))),
  #
  #      cex=2, lwd=1, main=plot_title, xlab="Time (day)", ylab="Measured (area ratio) [log-scale]",
  #      cex.lab=1.75, cex.axis=1.75,
  #      ylim = c(min_area_ratio*0.8,max_area_ratio*1.2))

  ##average results for each replicate for the same cell line and spike level
  plot_grouped_reps <- plot_fragment_ion_results %>% group_by(cell_line, spike_level)
  plot_ave_reps <- plot_grouped_reps %>%
    summarize(calculated_area_ratio_ave_reps = mean(calculated_area_ratio))

  ##set values that are zero to the smallest non-zero value divided by 2
  ##flag this as well
  ##actually don't need this since it will be on linear scale!
  # plot_ave_reps$zero_values <- FALSE
  # plot_ave_reps$zero_values[plot_ave_reps$calculated_area_ratio_ave_reps == 0] <- TRUE
  # plot_ave_reps$calculated_area_ratio_ave_reps[plot_ave_reps$calculated_area_ratio_ave_reps == 0] <-
  #   min(plot_ave_reps$calculated_area_ratio_ave_reps[plot_ave_reps$calculated_area_ratio_ave_reps > 0])/2
  #
  ##get minimum and maximum calculated area ratios
  min_area_ratio <- min(plot_ave_reps$calculated_area_ratio_ave_reps)
  max_area_ratio <- max(plot_ave_reps$calculated_area_ratio_ave_reps)

  plot_ave_reps$horiz_line <- min_area_ratio

  g <- ggplot(plot_ave_reps, aes(x=spike_level, y=calculated_area_ratio_ave_reps,
                                 color=cell_line, group=cell_line)) +
    geom_point(size=2.5) +
    geom_smooth(method="lm",linetype="dashed",se = FALSE) +
    ##geom_line(aes(group=cell_line),linetype="dashed") +
    ggtitle(plot_title) +
    ##coord_trans(y="log10") +
    scale_y_continuous(limits = c(min_area_ratio*0.8, max_area_ratio*1.2)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    guides(color = guide_legend(title="Cell line")) +
    xlab("Spike-in level") + ylab("Measured (area ratio) [linear scale]")

  # ##if there were any zeros, add in a horizontal line at the smallest non-zero value divided by 2
  # ##actually don't need this since it will be on linear scale!
  # if(sum(plot_ave_reps$zero_values) > 0)
  # {
  #   g <- g + geom_hline(aes(yintercept=horiz_line),color="grey",linetype="dashed")
  #   warning("Dashed line indicates values of 0 - They are added at the smallest non-zero value divided by 2.")
  # }

  g

  #x_axis_values <- c(1,2,3,4,5)
  #x_axis_values <- days

  #y_axis_values <- c(0.01,0.1,1,10,100,format(max_calculated_measured_concentration,digits=3))
  # y_axis_values <- c(format(min_calculated_area_ratio,digits=3),format(median(plot_fragment_ion_results$calculated_area_ratio),digits=3),format(max_calculated_area_ratio,digits=3))
  #
  # #format(y_axis_values,scientific=FALSE,digits=4)
  #
  # axis(2, y_axis_values, labels=format(y_axis_values,scientific=FALSE)) # draw y axis with required labels
  #
  # #par(xpd=TRUE)
  # #legend(x=4,y=max_calculated_area_ratio,legend=c("rep1","rep2","rep3","repX","Hi","Med","Lo"),pch=c(1,0,2,4,18,18,18), col=c("black","black","black","black","red","blue","green"), bty="n")
  #
  # legend(x=max(days)+0.2,y=max_calculated_area_ratio,legend=c("rep1","rep2","rep3","rep4","rep5","4C","FT","-70"),pch=c(1,0,2,5,6,18,18,18), col=c("black","black","black","black","black","red","blue","green"), cex=1, bty="n")
  #
  # Restore default clipping rect
  # par(mar=c(5, 4, 4, 2) + 0.1)



}


#
# Now back to your regularly scheduled program
#

input_peptide_sequence <- labkey.url.params$query.PeptideModifiedSequence
input_protein_name <- labkey.url.params$query.Protein
input_precursor_charge <- labkey.url.params$query.PrecursorCharge
curve_type <- labkey.url.params$curve_type

# *****************************************************************************



# ***** input parameters to Panorama from URL *****

#input_peptide_sequence <- labkey.url.params$input_peptide_sequence
#input_protein_name <- labkey.url.params$input_protein_name
#input_precursor_charge <- labkey.url.params$input_precursor_charge
#curve_type <- labkey.url.params$curve_type

# *************************************************



# ***** get data from Panorama *****

QC_set <- labkey.data[
            labkey.data$peptidemodifiedsequence == input_peptide_sequence &
            labkey.data$protein == input_protein_name &
            labkey.data$precursorcharge == input_precursor_charge,
          ]

#column names in Panorama query:
#protein
#peptidemodifiedsequence
#isotopelabel
#precursorcharge
#productcharge
#fragmention
#area
#replicatename
#replicatenumber
#analyteconcentration
#samplegroup
#dilutionfactor
#peptideconcentration
#id
#runid_folder

# sample row
#1 YARS.IPI00007074 VDAQFGGIDQR heavy 2 1 y6 3573011.0 GO_QCorig_Broad_1000ng_Interlab_092412_031 3 2 Med 1.0 2 22199 #66fba526-16af-1031-a003-


# rename columns in QC_set dataframe (replace Panorama names with new names used by R script)
colnames(QC_set)[colnames(QC_set) == "peptidemodifiedsequence"] <- "peptide"
colnames(QC_set)[colnames(QC_set) == "protein"] <- "protein_name"
colnames(QC_set)[colnames(QC_set) == "replicatename"] <- "replicate_name"
colnames(QC_set)[colnames(QC_set) == "precursorcharge"] <- "precursor_charge"
colnames(QC_set)[colnames(QC_set) == "productcharge"] <- "product_charge"
colnames(QC_set)[colnames(QC_set) == "fragmention"] <- "fragment_ion_only"
colnames(QC_set)[colnames(QC_set) == "analyteconcentration"] <- "analyte_concentration"      # three concentration levels (no spike and two analyte spikes)
colnames(QC_set)[colnames(QC_set) == "replicatenumber"] <- "replicate" # recplicate number
colnames(QC_set)[colnames(QC_set) == "exp3samplegroup"] <- "sample_group"     # individual samples of the matrix of interest
colnames(QC_set)[colnames(QC_set) == "isotopelabel"] <- "isotope_label_type"     # light, heavy
colnames(QC_set)[colnames(QC_set) == "area"] <- "area"

QC_set$fragment_ion <- paste(QC_set[, 'fragment_ion_only'], " (", QC_set[, 'product_charge'], "+)", sep = '')


# **********************************


now <- Sys.time()

plot_filename <- paste(input_peptide_sequence, "_", trunc(as.numeric(now)), ".png", sep='' )

ion_category <- 'error'


# convert some classes for some columns
QC_set[,'sample_group'] <- as.character(QC_set[,'sample_group'])

# remove factor version
QC_set[,'fragment_ion'] <- as.character(QC_set[,'fragment_ion'])
QC_set[,'isotope_label_type'] <- as.character(QC_set[,'isotope_label_type'])

# get a list of all unique fragment ions associated with current peptide
#fragment_ion_list <- sort(unique(QC_set[ , 'fragment_ion']))
fragment_ion_list <- unique(QC_set[ , 'fragment_ion'])

# get a list of all unique cell lines associated with current peptide
cell_lines <- sort(unique(QC_set[ , 'sample_group']))

# get a list of all unique spike levels (0,5,10) associated with current peptide
spike_levels <- sort(unique(QC_set[ , 'analyte_concentration']))

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

#image_frame_count <- length(fragment_ion_list)+1   # one for each fragment ion plus one more for sum

# only plot top 3 ions plus one more for sum
image_frame_count <- 4


# PNG using Panorama
# CairoPNG(filename="${imgout:QC_plot_png}", width=image_frame_count*400, height=400, bg="white", units="px", res=600)
CairoPNG(filename="${imgout:QC_plot_png}", width=image_frame_count*400, height=400, bg="white", units="px")
#CairoPNG(filename="${imgout:QC_plot_png}", width=800, height=600, bg="white")


# par(mfrow= c(1, image_frame_count))

# ********************************************





for (current_ion in fragment_ion_list) {

  for (current_cell_line in cell_lines) {

    for (current_spike_level in spike_levels) {

      for (current_rep in replicates) {


        current_set_count <- 0
        light_area <- 0
        heavy_area <- 0
        theoretical_area <- 0
        measured_area <- 0
        calculated_area_ratio <- 0

        current_set <- QC_set[QC_set$fragment_ion==current_ion & QC_set$sample_group==current_cell_line & QC_set$analyte_concentration==current_spike_level & QC_set$replicate==current_rep, ]

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
                                                                          fragment_ion = current_ion, cell_line = current_cell_line, spike_level = current_spike_level, replicate = current_rep, light_area=light_area, heavy_area=heavy_area,
                                                                          theoretical_area=theoretical_area, measured_area=measured_area,
                                                                          calculated_area_ratio=calculated_area_ratio, ion_category='individual') )

        }
        else {


        }


      } # end current_rep

    } # end current_sample

  } # end current_cell_line


} # end current_ion






# ***** repeat calculations for sum of ions *****




for (current_cell_line in cell_lines) {

  for (current_spike_level in spike_levels) {

    for (current_rep in replicates) {

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



        current_set <- QC_set[QC_set$fragment_ion==current_ion & QC_set$sample_group==current_cell_line & QC_set$analyte_concentration==current_spike_level & QC_set$replicate==current_rep, ]


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


      } # end current_ion



      if (skip_current_sample=='false'){

        if(sum_theoretical_area==0){
          calculated_area_ratio <- NA

        }
        else {
          calculated_area_ratio <- sum_measured_area/sum_theoretical_area


        }




        fragment_ion_results <-  rbind(fragment_ion_results, data.frame(Peptide = input_peptide_sequence, Protein_Name = input_protein_name, Precursor_Charge = input_precursor_charge,
                                                                        fragment_ion = 'all', cell_line = current_cell_line, spike_level = current_spike_level, replicate = current_rep, light_area=light_area, heavy_area=heavy_area,
                                                                        theoretical_area=theoretical_area, measured_area=measured_area,
                                                                        calculated_area_ratio=calculated_area_ratio, ion_category='all') )


      }

    } # end current_cell_line

  } # end current_spike_level

} # end current_rep




# ***** calculate CV (Coefficient of Variation) *****

# conditions <- c("4C", "-70", "FT")
# ##if condition happens to be something else, then replace 4C by that
# for(t in 5:10)
# {
#   temp <- paste(t,"C",sep="")
#   if(length(grep(temp, as.character(fragment_ion_results$analyte_concentration))) > 0)
#   {
#     conditions[1] <- temp
#   }
# }
ions <- c(fragment_ion_list, 'all')
# # days
# # replicates
# 
# 
# make summary data tables
summary_table <- data.frame(matrix(NA, ncol=length(cell_lines)+1, nrow=length(ions)))
colnames(summary_table) <- c("fragment_ion", cell_lines)
summary_table$fragment_ion <- ions

# make summary data tables in another format
summary_table2 <- data.frame(fragment_ion=c(),
                                      cell_line=c(),
                                      estimated_slope=c())

## also save the 3 data points in a file
values_for_spike_levels <- data.frame(fragment_ion=c(),
                                      cell_line=c(),
                                      spike_level=c(),
                                      calculated_rea_ratio_ave_rep=c())

##do this for each ion
for(current_ion in ions)
{
  plot_fragment_ion_results <- fragment_ion_results[!is.na(fragment_ion_results$calculated_area_ratio), ]
  plot_fragment_ion_results <- plot_fragment_ion_results[plot_fragment_ion_results$fragment_ion == current_ion, ]

  ##average results for each replicate for the same cell line and spike level
  plot_grouped_reps <- plot_fragment_ion_results %>% group_by(cell_line, spike_level)
  plot_ave_reps <- plot_grouped_reps %>%
    summarize(calculated_area_ratio_ave_reps = mean(calculated_area_ratio))
  plot_ave_reps$cell_line <- as.character(plot_ave_reps$cell_line)

  ##for each cell line, get the estimated slope
  for(current_cell_line in cell_lines)
  {
    cell_line_subset <- plot_ave_reps %>% dplyr::filter(cell_line == current_cell_line)
    estimated_slope <- coef(lm(calculated_area_ratio_ave_reps ~ spike_level, data=cell_line_subset))["spike_level"]
    summary_table[summary_table$fragment_ion == current_ion,current_cell_line] <- estimated_slope
    v_current_ion2 <- data.frame(fragment_ion=current_ion, cell_line=current_cell_line, estimated_slope=estimated_slope)
    summary_table2 <- rbind(summary_table2, v_current_ion2)
  }

  v_current_ion <- cbind(fragment_ion = rep(current_ion, nrow(plot_ave_reps)),
                         plot_ave_reps)
  values_for_spike_levels <- rbind(values_for_spike_levels,
                                   v_current_ion)
}

rownames(summary_table2) <- 1:nrow(summary_table2)

##save values for spike levels
#write.csv(values_for_spike_levels,
#          file=paste("../test_output/tables_values_for_points/Test_ave_values_for_spike_levels",current_number-1,".csv",sep=""))

write.table(values_for_spike_levels, file = "${tsvout:values_for_spike_levels}", sep = "\t", qmethod = "double", col.names=NA)
write.csv(values_for_spike_levels, file="${fileout:values_for_spike_levels.csv}")

ions_to_plot <- ions

# determine fragment ions to plot
if (length(ions) <= 4 ) {

  ions_to_plot <- ions

} else {

  ##get the ions that have the highest median area across the entire dataset
  ##first get the median area ratio by ion
  median_area_per_ion <- tapply(fragment_ion_results$calculated_area_ratio,
                                INDEX = fragment_ion_results$fragment_ion,
                                FUN=median,
                                na.rm=TRUE)
  ##remove the "all" entry (since that gets included anyway)
  median_area_per_ion <- median_area_per_ion[names(median_area_per_ion)!="all"]

  three_highest_area_ratios <- names(sort(median_area_per_ion, decreasing = TRUE))[1:3]

  ions_to_plot  <- c(three_highest_area_ratios, 'all')

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


  # do not make plots for fragment ions with no data
  if (nrow(plot_fragment_ion_results) != 0) {

  # make QC plot for current ion
  current_plot <- plot_QC(plot_fragment_ion_results, input_peptide_sequence, current_plot_ion)


  }

  all_plots[[current_plot_ion]] <- current_plot

}

multiplot(plotlist=all_plots, cols=image_frame_count)

dev.off()



# ***** output using Panorama *****

# output to web
write.table(summary_table, file = "${tsvout:summary_table}", sep = "\t", qmethod = "double", col.names=NA)
write.table(summary_table2, file = "${tsvout:summary_table2}", sep = "\t", qmethod = "double", col.names=NA)
#write.table(fragment_ion_results, file = "${tsvout:fragment_ion_results}", sep = "\t", qmethod = "double", col.names=NA)
#write.table(QC_set, file = "${tsvout:QC_set}", sep = "\t", qmethod = "double", col.names=NA)

# output to files
write.csv(summary_table, file="${fileout:summary_table.csv}")
write.csv(summary_table2, file="${fileout:summary_table2.csv}")
#write.csv(format(summary_table, digits=4),
#          file=paste("../test_output/tables_summary/Test_summary_table",current_number-1,".csv",sep=""))
#write.csv(QC_set, file="${fileout:cpatc_assay_portal_QC_input_set.csv}")
#write.csv(fragment_ion_results, file="${fileout:cpatc_assay_portal_QC_area_ratios.csv}")


#write.table(QC_set, file = "${tsvout:tsvfile}", sep = "\t", qmethod = "double", col.names=NA)

# *********************************




##dev.off()
