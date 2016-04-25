# This sample code returns the query data in tab-separated values format, which LabKey then
# renders as HTML. Replace this code with your R script. See the Help tab for more details.
#write.table(labkey.data, file = "${tsvout:tsvfile}", sep = "\t", qmethod = "double", col.names=NA)

# UPDATE: return top 3 plots plus sum of ions plot (Sept 2014)

# UPDATE: new sort to get top 3 plots (March 2016)

# UPDATE: in cases where there are medium labeled peptides instead of light labeled peptides, 
#   'medium' Isotope Label is converted to 'light' Isotope Label to calculate variable/constant ratios (April 2016)
# For a reverse curve: H/L = H/M. For a forward curve: L/H = M/H.

# UPDATE: convert all sample_group values to Lo, Med, Hi based on first character (L, M, H) of sample_group (April 2016)

library(Cairo) # need for producing PNG image using Panoramax



# ***** plot_QC function *****
plot_QC <- function(plot_fragment_ion_results, input_peptide_sequence, current_ion, days) {


  if (current_ion == 'all'){

      current_ion <- 'sum of ions'
  }


  plot_title <- paste(input_peptide_sequence, current_ion, sep='\n')

  min_calculated_area_ratio <- min(plot_fragment_ion_results$calculated_area_ratio)
  max_calculated_area_ratio <- max(plot_fragment_ion_results$calculated_area_ratio)



  # Expand right side of clipping rect to make room for the legend
  par(xpd=TRUE, mar=par()$mar+c(0,0,0,4))

  # bty="L",

  plot(plot_fragment_ion_results$day, plot_fragment_ion_results$calculated_area_ratio, log="y", yaxt="n", col=ifelse(plot_fragment_ion_results$sample_group=='Hi', 'red',
                                                                                                                     ifelse(plot_fragment_ion_results$sample_group=='Med', 'blue', 'green')),

       pch=ifelse(plot_fragment_ion_results$replicate==1 , 1,
                  ifelse(plot_fragment_ion_results$replicate==2 , 0,
                         ifelse(plot_fragment_ion_results$replicate==3 , 2,
                                ifelse(plot_fragment_ion_results$replicate==4 , 5,
                                       ifelse(plot_fragment_ion_results$replicate==5 , 6,
                                4))))),

       cex=2, lwd=1, main=plot_title, xlab="Time (day)", ylab="Measured (area ratio) [log-scale]", cex.lab=1.75, cex.axis=1.75)


  #x_axis_values <- c(1,2,3,4,5)
  x_axis_values <- days

  #y_axis_values <- c(0.01,0.1,1,10,100,format(max_calculated_measured_concentration,digits=3))
  y_axis_values <- c(format(min_calculated_area_ratio,digits=3),format(median(plot_fragment_ion_results$calculated_area_ratio),digits=3),format(max_calculated_area_ratio,digits=3))




  #format(y_axis_values,scientific=FALSE,digits=4)

  axis(2, y_axis_values, labels=format(y_axis_values,scientific=FALSE)) # draw y axis with required labels

  #par(xpd=TRUE)
  #legend(x=4,y=max_calculated_area_ratio,legend=c("rep1","rep2","rep3","repX","Hi","Med","Lo"),pch=c(1,0,2,4,18,18,18), col=c("black","black","black","black","red","blue","green"), bty="n")

  legend(x=max(days)+0.2,y=max_calculated_area_ratio,legend=c("rep1","rep2","rep3","rep4","rep5","Hi","Med","Lo"),pch=c(1,0,2,5,6,18,18,18), col=c("black","black","black","black","black","red","blue","green"), cex=1, bty="n")

  # Restore default clipping rect
  par(mar=c(5, 4, 4, 2) + 0.1)



} # end plot_QC function






# ***** input parameters *****

# test set 1
#input_peptide_sequence <- 'ASDAEGDVAALNR'
#input_protein_name <- 'TPM4.IPI00216975'
#input_precursor_charge <- '2'
#curve_type <- 'reverse'

# test set 2
#input_peptide_sequence <- 'VLGELWPLFGGR'
#input_protein_name <- 'ABCD3.IPI00002372'
#input_precursor_charge <- '2'
#curve_type <- 'reverse'

# ****************************


# ***** input parameters to Panorama from URL (data filtering from query) *****

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

QC_set <- labkey.data[labkey.data$peptidemodifiedsequence==input_peptide_sequence & labkey.data$protein==input_protein_name & labkey.data$precursorcharge==input_precursor_charge, ]

#column names in Panorama query:
#protein
#peptidemodifiedsequence
#isotopelabel
#precursorcharge
#productcharge
#fragmention
#area
#replicatename
#replicate
#concentration
#samplegroup
#dilutionfactor
#peptideconcentration
#id
#runid_folder

# sample row
#1 YARS.IPI00007074 VDAQFGGIDQR heavy 2 1 y6 3573011.0 GO_QCorig_Broad_1000ng_Interlab_092412_031 3 2 Med 1.0 2 22199 #66fba526-16af-1031-a003-


# rename columns in QC_set dataframe (replace Panorama names with new names used by R script)
colnames(QC_set)[colnames(QC_set)=="peptidemodifiedsequence"] <- "peptide"
colnames(QC_set)[colnames(QC_set)=="protein"] <- "protein_name"
colnames(QC_set)[colnames(QC_set)=="replicatename"] <- "replicate_name"
colnames(QC_set)[colnames(QC_set)=="precursorcharge"] <- "precursor_charge"
colnames(QC_set)[colnames(QC_set)=="productcharge"] <- "product_charge"
colnames(QC_set)[colnames(QC_set)=="fragmention"] <- "fragment_ion_only"
colnames(QC_set)[colnames(QC_set)=="concentration"] <- "concentration"      # day
colnames(QC_set)[colnames(QC_set)=="replicate"] <- "replicate"
colnames(QC_set)[colnames(QC_set)=="samplegroup"] <- "sample_group"     # Lo, Med, Hi
colnames(QC_set)[colnames(QC_set)=="isotopelabel"] <- "isotope_label_type"     # light, heavy
colnames(QC_set)[colnames(QC_set)=="area"] <- "area"

QC_set$fragment_ion <- paste(QC_set[ ,'fragment_ion_only'], " (", QC_set[ ,'product_charge'], "+)", sep='' )

# **********************************





now <- Sys.time()

plot_filename <- paste(input_peptide_sequence, "_", trunc(as.numeric(now)), ".png", sep='' )

ion_category <- 'error'


# convert columns from character to numeric
QC_set[,'concentration'] <- as.numeric(as.character(QC_set[,'concentration']))
QC_set[,'replicate'] <- as.numeric(as.character(QC_set[,'replicate']))
QC_set[,'area'] <- as.numeric(as.character(QC_set[,'area']))


# remove factor version
QC_set[,'fragment_ion'] <- as.character(QC_set[,'fragment_ion'])
QC_set[,'sample_group'] <- as.character(QC_set[,'sample_group'])
QC_set[,'isotope_label_type'] <- as.character(QC_set[,'isotope_label_type'])

# new: April 2016
# convert all sample_group values to Lo, Med, Hi based on first character of sample_group
QC_set[toupper(substr(QC_set$sample_group, 1, 1))=="L","sample_group"] <- "Lo"
QC_set[toupper(substr(QC_set$sample_group, 1, 1))=="M","sample_group"] <- "Med"
QC_set[toupper(substr(QC_set$sample_group, 1, 1))=="H","sample_group"] <- "Hi"


# get a list of all unique fragment ions associated with current peptide
#fragment_ion_list <- sort(unique(QC_set[ , 'fragment_ion']))
fragment_ion_list <- unique(QC_set[ , 'fragment_ion'])



# get a list of all unique days associated with current peptide
days <- sort(unique(QC_set[ , 'concentration']))


# get a list of all unique sample groups (Lo, Med, Hi) associated with current peptide
sample_groups <- sort(unique(QC_set[ , 'sample_group']))


# get a list of all unique replicates associated with current peptide
replicates <- sort(unique(QC_set[ , 'replicate']))


# new: April 2016
# *** for medium labled peptides ***
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
CairoPNG(filename="${imgout:QC_plot_png}", width=image_frame_count*400, height=400, bg="white", units="px")
#CairoPNG(filename="${imgout:QC_plot_png}", width=800, height=600, bg="white")


par(mfrow= c(1, image_frame_count))

# ********************************************





for (current_ion in fragment_ion_list) {


  for (current_day in days) {

    for (current_sample in sample_groups) {

      for (current_rep in replicates) {


        current_set_count <- 0
        light_area <- 0
        heavy_area <- 0
        theoretical_area <- 0
        measured_area <- 0
        calculated_area_ratio <- 0



        current_set <- QC_set[QC_set$fragment_ion==current_ion & QC_set$concentration==current_day & QC_set$sample_group==current_sample & QC_set$replicate==current_rep, ]


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
                                                                          fragment_ion = current_ion, day = current_day, sample_group = current_sample, replicate = current_rep, light_area=light_area, heavy_area=heavy_area,
                                                                          theoretical_area=theoretical_area, measured_area=measured_area,
                                                                          calculated_area_ratio=calculated_area_ratio, ion_category='individual') )

        }
        else {


        }


      } # end current_rep

    } # end current_sample

  } # end current_day




} # end current_ion






# ***** repeat calculations for sum of ions *****




for (current_day in days) {


  for (current_sample in sample_groups) {

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



        current_set <- QC_set[QC_set$fragment_ion==current_ion & QC_set$concentration==current_day & QC_set$sample_group==current_sample & QC_set$replicate==current_rep, ]


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
                                                                        fragment_ion = 'all', day = current_day, sample_group = current_sample, replicate = current_rep, light_area=sum_light_area, heavy_area=sum_heavy_area,
                                                                        theoretical_area=sum_theoretical_area, measured_area=sum_measured_area,
                                                                        calculated_area_ratio=calculated_area_ratio, ion_category='all') )

      }

    } # end current_rep

  } # end current_sample

} # end current_day




# ***** calculate CV (Coefficient of Variation) *****


LoMedHi <- c('Lo', 'Med', 'Hi')
ions <- c(fragment_ion_list, 'all')
# days
# replicates


# make CV summary data frame
CV_results <- data.frame()

for (current_ion in ions) {

  # define columns (make row for each fragment ion)
  CV_results <-  rbind(CV_results, data.frame(fragment_ion=current_ion, low_intra_CV=NA, med_intra_CV=NA, high_intra_CV=NA,
                                              low_inter_CV=NA, med_inter_CV=NA, high_inter_CV=NA,
                                              low_total_CV=NA, med_total_CV=NA, high_total_CV=NA,
                                              low_count=NA, med_count=NA, high_count=NA))
}




# *** intra-assay CV ***
for (current_ion in ions)  {

  for (current_LoMedHi in LoMedHi) {

    avg_intra_assay_CV <-0
    individual_intra_assay_CVs <- c()

    for (current_day in days) {

      #for (current_rep in replicates) {

      current_set <- fragment_ion_results[fragment_ion_results$fragment_ion==current_ion & fragment_ion_results$day==current_day & fragment_ion_results$sample_group==current_LoMedHi, ]



      # remove rows with a value of NA for calculated_area_ratio
      current_set <- current_set[complete.cases(current_set[ , 'calculated_area_ratio' ]),]



      if (nrow(current_set) <= 1 ){
        percent_CV <- NA
      }
      else{
        percent_CV <- (sd(current_set$calculated_area_ratio))/(mean(current_set$calculated_area_ratio)) * 100
        individual_intra_assay_CVs <- c(individual_intra_assay_CVs, percent_CV)

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





    if (current_LoMedHi=='Lo'){

      CV_results[CV_results$fragment_ion==current_ion, 'low_intra_CV'] <- round(avg_CV, digits=1)
    }
    else if(current_LoMedHi=='Med'){

      CV_results[CV_results$fragment_ion==current_ion, 'med_intra_CV'] <- round(avg_CV, digits=1)

    }
    else if(current_LoMedHi=='Hi'){

      CV_results[CV_results$fragment_ion==current_ion, 'high_intra_CV'] <- round(avg_CV, digits=1)

    }


  } # end current_LoMedHi

} # current_ion

# *** END: intra-assay CV ***



# *** inter-assay CV ***
for (current_ion in ions)  {

  for (current_LoMedHi in LoMedHi) {

    avg_inter_assay_CV <-0
    individual_inter_assay_CVs <- c()

    avg_CV  <- 0

    # for (current_day in days) {

    for (current_rep in replicates) {

      current_set <- fragment_ion_results[fragment_ion_results$fragment_ion==current_ion & fragment_ion_results$replicate==current_rep & fragment_ion_results$sample_group==current_LoMedHi, ]



      # remove rows with a value of NA for calculated_area_ratio
      current_set <- current_set[complete.cases(current_set[ , 'calculated_area_ratio' ]),]




      if (nrow(current_set) <= 1 ){
        percent_CV <- NA
      }
      else{
        percent_CV <- (sd(current_set$calculated_area_ratio))/(mean(current_set$calculated_area_ratio)) * 100
        individual_inter_assay_CVs <- c(individual_inter_assay_CVs, percent_CV)

      }

    } # end current_rep


    # } # end current_day




    if (length(individual_inter_assay_CVs)==0){
      avg_CV <- NA
      count <- 0

    }
    else{
      avg_CV <- mean(individual_inter_assay_CVs, na.rm = TRUE)
    }


    if (current_LoMedHi=='Lo'){

      CV_results[CV_results$fragment_ion==current_ion, 'low_inter_CV'] <- round(avg_CV, digits=1)
    }
    else if(current_LoMedHi=='Med'){

      CV_results[CV_results$fragment_ion==current_ion, 'med_inter_CV'] <- round(avg_CV, digits=1)

    }
    else if(current_LoMedHi=='Hi'){

      CV_results[CV_results$fragment_ion==current_ion, 'high_inter_CV'] <- round(avg_CV, digits=1)

    }


  } # end current_LoMedHi

} # current_ion

# *** END: inter-assay CV ***


# calculate total variability
CV_results[ , 'med_total_CV'] <- round(sqrt((CV_results[ , 'med_intra_CV'])*(CV_results[ , 'med_intra_CV']) + (CV_results[ , 'med_inter_CV'])*(CV_results[ , 'med_inter_CV'])), digits=1)
CV_results[ , 'high_total_CV'] <- round(sqrt((CV_results[ , 'high_intra_CV'])*(CV_results[ , 'high_intra_CV']) + (CV_results[ , 'high_inter_CV'])*(CV_results[ , 'high_inter_CV'])), digits=1)
CV_results[ , 'low_total_CV'] <- round(sqrt((CV_results[ , 'low_intra_CV'])*(CV_results[ , 'low_intra_CV']) + (CV_results[ , 'low_inter_CV'])*(CV_results[ , 'low_inter_CV'])), digits=1)


# determine counts
for (current_ion in ions){

  CV_results[CV_results$fragment_ion==current_ion, 'high_count'] <- nrow(fragment_ion_results[fragment_ion_results$sample_group=='Hi' & fragment_ion_results$fragment_ion==current_ion & !is.na(fragment_ion_results$calculated_area_ratio), ] )
  CV_results[CV_results$fragment_ion==current_ion, 'med_count'] <- nrow(fragment_ion_results[fragment_ion_results$sample_group=='Med' & fragment_ion_results$fragment_ion==current_ion & !is.na(fragment_ion_results$calculated_area_ratio), ] )
  CV_results[CV_results$fragment_ion==current_ion, 'low_count'] <- nrow(fragment_ion_results[fragment_ion_results$sample_group=='Lo' & fragment_ion_results$fragment_ion==current_ion & !is.na(fragment_ion_results$calculated_area_ratio), ] )

}


ions_to_plot <- c()


# determine fragment ions to plot
if (length(ions) <= 4 ) {

  ions_to_plot <- ions


} else {


  results_to_plot <- CV_results[CV_results$fragment_ion!='all' & !is.na(CV_results$low_total_CV) & !is.na(CV_results$med_total_CV) & !is.na(CV_results$high_total_CV) , ]
   
  # new sort to get Top 3 plots
  results_to_plot <- results_to_plot[order(results_to_plot$med_total_CV, results_to_plot$low_total_CV, results_to_plot$high_total_CV), ]
  
  three_lowest_med_total_CV <- head(results_to_plot, 3)


  three_lowest_ions <- as.character(three_lowest_med_total_CV[ , 'fragment_ion'])


  ions_to_plot  <- c(three_lowest_ions, 'all')

}



for (current_plot_ion in ions_to_plot ){


  plot_fragment_ion_results <- fragment_ion_results[!is.na(fragment_ion_results$calculated_area_ratio), ]
  plot_fragment_ion_results <- plot_fragment_ion_results[plot_fragment_ion_results$fragment_ion == current_plot_ion, ]


  plot_days <- sort(unique(plot_fragment_ion_results[ , 'day']))


  # do not make plots for fragment ions with no data
  if (nrow(plot_fragment_ion_results) != 0) {

  # make QC plot for current ion
  plot_QC(plot_fragment_ion_results, input_peptide_sequence, current_plot_ion, plot_days)


  }


}




# ***** output using Panorama *****

# output to web
write.table(CV_results, file = "${tsvout:CV_results}", sep = "\t", qmethod = "double", col.names=NA)
#write.table(fragment_ion_results, file = "${tsvout:fragment_ion_results}", sep = "\t", qmethod = "double", col.names=NA)
#write.table(QC_set, file = "${tsvout:QC_set}", sep = "\t", qmethod = "double", col.names=NA)

# output to files
write.csv(CV_results, file="${fileout:CV_results.csv}")
#write.csv(QC_set, file="${fileout:cpatc_assay_portal_QC_input_set.csv}")
#write.csv(fragment_ion_results, file="${fileout:cpatc_assay_portal_QC_area_ratios.csv}")


#write.table(QC_set, file = "${tsvout:tsvfile}", sep = "\t", qmethod = "double", col.names=NA)

# *********************************




dev.off()
