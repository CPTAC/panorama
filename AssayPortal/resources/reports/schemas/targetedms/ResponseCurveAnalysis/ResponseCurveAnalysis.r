# This sample code returns the query data in tab-separated values format, which LabKey then
# renders as HTML. Replace this code with your R script. See the Help tab for more details.

library(Cairo)
library(Rlabkey)
require(plyr)
require(ggplot2)
require(MASS)
require(reshape2)
myplotType <- labkey.url.params$plotType
if (length(myplotType) == 0){
	myplotType <- "linear"
}
mypeptideType <- tolower(labkey.url.params$peptideType)
if (length(mypeptideType) == 0){
	mypeptideType <- "purified"
}

#pivot table 

# Make sure a protein and peptide are selected. Throw an error otherwise.
if(length(labkey.url.params$query.PeptideModifiedSequence) == 0) {.
 	stop("PeptideModifiedSequence filter is missing. Please filter the grid to a single protein and peptide to view the results.")  
 }

# Columns we will get from the 'precursor' table to determine the internal standard type
# for the given peptide
cols <- c("PeptideId/PeptideGroupId/Label",      # Protein name
              "PeptideId/PeptideModifiedSequence", # Modified peptide sequence
              "IsotopeLabelId/Name",                     # name of the isotope label. e.g. 'heavy', 'light'
              "IsotopeLabelId/Standard"                  # true or false
             );


isotopelabels <- labkey.selectRows(
    baseUrl=labkey.url.base,
    folderPath=labkey.url.path,
    schemaName="targetedms",
    queryName="precursor",
    colSelect=paste(cols, sep=","),
    colFilter=makeFilter(c("PeptideId/PeptideModifiedSequence", "EQUAL", labkey.url.params$query.PeptideModifiedSequence))
)


isotopelabels <- unique(isotopelabels)
internal_standards <- isotopelabels[isotopelabels$Standard=='TRUE',]

if(nrow(internal_standards) > 1) {
   stop("Multiple internal standard types found")
}


labkey.data$donotuse <- as.character(labkey.data$donotuse)
labkey.data$donotuse[is.na(labkey.data$donotuse)] <- "FALSE"
labkey.data$donotuse[tolower(labkey.data$donotuse) == "x"] <- "TRUE"
labkey.data <- labkey.data[tolower(labkey.data$donotuse) != "true",];
labkey.data$concentration <- as.numeric(as.character(labkey.data$concentration))
labkey.data$peptideconcentration <- as.numeric(as.character(labkey.data$peptideconcentration))
labkey.data$multiplicationfactor <- as.numeric(as.character(labkey.data$multiplicationfactor))
labkey.data$peptideconcentrationis <- as.numeric(as.character(labkey.data$peptideconcentrationis))
labkey.data$area <- as.numeric(as.character(labkey.data$area))
labkey.data$background <- as.numeric(as.character(labkey.data$background))

if (is.na(labkey.data$concentration[1])){
	labkey.data$concentration <- labkey.data$peptideconcentration * labkey.data$multiplicationfactor;
}


if (is.na(labkey.data$isspike[1])){
	labkey.data$isspike <- labkey.data$peptideconcentrationis;
}

labkey.data$area[is.na(labkey.data$area)] <- 0 
labkey.data$background[is.na(labkey.data$background)] <- 0 

labkey.data$rawArea <- labkey.data$area + labkey.data$background

df <- dcast(labkey.data, protein + peptidemodifiedsequence + precursorcharge + productcharge + fragmention + replicate + concentration + samplegroup + isspike ~ isotopelabel, value.var="rawArea")

 
colnames(df) <-  c("ProteinName", "PeptideModifiedSequence", "PrecursorCharge", "ProductCharge", "FragmentIon", "Replicate", "Concentration", "SampleGroup", "ISSpike", "heavyArea", "lightArea")

df$lightArea[df$heavyArea ==0] <- NA;
df$heavyArea[df$lightArea ==0] <- NA;

if (internal_standards$Label[1]  != "heavy"){
   df$Ratio <- df$heavyArea/df$lightArea;
}
if (internal_standards$Label[1] == "heavy"){
   df$Ratio <- df$lightArea/df$heavyArea;
}
       
df <- df[is.finite(df$Ratio),]

if (internal_standards$Label[1] != "heavy"){
	TSum = ddply(df, .(ProteinName, PeptideModifiedSequence, Concentration, SampleGroup, Replicate, ISSpike), summarize, Ratio = sum(heavyArea, na.rm=TRUE)/sum(lightArea, na.rm=TRUE), heavyArea=sum(heavyArea, na.rm=TRUE), lightArea=sum(lightArea, na.rm=TRUE))
}
if (internal_standards$Label[1] == "heavy"){
	TSum = ddply(df, .(ProteinName, PeptideModifiedSequence, Concentration, SampleGroup, Replicate, ISSpike), summarize, Ratio = sum(lightArea, na.rm=TRUE)/sum(heavyArea, na.rm=TRUE), heavyArea=sum(heavyArea, na.rm=TRUE), lightArea=sum(lightArea, na.rm=TRUE))
}

       
TSum$FragmentIon <- "Sum.tr";
TSum$PrecursorCharge <- "";
TSum$ProductCharge <- "";
       
TSum <- TSum[, names(df)]


df <- rbind(df, TSum)

df$MeasuredConcentration <- df$Ratio * df$ISSpike

       
curveDataIndex <- with(df,  order(ProteinName, PeptideModifiedSequence, PrecursorCharge, FragmentIon, ProductCharge, Concentration, Replicate))
thisPeptide <- df[curveDataIndex,]
thisPeptide$PrecursorCharge <- substr(thisPeptide$PrecursorCharge,1,1)
thisPeptide$ProductCharge <- substr(thisPeptide$ProductCharge,1,1)
uniquePeptide <- unique(thisPeptide$PeptideModifiedSequence)
        

if (length(uniquePeptide) > 1){
  #output error message
}else{
#calculate LOD/LOQ
    thisPeptide$FragmentIon <- paste(thisPeptide$PrecursorCharge, thisPeptide$FragmentIon, thisPeptide$ProductCharge, sep=".") 
	lowConc <- sort(unique(thisPeptide$Concentration[thisPeptide$Concentration>0]))[1]
	usedData <- thisPeptide[,c("FragmentIon", "Concentration", "Replicate", "Ratio")]

   usedData <- usedData[(usedData$Concentration == 0 | usedData$Concentration == lowConc ),   ]
    LODData <- ddply(usedData,  .(FragmentIon, Concentration), 
	summarize, mmean=mean(Ratio, na.rm= TRUE), msd=sd(Ratio, na.rm=TRUE), mqt=qt(0.95,sum(!is.na(Ratio))-1))   

   
	methods <- c("blank+low_conc", "blank_only", "rsd_limit")
	#first method
	LODData[LODData$Concentration == lowConc,'mmean'] <- 0
	LOD1 <- ddply(LODData, .(FragmentIon), summarize, LOD=sum(mmean+msd*mqt))
	LOD1$LOQ <- 3*LOD1$LOD
    names(LOD1) <- c("FragmentIon", paste(methods[1],"LOD", sep="_"), paste(methods[1],"LOQ", sep="_"))
 	usedData <- usedData[(usedData$Concentration == 0) ,] 
    if (dim(usedData)[1] >0){
		LOD2 <- ddply(usedData, .(FragmentIon), 
		summarize, LOD=mean(Ratio, na.rm= TRUE)+ 3* sd(Ratio, na.rm= TRUE), LOQ=mean(Ratio, na.rm= TRUE)+ 10 * sd(Ratio, na.rm= TRUE) )
	
    	names(LOD2) <- c("FragmentIon", paste(methods[2],"LOD", sep="_"), paste(methods[2],"LOQ", sep="_"))
    }
    else{
       LOD2 = LOD1;
       LOD2[,2:3] <- NA
       names(LOD2) <- c("FragmentIon", paste(methods[2],"LOD", sep="_"), paste(methods[2],"LOQ", sep="_"))
    }
	usedData <- thisPeptide[,c("FragmentIon", "Concentration", "Replicate", "Ratio")]
	usedData <- usedData[(usedData$Concentration>0 & !is.na(usedData$Ratio)),  ]
	usedData$estimateConc <- usedData$Concentration * usedData$Ratio;
	LODData <- ddply(usedData, .(FragmentIon, Concentration), summarize, rsd=sd(estimateConc, na.rm= TRUE)/ mean(estimateConc))
	LOD3 <- NULL 
	rsd.max <- 0.15
	temp <- by(LODData, LODData[,c( "FragmentIon")],
			function(x){
				tfit <- try (fit <- lm (log(x$rsd) ~ log (x$Concentration)), silent=TRUE)
				if (!inherits (tfit, "try-error")) {
					# record regression results and LOD/LOQ
        			intercept <- fit$coefficients[1]
     		        slope <- fit$coefficients[2]
			        loq <- exp ( (log (rsd.max) - intercept) / slope )
        			lod <- loq/3
     			} else {
        			lod <- loq <- NA
     			}
     			LOD3 <<- rbind(LOD3, c(as.character(x[1,1]),lod, loq))
			})
	LOD3 <- data.frame(LOD3)
	names(LOD3) <- c("FragmentIon", paste(methods[3],"LOD", sep="_"), paste(methods[3],"LOQ", sep="_"))
	LOD3 <- LOD3[with(LOD3, order(FragmentIon)),]
    LOD3[,2] <- as.numeric(as.character(LOD3[,2]))
    LOD3[,3] <- as.numeric(as.character(LOD3[,3]))

LODTable <- merge(LOD1, LOD2, by="FragmentIon", all=T)
LODTable <- merge(LODTable, LOD3,  by="FragmentIon", all=T)           

result= ddply(df, .(ProteinName, PeptideModifiedSequence, PrecursorCharge, FragmentIon, ProductCharge, Concentration, SampleGroup), 
	summarize, MedianR=median(Ratio, na.rm= TRUE), MinR = min(Ratio, na.rm=TRUE), MaxR=max(Ratio, na.rm=TRUE), CVR=sd(Ratio, na.rm= TRUE)/mean(Ratio, na.rm= TRUE), MedianMeasuredC=median(MeasuredConcentration, na.rm= TRUE), SDMeasuredC=sd(MeasuredConcentration, na.rm= TRUE), MinMeasuredC=min(MeasuredConcentration, na.rm= TRUE), MaxMeasuredC=max(MeasuredConcentration, na.rm= TRUE)); 

result$Median[is.na(result$MedianR)] <- 0;

curveDataIndex <- with( result,  order(ProteinName, PeptideModifiedSequence, PrecursorCharge, FragmentIon, ProductCharge, Concentration))
thisPeptide <- result[curveDataIndex,]
thisPeptide$PrecursorCharge <- substr(thisPeptide$PrecursorCharge,1,1)
thisPeptide$ProductCharge <- substr(thisPeptide$ProductCharge,1,1)
uniquePeptide <- unique(thisPeptide$PeptideModifiedSequence)

}
if (length(uniquePeptide) > 1){
  #output error message
}else{

   
    thisPeptide$FragmentIon <- paste(thisPeptide$PrecursorCharge, thisPeptide$FragmentIon, thisPeptide$ProductCharge, sep=".")
   
    samePeptideLength <- dim(thisPeptide)[1];
    mProtein <- unlist(strsplit(as.character(thisPeptide$ProteinName[1]),"[.]"))[1]
	mTitle <- paste("Analyte: ", mProtein, ".", uniquePeptide[1], "\n", sep="")
    thisPeptide$MaxMeasuredC[thisPeptide$MaxMeasuredC > max(thisPeptide$MedianMeasuredC)*2] <- max(thisPeptide$MedianMeasuredC)
    if ( tolower(myplotType) == "linear") {
	    mxlabel <- "\nTheoretical Concentration (fmol/ug)"
  		mcolor <- "black"
       if (mypeptideType == "crude"){
			mxlabel <- "\nTheoretical Concentration (fmol/ug)\nEstimated from unpurified peptide"
          mcolor <- "red"
    	}
       CairoPNG(filename="${imgout:response_curve_png}", width=800, height=600, bg="white")
		p<- ggplot(data=thisPeptide, aes(x=Concentration, y=MedianMeasuredC, color=FragmentIon)) + geom_errorbar(aes(ymin=MinMeasuredC, 	ymax=MaxMeasuredC), width=.8) + geom_smooth(method=lm, se=FALSE) +geom_point(size=2) + xlab(mxlabel) + ylab("Measured Concentration (fmol/ug") + theme(title=element_text(size=18, colour="black"), axis.text=element_text(size=16), axis.text.x=element_text(colour=mcolor), axis.title=element_text(size=20), axis.title.x=element_text(colour=mcolor), legend.position = c(0.15, 0.85), legend.title = element_text(size=14), legend.text = element_text(size=14))+ labs(title=mTitle)+ scale_colour_discrete(name = "Transition")
		print(p)
		dev.off()
 	}
    if (tolower(myplotType) == "log"){
	    mxlabel <- "\nLog Theoretical Concentration (fmol/ug)"
  		mcolor <- "black"
       if (mypeptideType == "crude"){
			mxlabel <- "\nLog Theoretical Concentration (fmol/ug)\nEstimated from unpurified peptide"
          mcolor <- "red"
    	}
       CairoPNG(filename="${imgout:response_curve_png}", width=800, height=600, bg="white")
        pd <- position_dodge(.05)
		p<- ggplot(data=thisPeptide[thisPeptide$Concentration >0,], aes(x=log(Concentration,10), y=log(MedianMeasuredC,10), color=FragmentIon)) + geom_errorbar(aes(ymin=log(MinMeasuredC,10) , ymax=log(MaxMeasuredC, 10)),position=pd, width=.08)+geom_point(position=pd, size=2) + xlab(mxlabel) + ylab("Log Measured Concentration (fmol/ug") + theme(title=element_text(size=18, colour="black"), axis.text=element_text(size=16), axis.text.x=element_text( colour=mcolor), axis.title=element_text(size=20), axis.title.x=element_text(colour=mcolor), legend.position = c(0.15, 0.85), legend.title = element_text(size=14), legend.text = element_text(size=14))+ labs(title=mTitle)+ scale_colour_discrete(name = "Transition")
		print(p)
		dev.off()
	}
  
   #Calculate LOD/LOD 
   {
      	uniqueT <- unique(thisPeptide$FragmentIon)
        thisPeptide <- thisPeptide[thisPeptide$Concentration >0,]
		thisPeptide <- thisPeptide[with(thisPeptide, order(FragmentIon, Concentration)),]
        thisPeptideR <- thisPeptide[(thisPeptide$MedianMeasuredC != 0 & is.finite(thisPeptide$MedianMeasuredC)),c("FragmentIon", "Concentration", "MedianMeasuredC")]
        fitR <- NULL;                     
 		for (j in 1:length(uniqueT)){
 			thisFragmentIon <- thisPeptideR[(thisPeptideR$FragmentIon == uniqueT[j]),]		
 			x <- thisFragmentIon$Concentration;
 			y <- thisFragmentIon$MedianMeasuredC;
            w <- 1/(y)^2
 			tfit <- try(rlm(x~y, weights=w, method="MM", maxit=1000), silent=TRUE)
            if (!inherits ( tfit, "try-error")){
            	mCoef <- coef(summary(tfit))
                r2 <- (cor(x,y))^2
	           	fitR <- rbind(fitR, c(uniqueT[j], mCoef[2,1], mCoef[1,1], mCoef[2,2], mCoef[1,2],r2))
            }else{
               print("Regression fit failed")
            }
		}
        fitR <- as.data.frame(fitR)
        names(fitR) <-  c("FragmentIon", "Slope", "Intercept", "SlopeStdErr", "InterceptStdErr", "RSquare")  
   }
   for (i in 2: dim(fitR)[2]){
     fitR[,i] <- as.numeric(as.character(fitR[,i]))
   }
   
   #write.table(format(LODTable, digits=3), file = "${tsvout:LODTable.csv}", sep = "\t", qmethod = "double", col.names=NA)   
   write.table(format(fitR, digits=3), file = "${tsvout:fitTable.csv}", sep = "\t", qmethod = "double", col.names=NA)   
 
   write.csv(format(fitR, digits=3), file="${fileout:fitTable.csv}")
   
   mergeTable <- merge(LODTable, fitR, by="FragmentIon", all=T)
   mergeTable[,2:7] <-  (mergeTable[,2:7] - mergeTable[,9])  /(mergeTable[,8])
   mergeTable <- mergeTable[,c(1,2,4,6,3,5,7)]
   write.table(format(mergeTable, digits=3), file = "${tsvout:LODCTable.csv}", sep = "\t", qmethod = "double", col.names=NA)   
   write.csv(format(mergeTable[,1:7], digits=3), file="${fileout:LODCTable.csv}", row.names=F)
}