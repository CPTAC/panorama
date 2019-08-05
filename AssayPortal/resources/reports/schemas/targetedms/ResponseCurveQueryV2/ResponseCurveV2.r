# This sample code returns the query data in tab-separated values format, which LabKey then
# renders as HTML. Replace this code with your R script. See the Help tab for more details.
# write.table(labkey.data, file = "${tsvout:tsvfile}", sep = "\t", qmethod = "double", col.names=NA)
# This sample code returns the query data in tab-separated values format, which LabKey then
# renders as HTML. Replace this code with your R script. See the Help tab for more details.
# write.table(labkey.data, file = "${tsvout:tsvfile}", sep = "\t", qmethod = "double", col.names=NA)f

# renders as HTML. Replace this code with your R script. See the Help tab for more details.

library(Cairo)
library(Rlabkey)
require(plyr)
require(ggplot2)
require(MASS)
require(reshape2)

fitline <- function(xall, yall){
	TotalPoint <- length(xall);
	if (TotalPoint -2 >=1){
		usePoint <- TotalPoint-2;
		myRSquare <- 1;
		while (myRSquare > 0.98 && usePoint >=1){
			x <- xall[usePoint : TotalPoint];
			y <- yall[usePoint : TotalPoint];
			w = 1/y;
			fit <- lm(y~x, weights = w);
			myRSquare <- summary(fit)$r.squared;
			usePoint <- usePoint - 1;
		}
		if (myRSquare <= 0.98){
			if (usePoint == (TotalPoint -3)) {
				usePoint <- usePoint + 1
			}else{
				usePoint <- usePoint + 2;
			}
		}else{
			usePoint <- usePoint + 1;
		}
		bestUsePoint <- usePoint;
		bestEndPoint <- TotalPoint;
		x <- xall[bestUsePoint : bestEndPoint];
		y <- yall[bestUsePoint : bestEndPoint];
		w = 1/y;
		fit <- lm(y~x, weights = w);
		bestRSquare <- summary(fit)$r.squared;
		longest <- bestEndPoint - bestUsePoint + 1;
	}
	#if (myRSquare < 0.98 && usePoint == (TotalPoint-3)){
	if ( TotalPoint-3 >= 1){
		usePoint <- TotalPoint-3;
		myRSquare <- 1;
		while (myRSquare > 0.98 && usePoint >=1 ){
			x <- xall[usePoint : (TotalPoint-1)];
			y <- yall[usePoint : (TotalPoint-1)];
			w = 1/y;
			fit <- lm(y~x, weights = w);
			myRSquare <- summary(fit)$r.squared;
			usePoint <- usePoint - 1;
		}
		if (myRSquare <= 0.98){
			if (usePoint == (TotalPoint -4)) {
				usePoint <- usePoint + 1
			}else{
				usePoint <- usePoint + 2;
			}
		}else{
			usePoint <- usePoint + 1;
		}
		x <- xall[usePoint : (TotalPoint-1)];
		y <- yall[usePoint : (TotalPoint-1)];
		w = 1/y;
		fit <- lm(y~x, weights = w);
		myRSquare <- summary(fit)$r.squared;
		replace <- 0;
		if (bestRSquare < 0.98 && myRSquare > 0.98){
			replace <- 1;
		}else if ( ((bestRSquare > 0.98 && myRSquare > 0.98) || (bestRSquare < 0.98 && myRSquare < 0.98))&& longest < (TotalPoint-usePoint)) {
			replace <- 1;
		}
		if (replace == 1) {
			bestUsePoint <- usePoint;
			bestEndPoint <- TotalPoint-1;
			bestRSquare <- myRSquare;
			longest <- bestEndPoint - bestUsePoint + 1;	
		}
	}
	#if (myRSquare < 0.98 && usePoint == (TotalPoint-4)){
	if (TotalPoint-4 >= 1){	
		usePoint <- TotalPoint-4;
		myRSquare <- 1;
		while (myRSquare > 0.98 && usePoint >=1){
			x <- xall[usePoint : (TotalPoint-2)];
			y <- yall[usePoint : (TotalPoint-2)];
			w = 1/y;
			fit <- lm(y~x, weights = w);
			myRSquare <- summary(fit)$r.squared;			
			usePoint <- usePoint - 1;
		}
		if (myRSquare <= 0.98){
			if (usePoint == (TotalPoint -5)) {
				usePoint <- usePoint + 1
			}else{
				usePoint <- usePoint + 2;
			}
		}else{
			usePoint <- usePoint + 1;
		}
		x <- xall[usePoint : (TotalPoint-2)];
		y <- yall[usePoint : (TotalPoint-2)];
		w = 1/y;
		fit <- lm(y~x, weights = w);
		myRSquare <- summary(fit)$r.squared;
		replace <- 0;
		if (bestRSquare < 0.98 && myRSquare > 0.98){
			replace <- 1;
		}else if (((bestRSquare > 0.98 && myRSquare > 0.98) || (bestRSquare < 0.98 && myRSquare < 0.98)) && longest < (TotalPoint-usePoint-1)) {
			replace <- 1;
		}
		if (replace == 1) {
			bestUsePoint <- usePoint;
			bestEndPoint <- TotalPoint-2;
			bestRSquare <- myRSquare;
			longest <- bestEndPoint - bestUsePoint + 1;	
		}
	}
	if ( TotalPoint <3 ) {
		list();
	}else{
		x <- xall[bestUsePoint : bestEndPoint];
		y <- yall[bestUsePoint : bestEndPoint];
		w <- 1/y;
		fit <- lm(y~x, weights = w); 
		mCoef <- coef(summary(fit));
		mResidual <- c(TotalPoint);
		for (i in 1: TotalPoint){
			mResidual[i] <- (yall[i]-(xall[i]*mCoef[2,1]+mCoef[1,1]))/(yall[i]+(xall[i]*mCoef[2,1]+mCoef[1,1]));
		}
		usedPoints = bestEndPoint - bestUsePoint + 1;
		Rsquare <-summary(fit)$r.squared;
		list(mCoef=mCoef, mResidual=mResidual, Rsquare=Rsquare, usedPoints = usedPoints)
	}
}

       

#pivot table 

myplotType <- labkey.url.params$plotType
if (length(myplotType) == 0){
	myplotType <- "linear"
}

mypeptideType <- tolower(labkey.url.params$peptideType)
if (length(mypeptideType) == 0){
	mypeptideType <- "purified"
    #mypeptideType <- "crude"
}
       
# Make sure a protein and peptide are selected. Throw an error otherwise.
if( length(labkey.url.params$query.PeptideModifiedSequence) == 0) {.
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
    #colFilter=makeFilter(c("PeptideId/PeptideGroupId/Label", "EQUAL", labkey.url.params$query.Protein),c("PeptideId/PeptideModifiedSequence", "EQUAL", labkey.url.params$query.PeptideModifiedSequence))
    colFilter=makeFilter(c("PeptideId/PeptideModifiedSequence", "EQUAL", labkey.url.params$query.PeptideModifiedSequence))

)

isotopelabels <- unique(isotopelabels)
internal_standards <- isotopelabels[isotopelabels$Standard=='TRUE',]
if(nrow(internal_standards) > 1) {
   stop("Multiple internal standard types found")
}
       
if (sum(names(labkey.data$name) == "donotuse") > 0){ 
	labkey.data$donotuse <- as.character(labkey.data$donotuse)
	labkey.data$donotuse[is.na(labkey.data$donotuse)] <- "FALSE"
	labkey.data$donotuse[tolower(labkey.data$donotuse) == "x"] <- "TRUE"
	labkey.data <- labkey.data[tolower(labkey.data$donotuse) != "true",];
}

colnames(labkey.data)[colnames(labkey.data)=="analyteconcentration"] <- "concentration"
colnames(labkey.data)[colnames(labkey.data)=="internalstandardconcentration"] <- "isspike"
colnames(labkey.data)[colnames(labkey.data)=="concentrationmultiplier"] <- "multiplicationfactor"
colnames(labkey.data)[colnames(labkey.data)=="replicatenumber"] <- "replicate"

labkey.data$concentration <- as.numeric(as.character(labkey.data$concentration))
#labkey.data$peptideconcentration <- as.numeric(as.character(labkey.data$peptideconcentration))
labkey.data$multiplicationfactor <- as.numeric(as.character(labkey.data$multiplicationfactor))
#labkey.data$peptideconcentrationis <- as.numeric(as.character(labkey.data$peptideconcentrationis))
labkey.data$area <- as.numeric(as.character(labkey.data$area))
labkey.data$background <- as.numeric(as.character(labkey.data$background))

if (!is.na(labkey.data$multiplicationfactor[1])){
	labkey.data$concentration <- labkey.data$concentration * labkey.data$multiplicationfactor;
}

if(length(unique(labkey.data$concentration)) <2) {
   stop("More than one concentration level needed.")
}
       

#if (is.na(labkey.data$isspike[1])){
#	labkey.data$isspike <- labkey.data$peptideconcentrationis;
#}


labkey.data$area[is.na(labkey.data$area)] <- 0 
if (length(grep("background", names(labkey.data))) >0){
	labkey.data$background[is.na(labkey.data$background)] <- 0 
 	labkey.data$rawArea <- labkey.data$area + labkey.data$background 
}else{
	labkey.data$rawArea <- labkey.data$area 
}

df <- dcast(labkey.data, protein + peptidemodifiedsequence + precursorcharge + productcharge + fragmention + replicate + concentration + samplegroup + isspike ~ isotopelabel, value.var="rawArea")

colnames(df) <-  c("ProteinName", "PeptideModifiedSequence", "PrecursorCharge", "ProductCharge", "FragmentIon", "Replicate", "Concentration", "SampleGroup", "ISSpike", "heavyArea", "lightArea")

#write.table(df, file = "${tsvout:tsvfile}", sep = "\t", qmethod = "double", col.names=NA)


df$lightArea[is.na(df$heavyArea)] = NA;
df$heavyArea[is.na(df$lightArea)] = NA;
df$heavyArea[df$lightArea==0] <- NA;

if (internal_standards$Label != "heavy"){
   df$HLRatio <- df$heavyArea/df$lightArea;
}
if (internal_standards$Label == "heavy"){
   df$HLRatio <- df$lightArea/df$heavyArea;
}

df <- df[is.finite(df$HLRatio),]

if (internal_standards$Label[1] != "heavy"){
	TSum = ddply(df, .(ProteinName, PeptideModifiedSequence, PrecursorCharge, Concentration, SampleGroup, Replicate), summarize, HLRatio = sum(heavyArea, na.rm=TRUE)/sum(lightArea, na.rm=TRUE))
}


if (internal_standards$Label[1] == "heavy"){
	TSum = ddply(df, .(ProteinName, PeptideModifiedSequence, PrecursorCharge, Concentration, SampleGroup, Replicate), summarize, HLRatio = sum(lightArea, na.rm=TRUE)/sum(heavyArea, na.rm=TRUE))
}


TSum$FragmentIon <- "SUM";
TSum$ProductCharge <- "";

resultSum = ddply(TSum, .(ProteinName, PeptideModifiedSequence, PrecursorCharge, FragmentIon, ProductCharge, Concentration, SampleGroup), 
	summarize, Median=median(HLRatio, na.rm= TRUE), Min = min(HLRatio, na.rm=TRUE), Max=max(HLRatio, na.rm=TRUE), CV=sd(HLRatio, na.rm= TRUE)/mean(HLRatio, na.rm= TRUE), LOD=mean(HLRatio, na.rm= TRUE)+ 3* sd(HLRatio, na.rm= TRUE));       

if (internal_standards$Label[1] != "heavy"){
	temp = ddply(df, .(ProteinName, PeptideModifiedSequence, PrecursorCharge, FragmentIon, ProductCharge), summarize, medianArea=median(heavyArea, na.rm=TRUE))
}

if (internal_standards$Label[1] == "heavy"){
	temp = ddply(df, .(ProteinName, PeptideModifiedSequence, PrecursorCharge, FragmentIon, ProductCharge), summarize, medianArea=median(lightArea, na.rm=TRUE))
}

orderT <- with(temp,  order(ProteinName, PeptideModifiedSequence, -medianArea))

df <- merge(df, temp[orderT[1:3], ])


resultT= ddply(df, .(ProteinName, PeptideModifiedSequence, PrecursorCharge, FragmentIon, ProductCharge, Concentration, SampleGroup), 
	summarize, Median=median(HLRatio, na.rm= TRUE), Min = min(HLRatio, na.rm=TRUE), Max=max(HLRatio, na.rm=TRUE), CV=sd(HLRatio, na.rm= TRUE)/mean(HLRatio, na.rm= TRUE), LOD=mean(HLRatio, na.rm= TRUE)+ 3* sd(HLRatio, na.rm= TRUE));       

result <- rbind(resultT, resultSum)

result$Median[is.na(result$Median)] <- 0;

curveDataIndex <- with( result,  order(ProteinName, PeptideModifiedSequence, PrecursorCharge, FragmentIon, ProductCharge, Concentration))
thisPeptide <- result[curveDataIndex,]
thisPeptide$PrecursorCharge <- substr(thisPeptide$PrecursorCharge,1,1)
thisPeptide$ProductCharge <- substr(thisPeptide$ProductCharge,1,1)
uniquePeptide <- unique(thisPeptide$PeptideModifiedSequence)


if (length(uniquePeptide) > 1){
  #output error message
}else{

    thisPeptide$FragmentIon <- paste(thisPeptide$PrecursorCharge, thisPeptide$FragmentIon, thisPeptide$ProductCharge, sep=".")
	samePeptideLength <- dim(thisPeptide)[1];
    mProtein <- unlist(strsplit(as.character(thisPeptide$ProteinName[1]),"[.]"))[1]
	mTitle <- paste("Analyte: ", mProtein, ".", uniquePeptide[1], "\n", sep="")
    thisPeptide$Max[thisPeptide$Max > max(thisPeptide$Median)*2] <- max(thisPeptide$Median)
    if ( tolower(myplotType) == "linear") {
	   mxlabel <- "\nTheoretical Concentration (fmol/uL)"
	   mcolor <- "black"
       if (tolower(mypeptideType) == "crude"){
			mxlabel <- "\nTheoretical Concentration (fmol/uL)\nEstimated from unpurified peptide"
          mcolor <- "red"
    	}       
  		CairoPNG(filename="${imgout:response_curve_png}", width=800, height=600, bg="white")
		p<- ggplot(data=thisPeptide, aes(x=Concentration, y=Median, color=FragmentIon)) + geom_errorbar(aes(ymin=Min, ymax=Max), width=.8) + geom_smooth(method=lm, se=FALSE) +geom_point(size=2) + xlab(mxlabel) + ylab("Peak Area Ratio") + theme(title=element_text(size=18, colour="black"), axis.text=element_text(size=16), axis.text.x=element_text(colour=mcolor), axis.title=element_text(size=20) , axis.title.x=element_text(colour=mcolor), legend.position = c(0.15, 0.85), legend.title = element_text(size=14), legend.text = element_text(size=14))+ labs(title=mTitle)+ scale_colour_discrete(name = "Transition")
		print(p)
		dev.off()
    }
    if (tolower(myplotType) == "log"){
       mxlabel <- "\nLog Theoretical Concentration (fmol/uL)"
	   mcolor <- "black"
       if (tolower(mypeptideType) == "crude"){
			mxlabel <- "\nLog Theoretical Concentration (fmol/uL)\nEstimated from unpurified peptide"
          mcolor <- "red"
    	}       
       CairoPNG(filename="${imgout:response_curve_png}", width=800, height=600, bg="white")
       pd <- position_dodge(.05)
		p<- ggplot(data=thisPeptide[thisPeptide$Concentration >0,], aes(x=log(Concentration,10), y=log(Median,10), color=FragmentIon)) + geom_errorbar(aes(ymin=log(Min,10), ymax=log(Max,10)), position=pd, width=.08) +geom_point(position=pd, size=2) + xlab(mxlabel) + ylab("Log Peak Area Ratio") + theme(title=element_text(size=18, colour="black"), axis.text=element_text(size=16),   axis.text.x=element_text(colour=mcolor), axis.title=element_text(size=20) , axis.title.x=element_text( colour=mcolor), legend.position = c(0.2, 0.85), legend.title = element_text(size=14), legend.text = element_text(size=14))+ labs(title=mTitle)+ scale_colour_discrete(name = "Transition")    
		print(p)
		dev.off()
	}
    if ( tolower(myplotType) == "residual"){
       mxlabel <- "\nLog Theoretical Concentration (fmol/uL)"
       mcolor <- "black"
       if (tolower(mypeptideType) == "crude"){
			mxlabel <- "\nLog Theoretical Concentration (fmol/uL)\nEstimated from unpurified peptide"
          mcolor <- "red"
    	}       
	    uniqueT <- unique(thisPeptide$FragmentIon)
		thisPeptide <- thisPeptide[with(thisPeptide, order(FragmentIon, Concentration)),]
		thisPeptideR <- thisPeptide[(thisPeptide$Median != 0 & is.finite(thisPeptide$Median)),  ]      
 		for (j in 1:length(uniqueT)){
 			thisFragmentIon <- thisPeptideR[(thisPeptideR$FragmentIon == uniqueT[j]),]		
 			x <- thisFragmentIon$Concentration;
 			y <- thisFragmentIon$Median;
            if (length(x) > 2){
 	 			Lvalue <- fitline(x,y)       
 			    thisPeptideR$Residual[(thisPeptideR$FragmentIon == uniqueT[j])] <- Lvalue$mResidual;
            }else{
                thisPeptideR$Residual[(thisPeptideR$FragmentIon == uniqueT[j])] <- NA
            }
		}
		CairoPNG(filename="${imgout:response_curve_png}", width=800, height=600, bg="white")
 		p<- ggplot(data=thisPeptideR, aes(x=log(Concentration,10), y=Residual, color=FragmentIon)) +geom_point(size=2) + xlab(mxlabel) + ylab("Residual") + theme(title=element_text(size=18, colour="black"), axis.text=element_text(size=16) , axis.text.x=element_text(colour=mcolor), axis.title=element_text(size=20) , axis.title.x=element_text(colour=mcolor), legend.title = element_text(size=14), legend.text = element_text(size=14))+ labs(title=mTitle)+ scale_colour_discrete(name = "Transition")
 		print(p)
	 	dev.off()
    }
}



#CairoPNG(filename="${imgout:response_curve_png}", width=800, height=600, bg="white")

write.table(thisPeptide, file = "${tsvout:tsvfile}", sep = "\t", qmethod = "double", col.names=NA)