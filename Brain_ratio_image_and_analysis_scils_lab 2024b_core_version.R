setwd("D:/R/Marilena")

library(SCiLSLabClient)
library(spatstat)
library(data.table)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(dplyr)
library(pROC)
library(ggpubr)
library(car)
library(rstatix)
library(scales)
library(viridis)

# load Scils lab project version 2024b core
example_data <- SCiLSLabOpenLocalSession("D:/R/Marilena/06nedc1.slx", port = 8082)

# Get pixel data
regionTree<-getRegionTree(example_data)
getAttributeDF <- function(regionTree, duplicateValue = NA){
  
  if(!"RegionTree" %in% class(regionTree)){
    stop("RegionTree argument is not of class 'RegionTree'")
  }
  
  allRegions <- flattenRegionTree(regionTree)
  attributes <- new.env()
  for (region in allRegions){
    for (att in 1:nrow(region$attributes)) {
      attribute_name <- region$attributes[att, ]$name
      attribute_level <- region$attributes[att, ]$value
      if (!attribute_name %in% names(attributes)) {
        attributes[[attribute_name]] <- new.env()
      }
      
      if (!attribute_level %in% names(attributes[[attribute_name]])) {
        attributes[[attribute_name]][[attribute_level]] <-
          region$spots$spotId
      } else{
        attributes[[attribute_name]][[attribute_level]] <- union(
          attributes[[attribute_name]][[attribute_level]],
          region$spots$spotId)
      }
    }
  }
  # Create a return data.frame with all spots
  return_df <- allRegions[[1]]$spots
  
  for( attribute_name in names(attributes)){
    attr <- attributes[[attribute_name]]
    
    # Create a temporary data.frame with the spot ids and attribute level for each attribute
    
    tmp <- do.call(rbind.data.frame, lapply(names(attr), function (x){
      data.frame(attr[[x]], x )
    }))
    names(tmp) <- c("spotId", attribute_name)
    
    # Detect and indicate duplicated spots
    duplicates <- duplicated(tmp$spotId) | duplicated(tmp$spotId, fromLast = TRUE)
    tmp[duplicates, 2] <- duplicateValue
    
    tmp <- unique(tmp)
    
    # Merge the temporary data.frame with the return data.frame on the first column
    return_df <- merge(return_df, tmp, by.x = 1, by.y = 1, all.x = TRUE)
  }
  
  return_df
}

# overview of the result:
regionTree <- getRegionTree(example_data)
attribute_df <- getAttributeDF(regionTree)
str(attribute_df)

str(getAttributeDF(regionTree, duplicateValue = "#!duplicated") )

allRegions <- flattenRegionTree(regionTree)

do.call("rbind.data.frame", lapply(allRegions, "[", c(1:2)))

# Generate region name and region unique ID
SpotIndex<-do.call("rbind.data.frame", lapply(allRegions, "[", c(1:2)))

# Retrieve the peak list 
peakLists <- getFeatureLists(example_data)
str(peakLists)

# Get pixel data from feature list named "nga"
allPeaks<-getFeatures(example_data, mode = 'max',name='nga')

# Get the normalization ID for the normalization "Total Ion Count"
normalizations <- getNormalizations(example_data)
ticNorm <- normalizations$uniqueId[normalizations$description == "Total Ion Count"]
rmsNorm<-normalizations$uniqueId[normalizations$description == "Root Mean Square"]

# Get the peak intensity matrix
intensityMatrix <- do.call(
  cbind,
  lapply(allPeaks$id,
         function(x) getFeatureIntensities(
           example_data,
           x,
           mode = 'max',
           regionId = 'Regions',
           normId = ticNorm)$intensity
  )
)


colnames(intensityMatrix) <- makeFeatureNames(allPeaks)

# Merge intensity matrix with our attribute `data.frame` to get a table of all attributes and peak intensities:

spot_ids <- getRegionSpots(example_data)$spotId
intensityDF <- data.frame( spotId = spot_ids, intensityMatrix)

all_data_df <- merge(attribute_df, intensityDF, all.x = TRUE)

# Define raster region as character
all_data_df$raster<-as.character(all_data_df$raster)

# Annotate missing value as 1/5 of the minimum
metabolite_number<-ncol(all_data_df)-6

# Loop metabolite number for pixel data plot
for (i in 1:metabolite_number) {
  metabolite_number1<-i+6
  
  # Annotating missing values as 1/5 of the minimum
  
  all_data_df$Abundance<-all_data_df[,metabolite_number1]
  all_data_df$Abundance[all_data_df$Abundance== 0] = NA
  zhuix<-min(all_data_df$Abundance,na.rm=TRUE)
  all_data_df[,metabolite_number1][all_data_df[,metabolite_number1]==0] <-zhuix/5
}
all_data_df <- all_data_df[,c(1:(metabolite_number+6))]


# Get all subregion ROI coordinates,retrieve ROI ID to get ROI pixel data. Skip this if use whole sections for comparison
coords1 <- getRegionSpots(example_data, regionId = SpotIndex$uniqueId[81])
coords2 <- getRegionSpots(example_data, regionId = SpotIndex$uniqueId[82])

# Assign subregion spotId
subregion1<-coords1$spotId
subregion2<-coords2$spotId

# Filter subregion pixel data from all pixel data set
subregion1_data <- all_data_df %>%
  filter(spotId %in% subregion1)
subregion2_data <- all_data_df %>%
  filter(spotId %in% subregion2)

# Replace raster value with subregion ROI name from spotIndex$name 
subregion1_data$raster<-replace(subregion1_data$raster, subregion1_data$spotId %in% subregion1, SpotIndex$name[81])
subregion2_data$raster<-replace(subregion2_data$raster, subregion2_data$spotId %in% subregion2, SpotIndex$name[82])

# Row bind of subregion ROI dataframe into one subregion ROI dataframe or combine the region in scils lab
ROI_Data<-rbind(subregion1_data,subregion2_data)

# Pixel data file write to designated file path in excel csv file
fwrite(all_data_df, "D:\\R\\Marilena\\112624 Brain AA related whole section pixel data.csv")
fwrite(ROI_Data, "D:\\R\\Marilena\\112624 Brain AA related ROI pixel data.csv")

# Split wh0le section and ROI pixel data set into pixel data alone and spot ID raster name

metabolite_number<-ncol(all_data_df)

AllSectionRatio_Data <-all_data_df[,c(7:metabolite_number)]
AllSectionRatio_Name  <-all_data_df[,c(1:6)]

ROIRatio_Data<-ROI_Data[,c(7:metabolite_number)]
ROIRatio_name<-ROI_Data[,c(1:6)]


# Calculate any two metabolite/feature abundance ratio and create pixel ratio image data for all sections and ROIs
vars0 <- combn(colnames(AllSectionRatio_Data), 2)
vars0 <- cbind(vars0, vars0[2:1,])
done0 <- sapply(seq_len(ncol(vars0)), function(x) AllSectionRatio_Data[, vars0[1, x]] / AllSectionRatio_Data[, vars0[2, x]])
colnames(done0) <- paste(vars0[1, ], vars0[2, ], sep=" / ")
AllSectionRatio_Data<-cbind(AllSectionRatio_Name,done0)


# Ratio calculation for ROI subregion 
vars <- combn(colnames(ROIRatio_Data), 2)
vars <- cbind(vars, vars[2:1,])
done <- sapply(seq_len(ncol(vars)), function(x) ROIRatio_Data[, vars[1, x]] / ROIRatio_Data[, vars[2, x]])
colnames(done) <- paste(vars[1, ], vars[2, ], sep=" / ")
ROIRatio_Data<-cbind(ROIRatio_name,done)



# Ratio Pixel file write to designated file path in excel csv file
fwrite(AllSectionRatio_Data, "D:\\R\\Marilena\\112624 Brain NGA  ratio pixel data.csv")
fwrite(ROIRatio_Data, "D:\\R\\Marilena\\112624 Brain ROI NGA ratio pixel data.csv")

# If pixel data is already available, use the code below to directly import pixel data and proceed image plot and analysis

#  AllSectionRatio_Data <- read.csv("all section pixel data.csv", header = TRUE, stringsAsFactors = TRUE, sep = ",", dec = ".")
#  ROIRatio_Data <- read.csv("ROI oixel data.csv", header = TRUE, stringsAsFactors = TRUE, sep = ",", dec = ".")


# Image and ratio pixel data violin plot and ROC curve for ratios with AUC great than 0.6
pdf(file="112624 Brain all section and ROI ratio image and analysis.PDF", width=7.5, height=7.5)
n=1

# new result list
myresult1<-list()

# obtain the total metabolite number
metabolite_number<-ncol(ROIRatio_Data)-6

for (i in 1:metabolite_number) { 
  metabolite_number1<-i+6
  
  # Red abundnace for each metabolite  
  AllSectionRatio_Data$Abundance<-AllSectionRatio_Data[,metabolite_number1]
  ROIRatio_Data$Abundance<-ROIRatio_Data[,metabolite_number1]
  
  # remove mising value or infinite value
  AllSectionRatio_Data$Abundance[AllSectionRatio_Data$Abundance== Inf] = NA
  
  ROIRatio_Data$Abundance[ROIRatio_Data$Abundance== Inf] = NA
  
  # Create roc curve and calculate P value and AUC
  roc1 <- roc(ROIRatio_Data$raster, ROIRatio_Data$Abundance,na.rm=TRUE)
  
  #  Calculate P value with bonferroni adjust #
  stat.test1 <- ROIRatio_Data %>%
    t_test(Abundance ~ raster) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")
  
  Pvalue1<-as.data.frame(stat.test1)[, "p.adj"]
  pp1<-paste(c("P_KO vs WT",Pvalue1), collapse = "=")
  aa1<-paste(c("AUC_KO vs WT",auc(roc1)), collapse = "=")
  
  # fold change
  x <- tapply(ROIRatio_Data$Abundance, ROIRatio_Data$raster, mean)
  FC1<-x[1]/x[2]
  
  # Create result list for P value AUC value  group name of comparison groups#
  myresult1$Name[i]= colnames(ROIRatio_Data)[i+6]
  myresult1$T_test[i]=Pvalue1
  myresult1$AUC_value[i]=auc(roc1)
  myresult1$Fold_Change[i]=FC1
  myresult1$Group[i]="KO vs WT"
  
  
  # Remove outlier/hotspots
  Quan99<-quantile(AllSectionRatio_Data$Abundance, probs = 0.995,na.rm=TRUE)
  AllSectionRatio_Data$Abundance[AllSectionRatio_Data$Abundance>Quan99] = Quan99
  
  Quan991<-quantile(ROIRatio_Data$Abundance, probs = 0.995,na.rm=TRUE)
  ROIRatio_Data$Abundance[ROIRatio_Data$Abundance>Quan991] = Quan991
  
  #Plot all ratio image, ROI pixel data violin plot and roc curve in one page 
  if(myresult1$AUC_value[i]>0.6){
    
    # Plot all section ratio image with scale bar
    PixelImage0 <- ggplot(data = AllSectionRatio_Data)+
      geom_raster(aes(x = x, y = y, fill = Abundance),interpolate = FALSE,position = "identity",stat="identity") +
      scale_fill_viridis(option="turbo")+theme_void()+coord_equal()
    PixelImage0<-PixelImage0+labs(title=colnames(AllSectionRatio_Data[metabolite_number1]))+
      theme(plot.title=element_text(size=11, hjust=0.5, face="bold", colour="maroon"))
    
    PixelImage0<- PixelImage0+geom_segment(aes(x = min(AllSectionRatio_Data$x),
                                               y = min(AllSectionRatio_Data$y),
                                               xend = min(AllSectionRatio_Data$x)+5000,
                                               yend = min(AllSectionRatio_Data$y)),
                                           linewidth= 0.5, color="maroon")+geom_text(data=AllSectionRatio_Data,size=3,color="maroon",aes(x=min(AllSectionRatio_Data$x)+9500,y=min(AllSectionRatio_Data$y),hjust=1,vjust=0.5,label="5mm")) 
    # Plot ROI ratio image with scale bar 
    
    PixelImage1 <- ggplot(data = ROIRatio_Data)+
      geom_raster(aes(x = x, y = y, fill = Abundance),interpolate = FALSE,position = "identity",stat="identity") +
      scale_fill_viridis(option="turbo")+theme_void()+coord_equal()
    PixelImage1<-PixelImage1+labs(title=colnames(ROIRatio_Data[metabolite_number1]))+
      theme(plot.title=element_text(size=11, hjust=0.5, face="bold", colour="maroon"))
    
    PixelImage1<- PixelImage1+geom_segment(aes(x = min(ROIRatio_Data$x),
                                               y = min(ROIRatio_Data$y),
                                               xend = min(ROIRatio_Data$x)+5000,
                                               yend = min(ROIRatio_Data$y)),
                                           linewidth= 0.5, color="maroon")+geom_text(data=ROIRatio_Data,size=3,color="maroon",aes(x=min(ROIRatio_Data$x)+9500,y=min(ROIRatio_Data$y),hjust=1,vjust=0.5,label="5mm"))
    
    # Restore raster as character for ggplot  
    
    AllSectionRatio_Data$raster<-as.character(AllSectionRatio_Data$raster)    
    ROIRatio_Data$raster<-as.character(ROIRatio_Data$raster)
    
    # Violin and boxplot for ROI pixel ratio data
    P1 <- ggplot(ROIRatio_Data, aes(x=raster, y=Abundance, fill=raster))+geom_violin(trim = FALSE)+theme(plot.title=element_text(size=9, hjust=0.5, face="bold", colour="maroon"))
    P1<-P1 + geom_jitter(shape=16, position=position_jitter(0.1),size=0.05, color="black")
    P1<-P1+labs(x=" ", y="Abundance Ratio")
    P1<- P1+ theme(panel.background = element_blank())
    P1<-P1+theme(axis.line = element_line(colour = "black", size=0.2))
    P1<-P1+theme(text = element_text(size=9),axis.text.x = element_text(angle=15, hjust=0.5))
    
    # plot Roc curve 
    if (require(ggplot2)) {
      P2 <- ggroc(roc1, alpha = 0.5, colour = "red", linetype = 1, size = 2)
      
    }
    P2<-P2 + geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="black", linetype="dashed")
    P2<-P2+labs(title=pp1)+theme(plot.title=element_text(size=7, hjust=0.5,face="bold", color="maroon"))
    P2<-P2+labs(subtitle=aa1)+theme(plot.subtitle=element_text(size=7, hjust=0.5, face="italic", color="blue"))
    
    grid.arrange(PixelImage0,PixelImage1,P1,P2,nrow = 3,layout_matrix = rbind(c(1,1),c(2,2),c(3,4)))
    
    
  }
}
dev.off()

# result list excel output
myresult1_df<-as.data.frame(myresult1)
fwrite(myresult1_df, "D:\\R\\Marilena\\112624 ROI ratio image analysis P FC ROC value.csv")

close(example_data)
