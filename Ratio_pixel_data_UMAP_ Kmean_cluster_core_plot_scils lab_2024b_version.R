setwd("D:/R/Marilena")


library(SCiLSLabClient)
library(spatstat)
library(data.table)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(dplyr)
library(viridis)
library(magrittr)
library(umap)
library(RColorBrewer)

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

# Get pixel data from feature list named "LPER1"
allPeaks<-getFeatures(example_data, mode = 'max',name='LPER1')

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


# Pixel data file write to designated file path in excel csv file
fwrite(all_data_df, "D:\\R\\Marilena\\ LPE metabolite pixel data.csv")

# ratio calculation
metabolite_number<-ncol(all_data_df)

AllSectionRatio_Data <-all_data_df[,c(7:metabolite_number)]
AllSectionRatio_Name  <-all_data_df[,c(1:6)]


vars0 <- combn(colnames(AllSectionRatio_Data), 2)
vars0 <- cbind(vars0, vars0[2:1,])
done0 <- sapply(seq_len(ncol(vars0)), function(x) AllSectionRatio_Data[, vars0[1, x]] / AllSectionRatio_Data[, vars0[2, x]])
colnames(done0) <- paste(vars0[1, ], vars0[2, ], sep="__ /__")
#
AllSectionRatio_Data<-done0


#  Perform Umap and a k-means clustering using metabolite ratio pixel data. For k-means clustering the number of clusters need to be specified in advance, which usually means to try several settings. In this example, we use 9 clusters.

set.seed(42)
umap_results <- umap::umap(scale(AllSectionRatio_Data))

kmeans_umap <- kmeans(umap_results$layout, centers = 9)



# plot and print the K-mean results
pdf(file="umap9 LPER score and cluster pie plot.PDF", width=5, height=5)

w1<-umap_results$layout %>%
  as.data.frame() 
w1$clustergroup<-kmeans_umap$cluster


P1<-ggplot(w1, aes(x = V1, y = V2, col = as.factor(kmeans_umap$cluster))) +
  geom_point(pch = 20, alpha = 0.5) +
  labs(color = "Cluster") +
  ggtitle("UMAP9 LPER Images")
P1<-P1+scale_color_brewer(palette="Paired")

w2<-kmeans_umap$size %>%
  as.data.frame()
w2$size<-w2[1:9,]
w2$Clu<-c("1","2","3","4","5","6","7","8","9")

pairedcolor<-brewer.pal(n = 9, name = "Paired")

P2<-ggplot(w2,aes(x =" ", y = size, fill=Clu))+geom_bar(width = 1, stat = "identity")+
  coord_polar("y", start = 0)+scale_fill_manual(values = pairedcolor)


grid.arrange(P1,P2,nrow = 2,layout_matrix = rbind(c(1,1),c(2,2)))

dev.off()

#Write Umap score results to excel file
fwrite(w1, "D:\\Imaging\\Marilena\\UMAP9 LPER v1 v2 score export.csv")


# Write UMAP and clustering results back to SCiLS Lab need scils lab pro license
for(i in 1:ncol(umap_results$layout)){
  writeScoreSpotImage( example_data,
                       name = paste("Scores", i),
                       groupname = "LPER UMAP9 Scores",
                       values = data.frame(
                         spotId = spot_ids,
                         value = umap_results$layout[,i])
  )
}

writeLabel(example_data,
           "LPER UMAP9 Clusters",
           labels = data.frame(
             spotId = spot_ids,
             label = kmeans_umap$cluster)
)



close(example_data)
