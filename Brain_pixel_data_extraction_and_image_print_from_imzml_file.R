setwd("D:/R/Marilena")

library(Cardinal)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(data.table)
library(viridis)



# Read imzml data into cardinal using the format of MSProcessedImagingExperiment
pseudo1 <- as(readMSIData("0621Brain.imzML", resolution=5,
                          units="ppm", mass.range=c(80, 1500)),
              "MSProcessedImagingExperiment")

# Get pixel data xyz position and spot ID from one example metabolite image

pixel_d<-pixelData(pseudo1)
pixel_d<-as.data.frame(pixel_d)

 g<-image(pseudo1, mz=191.01973, contrast.enhance="histogram",colorscale=col.map("darkrainbow"))
    myData<-g$facets[[1]][[1]]$values
    myData<-t(myData)
   myDataD<-as.data.frame(myData)
   myDataD[is.na(myDataD)] = -1
  myDataM<-data.matrix(myDataD, rownames.force = NA)
    IndexM<-which(myDataM !=-1)
   mat_pos <- which(myDataM != -1,        
                   arr.ind = TRUE)
  
  #Create data-frame of matrix position
  mat_posD<-as.data.frame(mat_pos)
  
  # subset a dataframe with  x y value
  filteredD <- myDataD[myDataD != -1]
  filteredD <- myDataD[myDataD != -1]

  # Create a starting data frame with similar data structure as scils lab pixel data
 
  m0 <- as.data.frame(matrix(filteredD, ncol = 1, nrow = nrow(mat_posD)))
   m0$spotId<-(1:nrow(mat_posD))
  m0$raster<-0
  m0$x<-mat_posD$col
   m0$y<-mat_posD$row*(-1)
   m0$z<-0
   m0$Date<-0
  m0 = subset(m0, select = c(spotId,raster,x,y,z,Date))
    
   # Extract pixel data from target metabolite list 
   glutamine <- read.csv("metabolite plot.csv", header = TRUE, stringsAsFactors = TRUE, sep = ",", dec = ".")
      metabolites <- unique(glutamine$Compound.Name)
   n=1
   Name<-glutamine$ID
        for (i in 1:length(metabolites)) { 
  t=round(metabolites[i]*8.0e-6, digits=4)
     P1<-image(pseudo1, mz=metabolites[i],plusminus=t,contrast.enhance="histogram")
     testvalue<-P1$facets[[1]][[1]]$values
     myData<-t(testvalue)
    myDataD<-as.data.frame(myData)
    myDataD[is.na(myDataD)] = -1
  filteredD <- myDataD[myDataD != -1]
  
   # Create target list dataframe 
    m <- as.data.frame(matrix(filteredD, ncol = 1, nrow = nrow(mat_posD)))
       colnames(m)<-Name[i]
   m0<-cbind(m0,m)
    n=n+1
  }
 
    # Target list pixel data compilation and output
  
    all_data_df <-m0
    
    # Replace  xy value here with real X Y position for pixel data export
    all_data_df$xy<-paste(all_data_df$x,all_data_df$y)
    pixel_d$y<-(-1)*pixel_d$y
    pixel_d$xy<-paste(pixel_d$x,pixel_d$y)
    NewD1<-merge( pixel_d,all_data_df,by="xy", all.x = TRUE)
     
       NewD1 <- NewD1[, -c(1:5,11:13) ]
       names(NewD1)[names(NewD1) == "X3DPositionX"] <- "x"  
       names(NewD1)[names(NewD1) == "X3DPositionY"] <- "y"  
       names(NewD1)[names(NewD1) == "X3DPositionZ"] <- "z" 
     
       all_data_df<-NewD1
       all_data_df <-all_data_df[,!(names(all_data_df) %in% "xy")]
       all_data_df<-all_data_df[,c(4,5,1,2,3,6,7:ncol(all_data_df))]
       metabolite_number<-ncol(all_data_df)-6
    
    # Annotate missing value with 1/5 minimum 
        for (i in 1:metabolite_number) {
      metabolite_number1<-i+6
      
      all_data_df$Abundance<-all_data_df[,metabolite_number1]
      all_data_df$Abundance[all_data_df$Abundance== 0] = NA
      zhuix<-min(all_data_df$Abundance,na.rm=TRUE)
      all_data_df[,metabolite_number1][all_data_df[,metabolite_number1]==0] <-zhuix/5
      
    }
    all_data_df <- all_data_df[,c(1:(metabolite_number+6))]
  
     # Write pixel data into excel file 
    fwrite(all_data_df, "D:\\R\\Marilena\\120224 brain traget list pixel data.csv")
   
     #  Pixel data image plot with scale bar
        pdf(file="120224 Brain Pixel data image plot from imzml file.PDF", width=7.5, height=2.5)
      metabolite_number<-ncol(all_data_df)-6
  for (i in 1:metabolite_number) {
        metabolite_number1<-i+6
        all_data_df$Abundance<-all_data_df[,metabolite_number1]
        max1<-max(all_data_df$Abundance)
    
    # Remove hot spot outlier then image plot
    Quan99<-quantile(all_data_df$Abundance, probs = 0.9995,na.rm=TRUE)
    all_data_df$Abundance[all_data_df$Abundance>Quan99] = Quan99
        max1<-max(all_data_df$Abundance)
   
    PixelImage0 <- ggplot(data = all_data_df)+
      geom_raster(aes(x = x, y = y, fill = Abundance),interpolate = FALSE,position = "identity",stat="identity") +
      scale_fill_gradientn(colours=c("black","blue","cyan","green","yellow","red"),
                           limits=c(0,max1))+
      theme_void()+coord_equal()

    PixelImage0<-PixelImage0+labs(title=colnames(all_data_df[metabolite_number1]))+
      theme(plot.title=element_text(size=11, hjust=0.5, face="bold", colour="maroon"))

    PixelImage0<- PixelImage0+geom_segment(aes(x = min(all_data_df$x),
                                               y = min(all_data_df$y),
                                               xend = min(all_data_df$x)+5000,
                                               yend = min(all_data_df$y)),
                                           linewidth= 0.5, color="maroon")+geom_text(data=all_data_df,size=3,color="maroon",aes(x=min(all_data_df$x)+9000,y=min(all_data_df$y),hjust=1,vjust=0.5,label="5mm"))
    
    
    
    print(PixelImage0)
  }
  dev.off()
  
  