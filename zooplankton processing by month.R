################################################################################
###################Data processing for depth analyses###########################
########This code starts with Lizzie's zooplankton files in long################
########form and pulls the requested data out by depth, group, #################
########start date, species of interest, and/or lifestage ######################
########as specified at the top of the code. Caution should be##################
########used when looking at a particular species since species#################
########names are matched using the "grep" function that might #################
########return locations for species with similar names. In other###############
########words, it is necessary to inspect the species lists to make#############
########sure that the data summary is providing the desired output##############
################################################################################

######Output consists of these objects:
######1. "zooplankton": A time series by month (i.e. abundances are monthly averages)
######2. "species": a list of the unique species from the processed data
######Code requires the following packages: doBy, xts, zoo


########Before starting two files must be read in to the R environment: 

######## 1. The main data file containing all observations
# Check if main file has been loaded, if not read data into R
if(exists("zoodata")==FALSE){zoodata<-read.csv(file="zoopzeros250key1.csv", stringsAsFactors=F)}       #"c:/Users/gray/Dropbox/Postdoc/R data/zoopzeros250key1.csv")

######## 2. The zooplankton key
#Check if the key has been loaded yet, if not read it in 
if(exists("key")==FALSE){key<-read.csv(file="key.csv")}                            #"c:/Users/gray/Dropbox/Postdoc/R data/key.csv")

#######
#Enter depth of interest in meters (can be <, >, or range, e.g. "0to10")
#Depths available: 0to10, 10to25, 50to100, 100to150, 150to250, lessthan250 (previous 5 combined), "<251"=all depths less than 251m
depth<-"lessthan25"

#start date of data extracted
startdate<-"1955-01-01"

#group of interest ("all", "Copepod", "Cladoceran", "Rotifer")
taxa<-"Copepod"

#Enter species of interest ("all", or e.g. "Epischura baicalensis")
sppec<-"Epischura baicalensis"

#Enter lifestage of interest("all", "adult") #This has not been tested for anything except for "all" or "adult" although many other options are in the key
stage<-"all"

#remove doublecounts? ("yes","no")
dbl<-"no"

######################################################
#################start data processing################
######################################################

allzoo<-zoodata

if(taxa!="all"){
  allzoo<-allzoo[which(allzoo$Group_General==taxa),] }#take records for group of interest

if(stage!="all"){
  allzoo<-allzoo[which(allzoo$Lifestage_Cop==stage),] }#take records for lifestage of interest

#Take data only from species of interest (or take all data if an individual species wasn't specified)
if (sppec!="all"){
  allzoo<-allzoo[grep(sppec, paste(allzoo$Genus,allzoo$Species)),]}

#Take after requested date
allzoo<-allzoo[which(as.Date(allzoo$Date1)>=as.Date(startdate)),]

#create a volume column and total number of individuals column for later calculations
allzoo$Volume<-(as.numeric(allzoo$Nig_Gr)-as.numeric(allzoo$Ver_Gr))*(pi*(0.375/2)^2)*1000 # 0.375 diameter net, multiply by 1000 to convert to litres. Note NAs will be introduced due to some records having no depth information
allzoo<-allzoo[which(allzoo$Volume!="NA"),] #removes records with no volume
allzoo$Tot<-allzoo$Volume*allzoo$Count #Add a column that shows total number of individuals (volume sampled multiplied by Lyubov's densities (#/m3))

#create single depth column by merging Ver_Gr and Nig_Gr
allzoo$Depth<-paste(allzoo$Ver_Gr,allzoo$Nig_Gr, sep="-")

#create separate data sets for each depth category (0-10, 10-25, 50-100, 100-150, 150-250)
allzoo0to10<-allzoo[which(allzoo$Depth=="0-10"),]; allzoo0to10$Depth<-NULL
allzoo10to25<-allzoo[which(allzoo$Depth=="10-25"),]; allzoo10to25$Depth<-NULL
allzoolessthan25<-allzoo[which(allzoo$Depth=="10-25"|allzoo$Depth=="0-10"),]; allzoolessthan25$Depth<-NULL
allzoo50to100<-allzoo[which(allzoo$Depth=="50-100"),]; allzoo50to100$Depth<-NULL
allzoo100to150<-allzoo[which(allzoo$Depth=="100-150"),]; allzoo100to150$Depth<-NULL
allzoo150to250<-allzoo[which(allzoo$Depth=="150-250"),]; allzoo150to250$Depth<-NULL
allzoolessthan250<-allzoo[which(allzoo$Depth=="0-10"|allzoo$Depth=="10-25"|allzoo$Depth=="50-100"|allzoo$Depth=="100-150"|allzoo$Depth=="150-250"),]

#check which depth levels were requested
if (charmatch(">", depth, nomatch=0)>0|charmatch("<", depth, nomatch=0)>0){
  allzoocommondepths<-allzoo[which(as.numeric(allzoo$Nig_Gr)<as.numeric(gsub("\\D", "", depth))),] 
}

if (charmatch(">", depth, nomatch=0)>0|charmatch("<", depth, nomatch=0)==0){
  allzoocommondepths<-get(paste("allzoo", depth, sep=""))}

#remove doublecounts
if(dbl=="yes"){
  KODsnotdouble<-key$KOD[which(key$DoubleCount=="N")] #which codes are not doubles (according to Lizzie)
  allzoocommondepths<-allzoocommondepths[which(allzoocommondepths$KOD%in% KODsnotdouble),]
}

#create a list of species that entered analysis (useful for double-checking in the case where you only one species included)
species<-unique(paste(allzoo$Genus,allzoo$Species))

#create a column of genus +species
allzoocommondepths$genspp<-paste(allzoocommondepths$Genus,allzoocommondepths$Species,sep=" ")

#make all sample dates first of the month (makes it easier to process with the summaryBy function into monthly means)
allzoocommondepths$Date1<-paste(substr(allzoocommondepths$Date1,1,7),"-01",sep="") #convert all samples to 1st of month to make a regular time series

library(doBy)
allzoo3<-summaryBy(Tot+Volume~Date1+genspp, data=allzoocommondepths, FUN=sum, order=TRUE)
allzoo3$Density<-allzoo3$Tot.sum/allzoo3$Volume.sum #Figure out density by adding all individuals and dividing by total volume sampled
allzoo3<-data.frame("Date"=allzoo3$Date1, "Density"=allzoo3$Density, "species"=allzoo3$genspp) #Make a nice dataframe with Date, Density, and KOD
clustall<-reshape(allzoo3,idvar="Date",v.names="Density",timevar="species",direction="wide") #Convert from long form to wide form (species names as column names instead of row names)
names(clustall)<-(gsub("Density.", "", names(clustall)))  #take only numbers from column names. Necessary because text is introduced by the reshape command
row.names(clustall)<-unique(allzoo3$Date)[1:nrow(clustall)]; clustall[,1]<-NULL #Make rownames sample dates

if (ncol(clustall)!=1){
clustall<-clustall[,which(colSums(clustall,na.rm=T)>0)] }#get rid of columns that are all zeroes

zooplankton<-clustall
zooplankton<-data.frame("Date"=as.Date(row.names(zooplankton)), "Density"=zooplankton)

zooneed<-NAfill(zooplankton)
#write.csv(zooplankton, file=paste("zooplankton",depth,".csv"),row.names=F) #write the result to a csv if desired