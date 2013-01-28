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

########Output consists of these objects:
######1. "zooplankton": monthly time series or series by sample date depending on what is specified in "seasons" below
######2. "winter","spring","summer","fall": If seasonal time series is requested, there will be four objects "winter", "spring", "summer", "fall" that summarize abundance over time by seasons for the species requested
######3. "species": list of species found in the data
######4. "QC": 1000 randomly chosen cells from the dataframe for comparison with Lyubov's original data
######5. "anydoubles": For each species, were any codes doubles? Uses the anyDuplicated function to look for duplicated columns in dataframe. If NULL or all zeroes, then no duplicates detected. 
######6. "quarterly": Average abundance by quarter for all species (not constructed if data is summarized by sample date)
########Code requires the following packages: doBy, xts, zoo


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
interest<-c("150to250")

#start date of data extracted
startdate=c("1955-01-01")

#monthly, summarize by seasons=quarters (options= "seasons", "monthly", "sampledate")
seasons<-c("sampledate")

#group of interest ("all", "Copepod", "Cladoceran", "Rotifer")
taxa<-c("all")

#Enter species of interest ("all", or e.g. "Epischura baicalensis")
sppec<-c("all")

#Enter lifestage of interest("all", "adult") #This has not been tested for anything except for "all" or "adult" although many other options are in the key
stage=c("all")

#min number of observations per species by month (approx 720 months in dataset) species with fewer observations discarded)
minn=1

######################################################
#################start data processing################
######################################################

allzoo<-zoodata

if(taxa!="all"){
allzoo<-allzoo[which(allzoo$Group_General==taxa),] }#take records for group of interest

#create a volume column and total number of individuals column for later calculations
allzoo$Volume<-(as.numeric(allzoo$Nig_Gr)-as.numeric(allzoo$Ver_Gr))*(pi*(0.375/2)^2)*1000 # 0.375 diameter net, multiply by 1000 to convert to litres. Note NAs will be introduced due to some records having no depth information
allzoo<-allzoo[which(allzoo$Volume!="NA"),] #removes records with no volume
allzoo$Tot<-allzoo$Volume*allzoo$Count #Add a column that shows total number of individuals (volume sampled multiplied by Lyubov's densities (#/m3))
#create single depth column by merging Ver_Gr and Nig_Gr
allzoo$Depth<-paste(allzoo$Ver_Gr,allzoo$Nig_Gr, sep="-")

#create separate data sets for each depth category (0-10, 10-25, 50-100, 100-150, 150-250)
allzoo0to10<-allzoo[which(allzoo$Depth=="0-10"),]; allzoo0to10$Depth<-NULL
allzoo10to25<-allzoo[which(allzoo$Depth=="10-25"),]; allzoo10to25$Depth<-NULL
allzoo50to100<-allzoo[which(allzoo$Depth=="50-100"),]; allzoo50to100$Depth<-NULL
allzoo100to150<-allzoo[which(allzoo$Depth=="100-150"),]; allzoo100to150$Depth<-NULL
allzoo150to250<-allzoo[which(allzoo$Depth=="150-250"),]; allzoo150to250$Depth<-NULL
allzoolessthan250<-allzoo[which(allzoo$Depth=="0-10"|allzoo$Depth=="10-25"|allzoo$Depth=="50-100"|allzoo$Depth=="100-150"|allzoo$Depth=="150-250"),]


#check which depth levels were requested
if (charmatch(">", interest, nomatch=0)>0|charmatch("<", interest, nomatch=0)>0){
  allzoocommondepths<-allzoo[which(as.numeric(allzoo$Nig_Gr)<as.numeric(gsub("\\D", "", interest))),] 
}

if (charmatch(">", interest, nomatch=0)>0|charmatch("<", interest, nomatch=0)==0){
  allzoocommondepths<-get(paste("allzoo", interest, sep=""))}


#summarize by month
if (seasons!="sampledate"){allzoocommondepths$Date1<-paste(substr(allzoocommondepths$Date1,1,7),"-01",sep="")} #convert all samples to 1st of month to make a regular time series

library(doBy)
allzoo3<-summaryBy(Tot+Volume~Date1+KOD, data=allzoocommondepths, FUN=sum, order=TRUE)
allzoo3$Density<-allzoo3$Tot.sum/allzoo3$Volume.sum #Figure out density by adding all individuals and dividing by total volume sampled
allzoo3<-data.frame("Date"=allzoo3$Date1, "Density"=allzoo3$Density, "KOD"=allzoo3$KOD) #Make a nice dataframe with Date, Density, and KOD

clustall<-reshape(allzoo3,idvar="Date",v.names="Density",timevar="KOD",direction="wide") #Convert from long form to wide form (species names as column names instead of row names)
clustall<-clustall[1:(nrow(clustall)-1),] #remove the NA row
names(clustall)<-as.numeric(gsub("\\D", "", names(clustall)))  #take only numbers from column names. Necessary because text is introduced by the reshape command
row.names(clustall)<-unique(allzoo3$Date)[1:nrow(clustall)]; clustall[,1]<-NULL #Make rownames sample dates
clustall<-clustall[,which(colSums(clustall,na.rm=T)>0)] #get rid of columns that are all zeroes

#Figure out which KODs represent individual species (i.e. species may be represented by multiple KODs)
matchlocations<-which(names(clustall) %in% key$KOD)
result<-data.frame("Genus"=key$Genus[matchlocations], "Species"=key$Species[matchlocations], "Lifestage"=key$Lifestage_Gen[matchlocations], "Copepod stage"=key$Lifestage_Cop[matchlocations])
wholesplist<-paste(key$Genus,key$Species, sep=" ")
unspecies<-unique(wholesplist)

#This section combines data from different KODs if they represent the same species (but it eliminates double-counts as classified in Lizzie's key)
anydoubles<-c() #For quality control, see in loop below
if(stage=="all"){
  #combine adults and juveniles (combine columns with different KODs)
  for (i in 1:length(unspecies)){
    group<-key$KOD[which(wholesplist %in% unspecies[i])]
    KODsnotdouble<-key$KOD[which(key$DoubleCount=="N")] #which codes are not doubles (according to Lizzie)
    group<-group[which(group%in%KODsnotdouble==TRUE)] #include only codes that are not doubles
    
    #looks for doubles that could indicate KODs that refer to the same count
    if (length(which(duplicated(t(as.data.frame(clustall[,which(names(clustall)%in% group)])),MARGIN=1)==TRUE))>0){
      anydoubles[i]<-names(clustall)[which(duplicated(t(as.data.frame(clustall[,which(names(clustall)%in% group)])),MARGIN=1)==TRUE)]}
    
    if(is.data.frame(clustall[,which(names(clustall)%in% group)])==TRUE){
      if (i==1){sumdataframe<-data.frame(rowMeans(clustall[,which(names(clustall)%in% group)],na.rm=T)); names(sumdataframe)[1]<-as.character(unspecies[1])}
      if(i>1){sumdataframe[,i]<-rowMeans(clustall[,which(names(clustall)%in% group)],na.rm=T); names(sumdataframe)[i]<-as.character(unspecies[i])}}
    if(is.data.frame(clustall[,which(names(clustall)%in% group)])==FALSE){
      if (i==1){sumdataframe<-data.frame(clustall[,which(names(clustall)%in% group)]); names(sumdataframe)[1]<-as.character(unspecies[1])}
      if (i>1){sumdataframe[,i]<-clustall[,which(names(clustall)%in% group)]; names(sumdataframe)[i]<-as.character(unspecies[i])}
    }}
  sumdataframe<-sumdataframe[,which(colSums(sumdataframe, na.rm=T)!=0)] #REMOVE species with no observations
  #sumdataframe<-sumdataframe[,setdiff(1:ncol(sumdataframe), grep("unknown", names(sumdataframe)))] #get rid of species names with "unknown"
  sumdataframe[,which(names(sumdataframe)=="Golomyanka ")]<-NULL
  sumdataframe[,which(names(sumdataframe)=="Macrohectopus branickii")]<-NULL
}

if(stage!="all"){
  #Get data for a select species and separate by life stage
  for (i in 1:length(unspecies)){
    group<-key$KOD[which(wholesplist %in% unspecies[i])]
    
    #looks for doubles that could indicate KODs that refer to the same count
    if (length(which(duplicated(t(as.data.frame(clustall[,which(names(clustall)%in% group)])),MARGIN=1)==TRUE))>0){
      anydoubles[i]<-names(clustall)[which(duplicated(t(as.data.frame(clustall[,which(names(clustall)%in% group)])),MARGIN=1)==TRUE)]}
    
    KODsnotdouble<-key$KOD[which(key$DoubleCount=="N")] #which codes are not doubles (according to Lizzie's key)
    group<-group[which(group%in%KODsnotdouble==TRUE)] #include only codes that are not doubles
    lifestage<-as.vector(key$Lifestage_Cop[which(wholesplist %in% unspecies[i])])
    newgroup<-grep(stage, lifestage)
    if(is.data.frame(clustall[,which(names(clustall)%in% newgroup)])==TRUE){
      if (i==1){sumdataframe<-data.frame(rowMeans(clustall[,which(names(clustall)%in% newgroup)],na.rm=T)); names(sumdataframe)[1]<-as.character(unspecies[1])}
      if(i>1){sumdataframe[,i]<-rowMeans(clustall[,which(names(clustall)%in% newgroup)],na.rm=T); names(sumdataframe)[i]<-as.character(unspecies[i])}}
    if(is.data.frame(clustall[,which(names(clustall)%in% newgroup)])==FALSE){
      if (i==1){sumdataframe<-data.frame(clustall[,which(names(clustall)%in% group)]); names(sumdataframe)[1]<-as.character(unspecies[1])}
      if (i>1){sumdataframe[,i]<-clustall[,which(names(clustall)%in% group)]; names(sumdataframe)[i]<-as.character(unspecies[i])}
    }}}
row.names(sumdataframe)<-row.names(clustall)

#Take data only from species of interest (or take all data if an individual species wasn't specified)
if (sppec=="all"){series<-sumdataframe
                  series$Date<-as.Date(row.names(sumdataframe))}
if (sppec!="all"){
  series<-data.frame("Count"=sumdataframe[,grep(sppec, names(sumdataframe))], "Date"=as.Date(row.names(sumdataframe)))}
newnb<-series

if(seasons!="sampledate"){
#Identify months with no data and put an NA where data should be (analagous to joining a table with all dates to identify missing samples)
years<-substr(newnb$Date,1,4)
start = as.Date(newnb$Date[which.min(newnb$Date)])
full <- seq(start, by="1 month", length=length(unique(years))*12)
allsp<-data.frame("Date"=full)
for (i in 1:ncol(newnb)){
  ebaicalensis<-data.frame("Count"=with(newnb, newnb[,i][match(full, as.Date(newnb$Date))]))
  allsp<-cbind(allsp, ebaicalensis)}
allsp<-allsp[,2:ncol(allsp)]; names(allsp)<-names(newnb); allsp$Date<-full

nb<-allsp[which(allsp$Date>=startdate),]} #take data only after startdate

if(seasons=="sampledate"){
  nb<-newnb[which(newnb$Date>=startdate),]
}

zooplankton<-nb; row.names(zooplankton)<-zooplankton$Date; zooplankton$Date<-NULL; if(sppec=="all"){zooplankton<-zooplankton[,which(colSums(zooplankton,na.rm=T)>0)]}

#This removes species with observations below the "minn" threshold set at the top of the code
if(sppec=="all"){
belowmin<-c()
for (i in 1:ncol(zooplankton)){
  belowmin[i]<-length(which(zooplankton[,i]>0))}
zooplankton<-zooplankton[,which(belowmin>minn)]
}

if (sppec!="all"){names(zooplankton)<-sppec} #rename dataframe heading for single species
##########################
#summarize data by seasons?
##########################
if(sppec=="all"){
if (seasons=="seasons"){
  for(i in 1:(ncol(zooplankton))){
    library(xts)
    qtnb<-xts(x=zooplankton[,i], order.by=as.Date(row.names(zooplankton)), frequency=12)
    locations<-endpoints(qtnb, "quarters")-1
    locations<-locations[which(locations>0)]
    if(i==1){qtrnb<-data.frame("test"=rep(0,times=(length(locations)-1)))}
    qtrnb<-as.data.frame(period.apply(qtnb, INDEX=locations, FUN=function(x) mean(x,na.rm=T)))
    
    names(qtrnb)<-names(zooplankton)[i]
    if (i==1){quarterly<-data.frame(qtrnb)}
    if (i>1){quarterly[,i]<-qtrnb}
    
    f<-c(1,2,3,4); h<-rep(f, times=1000)
    startmnth<-as.numeric(substr(nb$Date[locations[2]],6,7))
    if (startmnth==12|startmnth==1|startmnth==2){startqtr<-1}
    if (startmnth==3|startmnth==4|startmnth==5){startqtr<-2}
    if (startmnth==6|startmnth==7|startmnth==8){startqtr<-3}
    if (startmnth==9|startmnth==10|startmnth==11){startqtr<-4}
    qtrnb$quarter<-h[(min(which(h==startqtr))):(nrow(qtrnb)+min(which(h==startqtr))-1)]
    
    if(i==1){winter<-data.frame("test"=rep(0,times=(length(which(qtrnb$quarter==1)))))}
    if(i==1){spring<-data.frame("test"=rep(0,times=(length(which(qtrnb$quarter==2)))))}
    if(i==1){summer<-data.frame("test"=rep(0,times=(length(which(qtrnb$quarter==3)))))}
    if(i==1){fall<-data.frame("test"=rep(0,times=(length(which(qtrnb$quarter==4)))))}
    
    winter[,i]<-qtrnb[which(qtrnb$quarter==1),1]; winter$quarter<-NULL
    spring[,i]<-qtrnb[which(qtrnb$quarter==2),1]; spring$quarter<-NULL
    summer[,i]<-qtrnb[which(qtrnb$quarter==3),1]; summer$quarter<-NULL
    fall[,i]<-qtrnb[which(qtrnb$quarter==4),1]; fall$quarter<-NULL
  }}

#Remove species with no observations and add species names as column headings
nonzerolocs<-which(colSums(winter,na.rm=T)>0)
winter<-winter[,nonzerolocs]
names(winter)<-names(zooplankton)[nonzerolocs]; winter2<-winter[,order(names(winter))]
winter<-as.data.frame(winter2); names(winter)<-names(winter2)
nonzerolocs<-which(colSums(spring,na.rm=T)>0)
spring<-spring[,nonzerolocs]; 
names(spring)<-names(zooplankton)[nonzerolocs]; spring2<-spring[,order(names(spring))]
spring<-as.data.frame(spring2); names(spring)<-names(spring2)
nonzerolocs<-which(colSums(summer,na.rm=T)>0)
summer<-summer[,nonzerolocs]
names(summer)<-names(zooplankton)[nonzerolocs]; summer2<-summer[,order(names(summer))] 
summer<-as.data.frame(summer2); names(summer)<-names(summer2)
nonzerolocs<-which(colSums(fall,na.rm=T)>0)
fall<-fall[,nonzerolocs]
names(fall)<-names(zooplankton)[nonzerolocs]; fall2<-fall[,order(names(fall))]
fall<-as.data.frame(fall2); names(fall)<-names(fall2)


#species lists
if (seasons=="seasons"){species<-names(winter)}
if (seasons=="monthly"){species<-names(zooplankton)}


#Choose 50 random cells for quality control
if (seasons=="seasons"){
  spp<-list(); num<-list(); qrt<-list()
  
  for (i in 1:1000){
    rowloc<-round(runif(1,1, nrow(winter)),0)
    colloc<-round(runif(1,1,ncol(winter)),0)
    spp[[i]]<-(names(winter)[colloc])
    num[[i]]<-winter[rowloc,colloc]
    qrt[[i]]<-row.names(winter)[rowloc]
  }
  
  QC<-data.frame("Species"=unlist(spp), "Abundance"=unlist(num), "Quarter"=unlist(qrt))
  QC<-orderBy(Abundance~., data=QC)  
}
}

if (seasons=="monthly"){
  
  spp<-list(); num<-list(); date<-list()
  
  for (i in 1:1000){
    rowloc<-round(runif(1,1, nrow(zooplankton)),0)
    colloc<-round(runif(1,1,ncol(zooplankton)),0)
    spp[[i]]<-(names(zooplankton)[colloc])
    num[[i]]<-zooplankton[rowloc,colloc]
    date[[i]]<-row.names(zooplankton)[rowloc]
  }
  QC<-data.frame("Species"=unlist(spp), "Abundance"=unlist(num), "Date"=unlist(date))
  QC<-orderBy(Abundance~., data=QC) }


write.csv(zooplankton, file="zooplankton150to250.csv")