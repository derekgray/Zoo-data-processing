
# Read data into R, convert date column to "Date" type rather than string, move "Date" column to be first
data<-read.csv(file.choose())
data[1:50,]
data$X<-data$Date1; data$Date<-NULL; colnames(data)[1]<-"Date"


#separate copepods, remove double-counts (KOD 15, 16, 128)
copepods<-data[which(data$Group_General=="Copepod"&data$Endemic_Cosmo=="endemic"),] #type in group of interest here
copepods<-copepods[which(copepods$KOD!=15|copepods$KOD!=16|copepods$KOD!=128),]


#create a volume column and total number of individuals column for later calculations
copepods$Volume<-as.numeric(levels(copepods$Nig_Gr)[copepods$Nig_Gr])-as.numeric(levels(copepods$Ver_Gr)[copepods$Ver_Gr])
copepods<-copepods[which(copepods$Volume!="NA"),] #removes records with no volume
copepods$Tot<-copepods$Volume*copepods$Count

#create single depth column by merging Ver_Gr and Nig_Gr
copepods$Depth<-paste(copepods$Ver_Gr,copepods$Nig_Gr, sep="-")
copepods$Ver_Gr<-NULL; copepods$Nig_Gr<-NULL

#create separate data sets for each depth category (0-10, 10-25, 50-100, 100-150, 150-250)
copepods0to10<-copepods[which(copepods$Depth=="0-10"),]; copepods0to10$Depth<-NULL
copepods10to25<-copepods[which(copepods$Depth=="10-25"),]; copepods10to25$Depth<-NULL
copepods50to100<-copepods[which(copepods$Depth=="50-100"),]; copepods50to100$Depth<-NULL
copepods100to150<-copepods[which(copepods$Depth=="100-150"),]; copepods100to150$Depth<-NULL
copepods150to250<-copepods[which(copepods$Depth=="150-250"),]; copepods150to250$Depth<-NULL

#alternatively, sum for all depths
copepodscommondepths<-copepods[which(copepods$Depth=="0-10"|copepods$Depth=="10-25"|copepods$Depth=="50-100"|copepods$Depth=="100-150"|copepods$Depth=="100-150"|copepods$Depth=="150-250"),]
library(reshape)
copmelt<-melt.data.frame(copepodscommondepths, id.vars=c("Date","KOD","Depth"),measure.vars=c("Tot","Volume"),na.rm=TRUE)
copepodsall<-cast(copmelt, formula=Date+KOD~variable,fun.aggregate=sum,margins=FALSE,add.missing=FALSE)
copepodsall$Count<-copepodsall$Tot/copepodsall$Volume

#for cluster analysis
library(reshape)
copepodscommondepths<-copepods[which(copepods$Depth=="0-10"|copepods$Depth=="10-25"|copepods$Depth=="50-100"|copepods$Depth=="100-150"|copepods$Depth=="100-150"|copepods$Depth=="150-250"),]
copmelt<-melt.data.frame(copepodscommondepths, id.vars=c("Date","KOD","Depth"),measure.vars=c("Tot","Volume"),na.rm=TRUE)
copepodsall<-cast(copmelt, formula=Date+KOD~variable,fun.aggregate=sum,margins=FALSE,add.missing=FALSE)
copepodsall$Count<-copepodsall$Tot/copepodsall$Volume
#allmelt<-melt.data.frame(copepodsall, id.vars=c("Date", "KOD"),measure.vars=c("Tot","Count","Volume"),na.rm=TRUE)
cope<-data.frame("Date"=copepodsall$Date, "Count"=copepodsall$Count, "KOD"=copepodsall$KOD)
clustcopes<-reshape(cope,idvar="Date",v.names="Count",timevar="KOD",direction="wide")
names(clustcopes)<-as.numeric(gsub("\\D", "", names(clustcopes)))  #take only numbers from column names


##################################################################
#########Creat a data frame summarizing by month###############
#####what depth range are you interested in (e.g. 0to10)?#########
##################################################################
library(zoo)
depthrange<-readline("enter depth range (0to10, 10to25, 50to100, 100to150, 150to250, all):           ")
interest<-paste("copepods",depthrange,sep="")

datofint<-get(interest)
datofint$Date<-paste(substr(datofint$Date,1,7),"-01",sep="") #convert all samples to 1st of month to make a regular time series
allcope<-ddply(datofint,~Date,summarise,sum=sum(Count)) #summarize count data by month (sum for all copepod species)
allcope$Date<-as.Date(allcope$Date)

#replace missing data with NA for processing (interpolation) in zoo package
years<-substr(allcope$Date,1,4)
start = as.Date(allcope$Date[1])
full <- seq(start, by="1 month", length=length(unique(years))*12)
allcopepods<-data.frame(Date=full, Count=with(allcope, sum[match(full, allcope$Date)]))
allcopepods<-allcopepods[which(allcopepods$Date>"1947-12-01"),] #take dates starting in 1948 to minimize missing values

#convert the data into a "zoo" time series object and replace NAs by interpolation
library(zoo)
nb<-zoo(allcopepods$Count,order.by=allcopepods$Date,frequency=12)
nb<-na.approx(nb, x = index(nb), na.rm = TRUE) #replaces NA by linear interpolation, na.rm=TRUE removes NA sequences on ends of series
plot(nb, xlab="Date", ylab="Abundance")
library(TSA)
periodogram(nb)

#remove linear trend from the time series
library(RSEIS)
newnb<-detrend(log(nb+1))
plot(newnb,col="red")


##########################
#follow with this code to# 
#convert data to QUARTERS#
##########################
library(xts)
locations<-endpoints(nb, "quarters")-1
locations<-locations[2:length(locations)]
copeq<-period.apply(nb, INDEX=locations, FUN=function(x) sum(x))
copepodquarters<-data.frame("Quarter"=seq(1:length(copeq)),"Count"=copeq)
plot(copepodquarters, type="l", xlab="Quarter",ylab="Abundance")


#####################################
#####process temperature data########
#####################################
library(plyr) #package for data manipulation

#read temp data in and get it in the right format 
temp<-read.csv(file="c:/Users/gray/Dropbox/Postdoc/R data/temperature19482002.csv"); names(temp)<-c("Date","Depth","Temp")
datofint<-temp ; datofint$Date<-as.Date(datofint$Date, origin="1899-12-30")
#datofint$Date<-paste(substr(datofint$Date,1,8),"01",sep="") #convert all samples to 1st of month to make a regular time series
datofint$Temp<-as.numeric(datofint$Temp)

#take only the surface temperatures (depth==0) #best sampled depths: 0, 25(OK), 50, 
library(plyr)
library(zoo)
surftemp<-datofint[which(datofint$Depth==0),]
alltemp<-ddply(surftemp,~Date,summarise,mean=mean(Temp)) #summarize by month by taking mean temp
alltemp$Date<-as.Date(alltemp$Date)


#replace missing data with NA for further processing 
years<-substr(alltemp$Date,1,4)
start = as.Date("1948-01-01")
full <- seq(start, by="1 month", length=length(unique(years))*12)
alltemperatures<-data.frame(Date=full, Count=alltemp[match(as.character(full), as.character(alltemp$Date)),])
alltemperatures<-data.frame("Date"=alltemperatures$Date, "Temp"=alltemperatures$Count.mean)
alltemperatures<-alltemperatures[which(alltemperatures$Date>"1955-01-01"),] #take dates after 1959 50

tp<-alltemperatures
missing<-which(is.na(tp$Temp)==TRUE)
for(i in 1:length(missing)){
  tp$Temp[missing[i]]<-mean(tp$Temp[which(substr(tp$Date,6,7)==substr(tp$Date,6,7)[missing[i]])],na.rm=TRUE)}
tp<-tp$Temp


#tp<-na.approx(tp, x = index(tp), na.rm = TRUE) #replaces NA by linear interpolation, na.rm=TRUE removes NA sequences on ends of series
plot(tp, xlab="Date", ylab="Temperature", type="l")
tp<-log10(tp+1);library(RSEIS); tp<-detrend(tp); nb<-tp
plot(tp~alltemperatures$Date, xlab="Date", ylab="log (temperature+1)", type="l")

tp0<-tp

tp50<-tp
stratindex<-(tp0/tp50)
nb<-stratindex

plot(tp0, type="l")
strantindex<-ts(stratindex,start=as.Date(alltemperatures$Date[1]),frequency=12)
plot(stratindex[12:24]~alltemperatures$Date[12:24], type="l", ylab="Stratification", xlab="Date")


#remove linear trend from the time series
library(RSEIS)
tp<-detrend(tp)
plot(tp,col="red")

#remove seasonality using traditional methods (if need be)
library(zoo)
tp<-ts(tp, start=index(tp)[1], frequency=12)
decomposedtp<-decompose(tp)
deseasonedtp<-ts(vv$x-vv$seasonal)


#######next plot missing months through time
x=(as.numeric(min(yearlevels))+1) : as.numeric(max(yearlevels))
y=unlist(lapply(missing, function(x) length(x)))
y<-y[2:length(y)]
plot(y~x, xlab="Sampling year", ylab="Number missing months")
barplot(y, names.arg=x, xlab="Year", ylab="Number missing months")

#most commonly observed species records
key<-read.csv(file.choose())
bob<-hist(as.numeric(data$KOD[which(data$Count!="0")]), breaks=length(unique(data$KOD)))
loc<-which(bob$counts>1500)
spp<-data.frame("Species"=key$English[loc], "Order"=key$Group_Fine[loc], "Status"=key$Endemic_Cosmo[loc])

#most commonly sampled depths
depths<-paste(data$Ver_Gr,data$Nig_Gr, sep="-")
barplot(table(depths), xlab="Depth (m)", ylab="Number of occurrences")





