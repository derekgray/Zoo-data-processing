require(doBy)

#Phytoplankton group of interest
phytgrp<-"diatom"

#should NA dates be filled?
fillNAs<-"no" 

#Start date
stdate<-as.Date("1974-01-01") #most consistent for depths and after preservation change

#End date
enddate<-as.Date("1999-12-31")

if(exists("phyto")==FALSE){phyto<-read.csv(file="phytoplankton.csv", stringsAsFactors=F)
                           phyto$DATE<-as.Date(phyto$DATE,"%m/%d/%Y")
                           names(phyto)<-c("Date", "Month","Year","Depth","Code","Count","Group","Genus","Species")
                           phyto<-orderBy(Date~.,data=phyto)
                           phyto$Group[which(phyto$Group=="chrysophtye")]<-"chrysophyte"
                           phyto$Group[which(phyto$Group=="Bacteria")]<-"cyano"
                           phyto$Group[which(phyto$Group=="Flagellates")]<-"Flagellate"
                           phyto$Group[which(phyto$PhytoGenus=="Romeria")]<-"PicoAlgae unID"
                           phyto$Group[which(phyto$PhytoGenus=="Synechococcus")]<-"PicoAlgae unID"
                           
                           #Read in phyto key
                           key<-read.csv("Baikal_phyto_Sp_Key_20120725.csv", stringsAsFactors=F)
                           
                           #replace species names with those from the key
                           for (i in 1:nrow(key)){
                             phyto$Species[which(phyto$Code==key$Code[i])]<-paste(key$Genus[i],key$Species[i])
                           }
                           
}



#take select depths that have been sampled regularly throughout the program
phyto2<-phyto[which(phyto$Depth==0|phyto$Depth==10|phyto$Depth==50|phyto$Depth==100|phyto$Depth==200),] #5, 250? |phyto$Depth==150

#Take data after the specified start date
phyto2<-phyto2[which(phyto2$Date>=stdate&phyto2$Date<=enddate),]

#Some functions needed for processing----------------------

#identify NAs and then fill with monthly average for whole time series
NAfill<-function(x){
  years<-substr(x$Date,1,4)
  start = as.Date(x$Date[which.min(x$Date)])
  full <- seq(start, by="1 month", length=((max(as.numeric(years))+1)-min(as.numeric(years)))*12)
  allsp<-data.frame("Date"=full)
  for (i in 1:ncol(x)){
    t1<-data.frame("Count"=with(x, x[,i][match(full, as.Date(x$Date))]))
    allsp<-cbind(allsp, t1)}
  allsp<-allsp[,2:ncol(allsp)]; names(allsp)<-names(x); allsp$Date<-full
  #names(allsp)<-c("Date","strat")
  missing<-allsp$Date[which(is.na(allsp[,2])==T)] 
  for (i in 1:length(missing)){
    allsp[,2][which(allsp$Date==as.Date(missing[i]))]<-mean(allsp[,2][which(as.numeric(substr(allsp$Date,6,7))==as.numeric(substr(missing[i],6,7)))],na.rm=T)
  }
  return(allsp)
}

NAdates<-function(x){
  years<-substr(x$Date,1,4)
  start = as.Date(x$Date[which.min(x$Date)])
  full <- seq(start, by="1 month", length=((max(as.numeric(years))+1)-min(as.numeric(years)))*12)
  allsp<-data.frame("Date"=full)
  for (i in 1:ncol(x)){
    t1<-data.frame("Count"=with(x, x[,i][match(full, as.Date(x$Date))]))
    allsp<-cbind(allsp, t1)}
  allsp<-allsp[,2:ncol(allsp)]; names(allsp)<-names(x); allsp$Date<-full
  return(allsp)
}

#---------------------------------------------------------

#Check which phyto group was requested
if (phytgrp!="all"){phyto2<-phyto2[which(phyto2$Group==phytgrp),]} #which phyto group was requested?
if (phytgrp=="all"){phyto2<-phyto2[which(phyto2$Group=="diatom"|phyto2$Group=="cyano"|phyto2$Group=="green"|phyto2$Group=="chrysophtye"|phyto2$Group=="cryptophyte"|phyto2$Group=="dinoflagellate"|phyto2$Group=="PicoAlgae unID"),]} #which phyto group was requested?

#summarize data by density according to date and depth
phyto2$Date<-paste(substr(phyto2$Date,1,7),"-01",sep="") #make all dates 1st of month for easy monthly averages using summaryBy

#check if all depths were sampled on a particular date and exclude dates that have missing depths
exc<-c()
da<-c()
cat<-list()
undates<-unique(phyto2$Date)
for (i in 1:length(undates)){
  this<-phyto2[which(phyto2$Date==undates[i]),]
  if(nrow(this)>0){
    cat[[i]]<-(unique(this$Depth))
    da[i]<-length(unique(this$Depth))
    if(da[i]==length(unique(phyto2$Depth))){exc[i]<-"N"}; if(da[i]<length(unique(phyto2$Depth))){exc[i]="Y"}}
}
gooddates<-undates[which(exc=="N")]
phyto2<-phyto2[which(phyto2$Date%in%gooddates),]

#Figure out densities throught ime
mnfun<-function(x){mean(x,na.rm=T)}
phytdens<-summaryBy(Count~Date+Species,data=phyto2, FUN=mnfun)
phytdens$Date<-as.Date(phytdens$Date)
names(phytdens)<-c("Date","Species","Density")

#Reorganize phytdens dataframe so that dates are in the first row and species across the top
library(reshape)
temp<-melt(phytdens, id=c("Species","Date"), measure.vars="Density")
clustall<-cast(temp, Date~Species)
row.names(clustall)<-clustall$Date; clustall$Date<-NULL

#Divide density data into quarters using the xts package
library(xts)
for (i in 1:ncol(clustall)){
  qtnb<-xts(x=clustall[,i], order.by=as.Date(row.names(clustall)))
  locations<-endpoints(qtnb, "quarters")
  if(i==1){
    qtrnb<-as.data.frame(period.apply(qtnb, INDEX=locations, FUN=function(x) mean(x,na.rm=T)))
  }
  if(i>1){
    qtrnb[,i]<-period.apply(qtnb, INDEX=locations, FUN=function(x) mean(x,na.rm=T))
  }
}
names(qtrnb)<-names(clustall)
  
#Make numbers for the quarters (1-4)
quarters<-(as.numeric(substr(row.names(qtrnb),6,7)))
quarters[which(quarters==1|quarters==2|quarters==3)]<-1
quarters[which(quarters==4|quarters==5|quarters==6)]<-2
quarters[which(quarters==7|quarters==8|quarters==9)]<-3
quarters[which(quarters==10|quarters==11|quarters==12)]<-4
qtrnb$quarter<-quarters

#insert a column that displays the quarters
row.names(qtrnb)<-paste(substr(row.names(qtrnb),1,4),qtrnb$quarter)

#separate the dataset into four dataframes: winter, spring, summer, and fall
winter<-as.data.frame(qtrnb[which(qtrnb$quarter==1),]); winter$year<-substr(row.names(winter),1,4)
spring<-qtrnb[which(qtrnb$quarter==2),];spring$year<-substr(row.names(spring),1,4)
summer<-qtrnb[which(qtrnb$quarter==3),];summer$year<-substr(row.names(summer),1,4)
fall<-coredata(qtrnb[which(qtrnb$quarter==4),]);fall$year<-substr(row.names(fall),1,4)

library(vegan)
prop.peryear<-decostand(summer[,1:(ncol(summer)-2)], method="total")
byspecies<-data.frame("Species"=names(summer)[1:(ncol(summer)-2)],"Proportion"=colMeans(prop.peryear))
byspecies<-byspecies[which(byspecies$Proportion>0),]

byspecies<-orderBy(Proportion~., data=byspecies)

par(oma = c(3,0,0,0))
par(mar=c(6,8,3,3))
mp<-barplot(byspecies$Proportion[which(byspecies$Proportion>0.02)], ylab="Proportion", space=0, main="Summer diatoms")

text(mp, par("usr")[3] - 0.001, srt = 20, adj = 1, 
     labels = byspecies$Species[which(byspecies$Proportion>0.02)], 
     xpd = TRUE, font = 0.25) 
abline(h=0.01)

abundant.diatoms<-as.character(byspecies$Species[which(byspecies$Proportion>0.02)])
colr<-c("black","red","blue","orange","green","yellow","purple","cyan","magenta","pink")

par(mfrow=(c(3,3)))
for (i in 1:length(abundant.diatoms)){
plot(summer[,which(names(summer)==abundant.diatoms[i])], type="l", ylab="Density", main=abundant.diatoms[i], xlab="Year")
}
