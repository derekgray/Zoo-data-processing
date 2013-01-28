######################################################################
##############Alternative way to make data frame######################
######################################################################

overall=list()
missing=list()

years=substr(get(interest)$Date,1,4) #strip month and day and leave only year
years<-factor(years) #determine how many years the data spans
yearlevels=levels(years)


for (j in 2:length(yearlevels)){ #loop over the number of years
previousyear=get(interest)[which(years==yearlevels[j-1]),]
currentyear=get(interest)[which(years==yearlevels[j]),]

#create a vector of months for the current year
month=substr(currentyear$Date,6,7)#isolates months
month=as.numeric(month)

#determine which months are missing from the data set
miss<-setdiff(seq(1:12),unique(month)); try(miss<-miss[which(miss!=12)])
if (length(which(substr(previousyear$Date,6,7)=="12"))==0){miss<-append(miss,12,after=length(miss))}


#make dataframe for winter
wintera<-previousyear[which(substr(previousyear$Date,6,7)=="12"),]
winterb<-currentyear[which(substr(currentyear$Date,6,7)=="01"),]
winterc<-currentyear[which(substr(currentyear$Date,6,7)=="02"),]
winter<-rbind(wintera,winterb,winterc)
library(doBy)
if (nrow(winter)>0){winter<-summaryBy(Count~KOD, data = winter, id = NULL, FUN = sum, order=FALSE)}

#make dataframe for spring, summer, fall
spring<-currentyear[which(month==3|month==4|month==5),]
if (nrow(spring)>0){spring<-summaryBy(Count~KOD, data = spring, id = NULL, FUN = sum, order=FALSE)}

summer<-currentyear[which(month==6|month==7|month==8),]
if (nrow(summer)>0){summer<-summaryBy(Count~KOD, data = summer, id = NULL, FUN = sum, order=FALSE)}

fall<-currentyear[which(month==9|month==10|month==11),]
if (nrow(fall)>0){fall<-summaryBy(Count~KOD, data = fall, id = NULL, FUN = sum, order=FALSE)}


#combine quarters into single dataframe

#all possibilities for missing first and second quarters
if (nrow(winter)>0 & nrow(spring)>0) {alls1<-merge(winter, spring, by="KOD",all.y=TRUE, all.x=TRUE) ; names(alls1)<-c("KOD",paste(yearlevels[j],"WINTER"),paste(yearlevels[j],"SPRING"))}
if (nrow(winter)>0 & nrow(spring)==0) {alls1<-winter; names(alls1)<-c("KOD",paste(yearlevels[j],"WINTER"))}; if (nrow(winter)==0 & nrow(spring)>0) {alls1<-spring; names(alls1)<-c("KOD",paste(yearlevels[j],"SPRING"))}
if (nrow(winter)==0 & nrow(spring)==0) {alls1<-c()}

#all possibilities for missing third and fourth quarters
if (nrow(summer)>0 & nrow(fall)>0) {alls2<-merge(summer, fall, by="KOD",all.y=TRUE, all.x=TRUE); names(alls2)<-c("KOD",paste(yearlevels[j],"SUMMER"),paste(yearlevels[j],"FALL"))}
if (nrow(summer)>0 & nrow(fall)==0) {alls2<-summer; names(alls2)<-c("KOD",paste(yearlevels[j],"SUMMER"))}; if (nrow(summer)==0 & nrow(fall)>0) {alls2<-fall; names(alls2)<-c("KOD",paste(yearlevels[j],"FALL"))}
if (nrow(summer)==0 & nrow(fall)==0) {alls2<-c()}


if (length(alls1)>0 & length(alls2)>0) {allseasons<-merge(alls1, alls2, by="KOD",all.y=TRUE, all.x=TRUE)}
if (length(alls1)>0 & length(alls2)==0) {allseasons<-alls1} ; if (length(alls1)==0 & length(alls2)>0) {allseasons<-alls2}

missing[[j]]<-miss  #list of missing months sampled
overall[[j]]<-allseasons #list of species IDs and counts by quarter

miss<-c(); WINTER<-c(); SPRING<-c(); SUMMER<-c(); FALL<-c(); miss<-c(); allseasons<-c(); winter<-c(); spring<-c(); summer<-c(); fall<-c()
} 


#merge lists for all years
h1<-data.frame("species"=overall[[2]]$KOD, overall[[2]][2])

for (g in 2:(length(overall)-1)){
h2<-data.frame("species"=overall[[(g+1)]]$KOD, overall[[(g+1)]][2:length(overall[[(g+1)]])])
h1<-merge(h1, h2, by="species",all.y=TRUE, all.x=TRUE) }
data<-h1
data[is.na(data)] <- 0
rownames(data)<-data$species
data$species<-NULL
data<-data[order(row.names(data)),]
names(data)<-sub("X", "", names(data), ignore.case = FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE)

########Make a key to identify unique quarters############
yearss<-rep(yearlevels[1],4)
seasons<-c("WINTER","SPRING","SUMMER","FALL")
for (i in 2:length(yearlevels)){
te<-rep(yearlevels[i],4)
yearss<-append(yearss,te,after=length(yearss))}
quarterkey<-data.frame("season"=rep(seasons,(length(yearlevels))),"quarter"=seq(1:(length(yearlevels)*4)), "Year"=yearss)
quarterkey$combo<-paste(quarterkey$Year,quarterkey$season,sep=".")
quarterkey$Year<-NULL; quarterkey$season<-NULL

#########replace column headings with quarter number#######
for (i in 1:length(names(data))){
names(data)[i]<-as.numeric(quarterkey$quarter[which(quarterkey$combo==names(data)[i])])}


#########linearly interpolate missing values################
intdata<-approx(as.numeric(names(data)), y = colSums(data), xout=(1:max(quarterkey$quarter)),method="linear", ties = ordered)


##########plot interpolated data versus data with missing points##########
par(mfrow = c(2,1))
plot(approx(as.numeric(names(data)), y = colSums(data), xout=(1:max(quarterkey$quarter)),method="linear", ties = ordered),type="l")


#plot total copepod abundance through time
plot(colSums(data)~ as.numeric(names(data)), xlab="Year", ylab="Number individuals per litre", type="l")
