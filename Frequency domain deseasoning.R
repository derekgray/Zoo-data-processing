#get zooplankton data

# Read data into R
data<-read.csv(file.choose())
data[1:50,]
head(data)

#remove seasonality using traditional methods (if need be)
library(zoo)
series<-ts(data$SurfTemp,start=1948,frequency=4)
tp<-zoo(series,order.by=index(series),frequency=4)
tp<-na.approx(tp, x = index(tp), na.rm = TRUE) #replaces NA by linear interpolation, na.rm=TRUE removes NA sequences on ends of series
plot(tp, xlab="Date", ylab="Temperature")
vv<-decompose(as.ts(tp))
deseasonedtp<-ts(vv$x-vv$seasonal)
plot(deseasonedtp)
deseasonedtp<-ts(deseasonedtp, start=1948, frequency=4)

#smooth the time series(if desired)
library(TSA)
bob<-lowess(deseasonedtp~index(deseasonedtp), f=.03)
plot(bob,type="l", xlab="Time",ylab="Abundance", ylim=c(0,7))



#########################################
####Frequency domain deseasoning#########
#########################################

series<-ts(data$Copepod,start=1948,frequency=4)

series<-ts(data$SurfTemp,start=1948,frequency=4)
          
plot(log(series+1), type="l", xlab="Date",ylab="Abundance")

freq= 1/4 #annual harmonic

#divide data into windows of length "windowsize"
windows<-list()
windowsize = 40
winlength=windowsize-1 #this number should be desired window length -1
increment=1
step=0

#check if windows need to be padded
g<-windowsize:(windowsize+11)
reqpad<-as.numeric(which(is.wholenumber(g*(freq))==TRUE)) # How many zeroes to add

for (k in 1:(length(series)-winlength)){
  step=step+increment
  windows[[k]]<-series[step:(step+winlength)]
  windows[[k]]<-append(windows[[k]],rep(mean(series),(reqpad[1]-1)),after=length(windows[[k]])) #append zeroes at end of vector
}


deseasonedwindow<-list()
for (wind in 1:(length(windows))){
  
  bin= (windowsize*(freq))+1
  windowharmonicFD<-fft(windows[[wind]]) 
    
  magnitudewindow<-as.numeric(Mod(windowharmonicFD[bin])) #magnitude of annual harmonic
  phase<-as.numeric(Arg(windowharmonicFD[bin])) #phase of annual harmonic
  magnitudes<-as.numeric(Mod(windowharmonicFD[1:(length(windowharmonicFD)/2)]))
  avgmagnitude<-mean(magnitudes[which(magnitudes!=magnitudes[bin])]) #not sure whether to include the first magnitude?
  #avgmagnitude<-mean(magnitudes[(bin+1):length(magnitudes)]) #not sure whether to include the first magnitude?
  avgmagnitude<-0
  #magnitudes<-magnitudes[2:length(magnitudes)]
  #avgmagnitude<-mean(magnitudes[which(magnitudes!=magnitudes[bin])]) #not sure whether to include the first magnitude?
  
  
  realwindowharmonicFD<-fft(windows[[wind]])[1:((length(windows[[wind]])/2)+1)]
  realwindowharmonicFD[bin]<-complex(real=avgmagnitude*cos(phase), imaginary=avgmagnitude*sin(phase))
  
  #start with the real fft spectrum, convert to complex fft representation and then take the inverse
 
  REX<-c()
  IMX<-c()
  
  REX<-as.numeric(Re(realwindowharmonicFD))
  IMX<-as.numeric(Im(realwindowharmonicFD))
  
  for (k in (length(windows[[wind]])/2+2):(length(windows[[wind]]))){
    REX[k] =  REX[(length(windows[[wind]])+2)-k]
    IMX[k] = -IMX[(length(windows[[wind]])+2)-k]
  }
  
  aas<-complex(real=REX, imaginary=IMX)
  deseasonedwindow[[wind]]<-Re((fft(aas, inverse=TRUE)))/length(aas)
}

#overlap save to reconstruct the whole series
reconstructed<-deseasonedwindow[[1]][1:(windowsize/2)]
reconstructed<-append(reconstructed,deseasonedwindow[[2]][(windowsize/2):windowsize],after=length(reconstructed))
for (j in 3:length(deseasonedwindow)){
  lng<-length(reconstructed)
  reconstructed<-reconstructed[1:(lng-(windowsize/2))]
  reconstructed<-append(reconstructed,deseasonedwindow[[j]][(windowsize/2):windowsize],after=length(reconstructed))
}
reconstructed<-ts(reconstructed,start=1948,frequency=4)
plot(log(reconstructed+1), type="p", pch=16, ylab="ln(deaseasoned copepod abundance +1)", ylim=c(0,5), xlim=c(1945,2005))
lines(log(series+1),col="grey")

#Rectangles for the PDO
rect(1957,-70,1961,80, density=5, col="green")
rect(1976,-70,1988,80, density=5, col="green")
rect(1992,-70,1998,80, density=5, col="green")
