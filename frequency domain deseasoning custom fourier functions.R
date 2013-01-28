
###################################################################################
####Frequency domain deseasoning with real discrete time fourier transform#########
###################################################################################

#testing code
x1<-2*cos(2*pi*1:100*6/100)+3*sin(2*pi*1:100*6/100)
x2<-4*cos(2*pi*1:100*10/100)+ 5*sin(2*pi*1:100*10/100)
x3<-6*cos(2*pi*1:100*40/100)+7*sin(2*pi*1:100*40/100)
x=x1+x2+x3
library(TSA)
uperiodogram(x)

#read data in, convert to time series, plot the data
data<-read.csv(file.choose()) #read in quarterly data from the GCB paper
series<-ts(data$Copepod,start=1948,frequency=4)
plot(log(series+1), type="l", xlab="Date",ylab="Abundance")

series<-x
plot(series, type="l")
#set window size, amount to slide window, and sampling frequency
windowsize = 10
increment=1 #how much to slide each window
freq=10 #frequency(series) #set annual harmonic frequency

#check if windows need to be zero padded
g<-windowsize:(windowsize+11)
reqpad<-as.numeric(which(is.wholenumber(g*(1/freq))==TRUE)) # How many zeroes to add

#break series into windows and add zeroes if necessary
windows<-list()
step=0
for (k in 1:(length(series)-(windowsize-1))){
  step=step+increment
  windows[[k]]<-series[step:(step+(windowsize-1))]
  windows[[k]]<-append(windows[[k]],rep(mean(series),(reqpad[1]-1)),after=length(windows[[k]])) #append zeroes at end of vector
}


deseasonedwindow<-list()
for (wind in 1:(length(windows))){
  
  bin= (windowsize*(1/freq))+1 #bin where the magnitude and phase of annual harmonic can be found
  
 #custom function for the discrete fourier transform
  realfft<-function(x){
    REX<-rep(0, (length(x)/2+1))
    IMX<-rep(0,(length(x)/2+1))
    
    for (k in 1:((length(x)/2)+1)){
      for (i in 1:length(x)){
        REX[k]<-REX[k]+x[i]*cos(2*pi*(k-1)*(i-1)/length(x))
        IMX[k]<-IMX[k]-x[i]*sin(2*pi*(k-1)*(i-1)/length(x))}}
    result<-complex(real=REX, imaginary=IMX)
    result
  }
  
  windowharmonicFD<-realfft(windows[[wind]])
  
  phase<-as.numeric(Arg(windowharmonicFD[bin])) #phase of annual harmonic
  magnitudes<-as.numeric(Mod(windowharmonicFD))[2:length(windowharmonicFD)] #magnitudes of all frequency bins except 0
  avgmagnitude<-mean(magnitudes[which(magnitudes!=magnitudes[bin-1])]) #mean of all magnitudes except for the bin with the annual harmonic
  
  ##Alternative ways of calculating the mean noise level
  #avgmagnitude<-mean(magnitudes[which(magnitudes!=magnitudes[bin])]) #not sure whether to include the first magnitude?
  #avgmagnitude<-mean(magnitudes[(bin+1):length(magnitudes)]) #not sure whether to include the first magnitude?
  
  #replace values for the annual harmonic with the mean value "avgmagnitude"
  windowharmonicFD[bin]<-complex(real=avgmagnitude*cos(phase), imaginary=avgmagnitude*sin(phase))
  
  
  #custom function for inverse discrete fourier transform
  irealfft<-function(y){
    REX<-c()
    IMX<-c()
    
    REX<-Re(y)
    IMX<-Im(y)
    
    REX1<-c()
    IMX1<-c()
    
    for (i in 1:(length(REX))){
      REX1[i]<-REX[i]/(length(REX)-1)
      IMX1[i]<-(-1*IMX[i]/(length(REX)-1))}
    REX1[1]<-REX[1]/((length(REX)-1)*2)
    REX1[length(REX1)]<-REX[length(REX)]/((length(REX)-1)*2)
    
    result<-c()  
    
    for (i in 1:((length(REX)-1)*2)){
      xx=0
      for (k in 1:(length(REX))){
        xx=xx+(REX1[k]*cos((2*pi*(k-1)*(i-1))/((length(REX)-1)*2)))+ (IMX1[k]*sin((2*pi*(k-1)*(i-1))/((length(REX)-1)*2)))}
      result[i]<-xx
    }
    result
  }
    
  deseasonedwindow[[wind]]<-irealfft(windowharmonicFD) #inverse fourier transform to return to the time domain
}

#overlap save to reconstruct the whole series
reconstructed<-deseasonedwindow[[1]][1:(windowsize/2)]
reconstructed<-append(reconstructed,deseasonedwindow[[2]][(windowsize/2):windowsize],after=length(reconstructed))
for (j in 3:length(deseasonedwindow)){
  lng<-length(reconstructed)
  reconstructed<-reconstructed[1:(lng-(windowsize/2))]
  reconstructed<-append(reconstructed,deseasonedwindow[[j]][(windowsize/2):windowsize],after=length(reconstructed))
}

#make the reconstructed data into a time series and plot the deseasoned data (points) and original data (lines)
reconstructed<-ts(reconstructed,start=1948,frequency=4)
plot(log(reconstructed+1), type="p", pch=16, ylab="ln(deaseasoned copepod abundance +1)", ylim=c(0,5), xlim=c(1945,2005))
lines(log(series+1),col="grey")

#Add rectangles for the PDO to the plot
rect(1957,-10,1961,8, density=5, col="green")
rect(1976,-10,1988,8, density=5, col="green")
rect(1992,-10,1998,8, density=5, col="green")
