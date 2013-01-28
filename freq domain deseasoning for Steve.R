#########################################
####Frequency domain deseasoning#########
#########################################

newnb<-ts(data$Copepod,start=1948,frequency=4)

plot(log(newnb+1), type="l", xlab="Date",ylab="Abundance")

freq=4 #set sampling frequency

#divide data into windows of length "windowsize"
windows<-list()
windowsize = 40
winlength=windowsize-1 #this number should be desired window length -1
increment=1
step=0

#check if windows need to be padded
g<-windowsize:(windowsize+11)
reqpad<-as.numeric(which(is.wholenumber(g*(1/freq))==TRUE)) # How many zeroes to add

for (k in 1:(length(newnb)-winlength)){
  step=step+increment
  windows[[k]]<-newnb[step:(step+winlength)]
  windows[[k]]<-append(windows[[k]],rep(0,(reqpad[1]-1)),after=length(windows[[k]])) #append zeroes at end of vector
}

#loop though all the windows and reduce the annual harmonic to the mean noise level
deseasonedwindow<-list()
for (wind in 1:(length(windows))){
  
  bin= (windowsize*(1/freq))+1
  windowharmonicFD<-fft(windows[[wind]]) 
  
  
  magnitudewindow<-as.numeric(Mod(windowharmonicFD[bin])) #magnitude of annual harmonic
  phase<-as.numeric(Arg(windowharmonicFD[bin])) #phase of annual harmonic
  magnitudes<-as.numeric(Mod(windowharmonicFD[1:(length(windowharmonicFD)/2)]))
  avgmagnitude<-mean(magnitudes[which(magnitudes!=magnitudes[bin])]) #not sure whether to include the first magnitude?
    
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
  deseasonedwindow[[wind]]<-Re((fft(aas, inverse=TRUE)))/length(aas) #inverse fft to return to time domain
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

#plot the result
plot(log(reconstructed+1), type="p", pch=16, ylab="ln(deaseasoned copepod abundance +1)", ylim=c(0,5), xlim=c(1945,2005))
lines(log(newnb+1),col="grey")

#Rectangles for the PDO
rect(1957,-10,1961,8, density=5, col="green")
rect(1976,-10,1988,8, density=5, col="green")
rect(1992,-10,1998,8, density=5, col="green")
