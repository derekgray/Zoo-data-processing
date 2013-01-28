


Nproc = 512;
xproc = rep(0, 2*Nproc) # %# initialize temp buffer
idx = 1:Nproc #; %# initialize half-buffer index
ycorrect = zeros(2*Nproc,1); %# initialize destination

for ctr = 1:(length(x)/Nproc) %# iterate over x 512 points at a time
    xproc(1:Nproc) = xproc((Nproc+1):end); %# shift 2nd half of last iteration to 1st half of this iteration
    xproc((Nproc+1):end) = x(idx); %# fill 2nd half of this iteration with new data
    yproc = ifft(fft(xproc) .* fft(h,2*Nproc)); %# calculate new buffer
    ycorrect(idx) = real(yproc((Nproc+1):end)); %# keep 2nd half of new buffer
    idx = idx + Nproc; %# step half-buffer index
end