function [PSDRAW,freqRAW]=PSDofRAW(DATA,POINTS,delt,UPPERfLIM)
%PSDRAW=Output PSD
%freqRAW=Corresponding PSD frequencies
%DATA=Input raw wind data array
%POINTS=Number of frequencies to calculare the PSD at between 0 Hz and
%the
%Nyquist Frequency. For speed, will be trimmed to the nearest value
%that
%satisfies POINTS=2^p where p=integer. Generally, 128, 256, and 512
%are
%adequate values of POITNS
%delt=Time between samples of DATA
% UPPERfLIM=Highest limit at which to calculate the PSD, if UPPERfLIM=0
%Hz,
% all frequencies will be calculated from 0 Hz to the Nyquist Frequency
%% Orient Raw Data Arrays
if size(DATA,2) > 1
DATA=DATA';
end
%% Break Data Into Segments
N=2^nextpow2(POINTS); %Transfer POINTS to POINTS=n^p where p=Integer
nd=floor(length(DATA)/N);%Determine number of segments within DATA
%array
WINDDATA=DATA(1:nd*N);%Remove excess data points that cannot make up a
%full segment
%Win= length(WINDDATA)
WINDOW=(1-cos(pi*(0:N-1)/N).^2)'; %Create Hanning Window
SEGMENTS=zeros(nd,POINTS);%Initialize Segment vector
%% Separate Segments into a 2D matrix of (nd x N) dimensions while also
%applying the window to the segments
for j1=1:nd
CHUNK=(WINDDATA(1+N*(j1-1):j1*N));
SEGMENTS(j1,:)=CHUNK.*WINDOW;
end 
%% Create frequency vector
k=0:N/2;
n=0:N-1;
freqRAW=k/(N*delt); %frequency values to calculate the PSD
FFT=zeros(1,length(n));
power=zeros(1,length(k));
POWTOT=zeros(nd,length(k)); %Initialize Vectors
PSDRAW=zeros(1,length(k));
%% Calculate the PSD of each Segment
for f=1:nd
for kk=1:length(k)
for nn=1:length(n)
FFT1(nn)=SEGMENTS(f,nn)*exp(-i*2*pi*k(kk)*n(nn)/N);
end
FFT(kk)=sqrt(8/3)*sum(FFT1)*delt;
power(kk)=abs(FFT(kk))^2;
end
POWTOT(f,:)=power;
end
%% Average all Segments Together
for h=1:length(k)
PSDRAW(h)=2*sum(POWTOT(:,h))/(nd*N*delt);
end
%% Reorient Arrays
freqRAW=freqRAW';
PSDRAW=PSDRAW';
if isnan(UPPERfLIM) == 0
PSDRAW(freqRAW>UPPERfLIM)=[];
freqRAW(freqRAW>UPPERfLIM)=[];
end