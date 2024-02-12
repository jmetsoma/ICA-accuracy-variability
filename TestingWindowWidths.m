% .........................................................................
% 12.2.2024 Johanna Metsomaa, NBE, Aalto university
% email johanna.metsomaa@gmail.com
% .........................................................................


%Needed variables

%Xorig3 = Another clean data set of size channels X times X trials
%A = a mixing matrix from which you have selected one artifact topography 
%%
%time window width where artifact occurs
t0_intervals=[0 2 5 8 12 16 25 50 100 200 500 1001]; 

clear S Fmeasure FmeasureAve measVariability1 measVariability2 
Xorig=reshape(Xorig3,size(Xorig3,1),[]); %concatenate trials of the data into 2D

alpha=0; %alpha parameter is 0 (deterministic waveform appearance).

R=1000;%number of trials 
N=1001% Number of time points
t=1:N;%time axis

S=zeros(1,N,R); %memory for artifact component
tm0=400; %average latency (midpoint of the artifact time window)
ht=50; %artifact 'width'
randLatency=true; %random latency within the given window

%for loops over parameters and iterations
for n = 1:100 %complete simulation rounds 
for j=1:length(t0_intervals) %widow widths
    t0_interval=t0_intervals(j);
for k=1:R %trials
    if randLatency 
        % randomize peak latency within the time window 
        t1=max(1, round(tm0-(t0_interval/2))); 
        deltat=1-round(tm0-t0_interval/2);
        t2=min(N, round((tm0+t0_interval/2)+deltat*(deltat>0) ));
        t0=t1+rand(1)*(t2-t1);
        
    end

phi0=pi/2-t0*pi/100; %phase at 0 ms

%simulated artifact waveform
Stemp(t)=(1-alpha)*sin(pi/100*t+phi0)+alpha*sin(pi/100*t+rand(1)*2*pi);
%multiply by envelope (Gaussian)
Stemp=Stemp.*exp(-(t-t0).^2/ht.^2);

%normalize to std = 1
Stemp=Stemp./std(Stemp,1); 
S(1,:,k)=Stemp; %save the component in 3D matrix
end

amplitudeArt=200;% artifact scaling factor
Aart=A(:,13); % predefined artifact topography
Aart=Aart*diag(1./sqrt(sum(Aart.^2)))*amplitudeArt; %normalize the topography and scale

%simulate artifactual EEG by to the clean EEG (Xorig), the simulated
%artifact
X=Aart*reshape(S,1,[])+Xorig; 

%Run ICA. 
[A2, W, icasig, Proj] = fastica_nodemean_autoCompression(X, 200, 45); %max 200 iterations, 45 components
    
icaS=icasig./std(icasig,1,2); %normalize
Sorig=reshape(S,1,[]); %concatenate simulated artifact waveforms trials

crossC=icaS*Sorig'/(N*R); %cross-correlation vector betweeen the simulated and all the estimated components

%sort by magnitute of correlation
[~, isort]=sort(abs(crossC), 'descend') ;

% best-matching estimated component represent the artifact, and will be
% removed to get clean data X1
X1=X-A2(:,isort(1))*icasig(isort(1),:);


%time window of artifact
tWindow= (t1-ht):(t2+ht); 
tWindow=tWindow(tWindow>0 & tWindow <=N);


%Original clean data
XorigAve=mean(reshape(Xorig, size(A,1), N,R),3); 
% estimated clean data.
X1ave=mean(reshape(X1, size(A,1), N,R),3);
%Difference
XdeltaAve=XorigAve-X1ave;
%Relative error. Frobeniusnorm of the difference relative to that of the original clean data
FmeasureAve(j,n)= sqrt(sum(sum((Proj*XdeltaAve(:,tWindow)).^2)))...
    /sqrt(sum(sum((Proj*XorigAve(:,tWindow)).^2)))*100;

%variability measures
swave=reshape(icasig(isort(1),:), N, R); % estimated artifact waveform
ssave=mean(swave,2); 
measVariability2(j,n)=sqrt(mean(mean((swave-repmat(ssave, [1 R])).^2))); %variability measure 1

swave=swave-mean(swave,1);% demeaning 
measVariability1(j,n)=sqrt(mean((ssave-mean(ssave)).^2))/sqrt(mean(swave(:).^2));%variability measure 2
end
end

%%
