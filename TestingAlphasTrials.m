%Needed variables

%Xorig3 = Another clean data set of size channels X times X trials
%A = a mixing matrix from which you have selected one artifact topography 


%%
%Setting values for parameter alpha. Eq (6)
alphas=[0:0.1:1]; 

clear S Fmeasure FmeasureAve measVariability1 measVariability2 ErrorMeasureAve topoCC 

%different numbers of trials
Rnumbers=[50 150 300 500 750 1000];

%loops over parapemeter and repetitions
for jR=1:length(Rnumbers)
for j=1:length(alphas)
    
for n=1:100 %complete simulation rounds 

% Simulation begins

alpha=alphas(j);

R=Rnumbers(jR);
N=1001%time points 
t=1:N;%time axis 

S=zeros(1,N,R); %the artifact component over trials and time
t0=400; %artefact peak time
ht=50; %artifact 'width'

for k=1:R
    
phi0=pi/2-t0*pi/100; %phase at time zero

%artifact waveform in trial k
Stemp(t)=(1-alpha)*sin(pi/100*t+phi0)+alpha*sin(pi/100*t+rand(1)*2*pi);
%multiply by Gaussian envelope
Stemp=Stemp.*exp(-(t-t0).^2/ht.^2);

Stemp=Stemp./std(Stemp,1); %normalize to std=1
S(1,:,k)=Stemp; 

end

amplitudeArt=200;% scaling amplitude of artifact
Aart=A(:,13); %preselected artifact topography from another mixing matrix A
Aart=Aart*diag(1./sqrt(sum(Aart.^2)))*amplitudeArt; %normalize magnitude and scale

%EEG:n simulointi
indsR=randperm(size(Xorig3,3), R); %randomize trial indices

Xorig=reshape(Xorig3(:,:,indsR),size(Xorig3,1),[]); %concatenate trials

X=Aart*reshape(S,1,[])+Xorig; %create artifact data and superimpose onto clean data

%run ICA. 
[A2, W, icasig, Proj] = fastica_nodemean_autoCompression(X, 200, 45); %max 200 iterations, find 45 components
    
icaS=icasig./std(icasig,1,2); %normalize
Sorig=reshape(S,1,[]); %concatenate simulated artifact waveforms trials

crossC=icaS*Sorig'/(N*R); %cross-correlation vector betweeen the simulated and all the estimated components

%sort by magnitute of correlation
[~, isort]=sort(abs(crossC), 'descend') ;

% best-matching estimated component represent the artifact, and will be
% removed to get clean data X1
X1=X-A2(:,isort(1))*icasig(isort(1),:);

%define the window where the artifact exists.
tWindow=t0+(-ht:ht);
tWindow=tWindow(tWindow>0 & tWindow <=N);


%Original clean data
XorigAve=mean(reshape(Xorig, size(A,1), N,R),3); 
% estimated clean data.
X1ave=mean(reshape(X1, size(A,1), N,R),3);
%Difference
XdeltaAve=XorigAve-X1ave;
%Relative error. Frobeniusnorm of the difference relative to that of the original clean data
FmeasureAve(j,jR,n)= sqrt(sum(sum((Proj*XdeltaAve(:,tWindow)).^2)))...
    /sqrt(sum(sum((Proj*XorigAve(:,tWindow)).^2)))*100;

%variability measures
swave=reshape(icasig(isort(1),:), N, R); % estimated artifact waveform
ssave=mean(swave,2); % 
measVariability2(j,jR,n)=sqrt(mean(mean((swave-repmat(ssave, [1 R])).^2))); %variability measure 1

swave=swave-mean(swave,1);% demeaning 
measVariability1(j,jR,n)=sqrt(mean((ssave-mean(ssave)).^2))/sqrt(mean(swave(:).^2)); %variability measure 2


end
end
end
%
