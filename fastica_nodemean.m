function [A, W, icasig, P, Mwhite, Wwhite] = fastica_nodemean(mixedsig, maxNumIterations)
%
% FASTICA - Fast Independent Component Analysis without centering/demeaning
% step
% note fastica toolbox is needed to run this function
%
% input:
% mixedsig: all data (channels x times x epochs)
% maxNumIterations: max number of iterations for optimizing fastica
%
% output: (all component variable are sorted from highest to smallest
% variance)
% A: mixing matrix
% W: demixing
% icasig: component waveforms
% P: projection matrix (if first pca component was > 1, this is important
% Mwhite: whitening matrix
% White: demixing matrix for whitened data
% .........................................................................
% 29 March 2021 : Johanna Metsomaa, BNP, University of Tübingen  
% .........................................................................
%

[C,T,R]=size(mixedsig);

mixedsig=reshape(mixedsig,C,[]);
Cov=mixedsig*mixedsig'/(T*R);

[u, d, ~]=svd(Cov);

d=sqrt(diag(d));
figure
bar(d);
title('Set the compression level. Default 35 biggest dimensions.')
hold on;
plot([35 35], [0, max(d)], 'k--', 'LineWidth', 1)

ylabel('Variance');
xlabel('Singular value')

ind1 = 1;
ind2=input('Last included dimension: ');

if isempty(ind2)
    ind2 = [35];
end

disp(['Compression level set to ',num2str(ind2)]);


Mwhite=diag(1./d(ind1:ind2))*u(:,ind1:ind2)'; % whitening or sphering matrix
Mdewhite=u(:,ind1:ind2)*diag(d(ind1:ind2));
P=u(:,ind1:ind2)*u(:,ind1:ind2)';

mixedsig=Mwhite*mixedsig;
N=ind2-ind1+1; 

% perform fast ICA
[A, Wwhite] = fpica(mixedsig, eye(N), eye(N), 'symm', ...
			N, 'tanh', 'off', 1, 1, 1, 'on', ...
			0.0001, maxNumIterations, 100, 'rand', ...
			eye(N), 1, 'off', 1e4, ...
			'on');
icasig= Wwhite*mixedsig; %  EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:); Wwhite = icaweights*icasphere, ica to sphere matrix
A=Mdewhite*A;
[~, isort]=sort(sum(A.^2), 'descend'); % sort the orders of components
A=A(:, isort);
Wwhite=Wwhite(isort,:);
icasig=icasig(isort,:,:);
W=Wwhite*Mwhite;

icasig=reshape(icasig,N,T,R);