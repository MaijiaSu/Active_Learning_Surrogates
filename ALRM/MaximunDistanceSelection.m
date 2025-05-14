function [DoEX,ID] = MaximunDistanceSelection(CandidatePool,NofDoE,StartPoint)
% normalizaiton method

    NormalizationMethod = 3;

if NormalizationMethod == 1
    MeanC = mean(CandidatePool);
    stdC = std(CandidatePool);
    CandidatePool_norm = (CandidatePool-MeanC)./stdC;
elseif NormalizationMethod == 2
    MeanC = mean(CandidatePool);
    CandidatePool_norm = (CandidatePool./MeanC);
elseif NormalizationMethod == 3
    maxC = max(CandidatePool);
    minC = min(CandidatePool);
    CandidatePool_norm = CandidatePool./(maxC-minC);
else
    CandidatePool_norm = CandidatePool;
end

if nargin == 2
    DoEX_norm = [];
    NEP = size(CandidatePool,1);
    ID = zeros(1,NofDoE);
    ID(1) = randperm(NEP,1);
    DoEX_norm(1,:) = CandidatePool_norm(ID(1),:);
    distance = [];
elseif nargin == 3
    DoEX_norm = [];
    tempDist = pdist2(CandidatePool,StartPoint,'euclidean');
    [~,id] = min(tempDist);
    ID(1) = id;
    DoEX_norm(1,:) = CandidatePool_norm(ID(1),:);
    distance = [];
end

for ii = 2:NofDoE
    tempDist = pdist2(CandidatePool_norm,DoEX_norm(ii-1,:),'euclidean');
    distance = [distance,tempDist];
    if ii == 2
       [~,id] = max(distance);
    else
        [~,id] = max(min(distance')');
    end
    
    DoEX_norm(ii,:) = CandidatePool_norm(id(1),:);
    ID(ii) = id;
end

% Cancel normalizaiton
if NormalizationMethod == 1
    DoEX = (ones(NofDoE,1)*stdC).*DoEX_norm+...
                     ones(NofDoE,1).*MeanC;
elseif NormalizationMethod == 2
    DoEX = (ones(NofDoE,1)*MeanC).*DoEX_norm;
elseif NormalizationMethod == 3
    DoEX = (ones(NofDoE,1)*(maxC-minC)).*DoEX_norm;
else
    DoEX = DoEX_norm;
end

return

% Maximun distance selction
% clc,clear,close all
% TestExample = 'eg30';  
% LimtStateFunction_select
% NEP = 1e5;
% Ndim = 2;
% ThinChain = ProSys.MCMC.ThinChain;
% mccount = NEP*ThinChain*1.25;   % 1.25 = 1/(1-20%), 20% as burn-in
% NofChain = ProSys.MCMC.NofChain;               %number of Chain
% logPfuns = @(x) log(ProSys.MCMC.TargetPdf(x));
% [CandidatePool,logP]=gwmcmc(randn(Ndim,NofChain),logPfuns,mccount,'ThinChain',ThinChain);
% CandidatePool(:,:,1:floor(size(CandidatePool,3)*0.2))=[];  %remove 20% as burn-in
% CandidatePool=CandidatePool(:,:)';

% DoE = DoE(end-NofDoE+1:end,:);
% plot(CandidatePool(:,1),CandidatePool(:,2),'.')
% hold on
% plot(DoE(:,1),DoE(:,2),'k*')

% hold on
% bound = [-5,-2,5,8];
% gap = 300;
% [X1, X2, YX] = contour2D(bound,gap,ProSys.MCMC.TargetPdf);


