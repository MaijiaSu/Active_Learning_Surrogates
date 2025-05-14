function DoE = IniDoE_MCMC(TargetPdf,SamplePool,NofDoE,Ndim)

% TargetPdf = ProSys.MCMC.TargetPdf;

PoolPDF = TargetPdf(SamplePool);
p1 = sort(PoolPDF,'descend');
p0 = p1(floor(0.999*size(SamplePool,1)));

logPfuns_IX = @(x) log(MyStepFun(TargetPdf,x,p0));
NofSamples = size(SamplePool,1);
ThinChain = 10;
ExpextedNofSamples = 1e4;
mccount = ExpextedNofSamples*ThinChain;  % 2 = 1/(1-50%)
seed = SamplePool(randperm(1e4,Ndim*2),:);

DoE = gwmcmc(seed',logPfuns_IX,mccount,'ThinChain',ThinChain);
DoE=DoE(:,:)';
DoE = DoE(end-NofDoE+1:end,:);

return

% clc,clear,close all
% TestExample = 'eg30';  
% LimtStateFunction_select
% ExpextedNofSamples = 1e4;
% Ndim = 2;
% ThinChain = ProSys.MCMC.ThinChain;
% mccount = ExpextedNofSamples*ThinChain*1.25;   % 1.25 = 1/(1-20%), 20% as burn-in
% NofChain = ProSys.MCMC.NofChain;               %number of Chain
% logPfuns = @(x) log(ProSys.MCMC.TargetPdf(x));
% 
% [SamplePool,logP]=gwmcmc(randn(Ndim,NofChain),logPfuns,mccount,'ThinChain',ThinChain);
% SamplePool(:,:,1:floor(size(SamplePool,3)*0.2))=[];  %remove 20% as burn-in
% SamplePool=SamplePool(:,:)';

%%
% inital DoE
% input: SamplePool, TargetPDF
% 
% NofDoE = 300;
% 
% figure





plot(SamplePool(:,1),SamplePool(:,2),'.')
hold on
plot(DoE(:,1),DoE(:,2),'r*')

hold on
bound = [-5,-2,5,8];
gap = 300;
[X1, X2, YX] = contour2D(bound,gap,ProSys.MCMC.TargetPdf);
hold on
[X1, X2, YX] = contour2D(bound,gap,@(x) MyStepFun(TargetPdf,x,p0));

%%

% PoolPDF = ProSys.MCMC.TargetPdf(SamplePool);
% p1 = sort(PoolPDF,'descend');
% p0 = p1(floor(0.95*ExpextedNofSamples));
% 
% 
% bound = [-4,-2,4,7];
% gap = 300;
% Myfun = @(x) tempfun(ProSys.MCMC.TargetPdf,x,p0);
% [X1, X2, YX] = contour2D(bound,gap,Myfun);
% 
% hold on
% [X1, X2, YX] = contour2D(bound,gap,ProSys.MCMC.TargetPdf);






