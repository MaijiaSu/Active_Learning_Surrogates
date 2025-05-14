
% Self-input Candidate Pool
clc,clear,close all
TestExample = 'eg17_UQ';                             % select the example from LSF_library, LimtStateFunction_select.m;
SBM.method = 'MCS_UQ';                            % numerical simulation method，= 'MCS'、'IS'、'SDMCS' for Pf; = 'MCS_UQ' for full distribution
LimtStateFunction_select
SBM.iniNoS = 1e5;                                %
[SBM.SamplePool,~] = RNgeneratorV2...
    (ProSys.muX,ProSys.sigmaX,ProSys.Distri,SBM.iniNoS);

SBM.G_ture = ProSys.fname(SBM.SamplePool,ProSys.fun_par);

[PDF,xc] = histcounts(SBM.G_ture,100,'Normalization','pdf');
[CDF,xc] = histcounts(SBM.G_ture,100,'Normalization','cdf');
xc = (xc(1:end-1)+xc(2:end))/2;
plot(xc,PDF,'--')
figure
subplot(1,2,1)
plot(xc,CDF,'--')
subplot(1,2,2)
CCDF = 1-CDF;
plot(xc,CCDF,'--');


ymin = min(SBM.G_ture)
ymax = max(SBM.G_ture)
%% Plot the performance surface
vm =5;
bound = [ProSys.muX-vm*ProSys.sigmaX,ProSys.muX+vm*ProSys.sigmaX]
density = 10;
Myfun = @(X) ProSys.fname(X,ProSys.fun_par)
data = PlotHighDimSurface(Myfun,bound,density);

%%
% iniState.iniSBM.SamplePool = SBM.SamplePool;
% iniState.iniSBM.NofSamples = SBM.iniNoS;
% iniState.iniSBM.G_ture = SBM.G_ture;
% save eg17_iniState iniState