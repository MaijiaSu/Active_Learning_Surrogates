clc,clear,close all
TestExample = 'eg18_UQ';                             % select the example from LSF_library, LimtStateFunction_select.m;
SBM.method = 'MCS_UQ';                            % numerical simulation method，= 'MCS'、'IS'、'SDMCS' for Pf; = 'MCS_UQ' for full distribution
SBM.iniNoS = 1e5;                                 % the size of the population of initial sampling-based-method 
SBM.NofInterval = 100;

% Metamodel
SurrModelPar.Type = 'UQLab_Kriging';
SurrModelPar.MetaOpts.Trend.Type = 'ordinary';     % 'ordinary'(default), 'linear', and 'quadratic'; 'polynomial'
SurrModelPar.MetaOpts.Corr.Type =  'Ellipsoidal';   % 'separable', 'Ellipsoidal'(default)
SurrModelPar.MetaOpts.Corr.Family = 'Gaussian';    % Linear','exponential','Gaussian', 'Matern-3_2', 'Matern-5_2'(default)

% AL strategy
ALSMPar.IniDoE.GenType = 'iniMDS';                  % = 'iniLHS','iniRandom','iniMDS'
ALSMPar.IniDoE.N0 = 15;                     
ALSMPar.LF_type = 'TwoStepLF';                      % lerarning function， = 'MoV', 'TwoStepLF','TwoStepLF_modified'
ALSMPar.LF_Par.Kernel = 'TwoStepLF_GaussianKernel'; %  'TwoStepLF_GaussianKernel','TwoStepLF_DiracKernel'
ALSMPar.Stopcon.type='SC_FPD';                      % Stop condition for full probability distribution (FPD), ='SC_FPD_Stability','SC_FPD'
ALSMPar.Stopcon.etol= 0.2;

% Self-input Candidate Pool
ALSMTimeHis.File_iniState = 'eg18_iniState';
ALSMTimeHis.false_ture_statistics = 1;
%% single run

% ALRMResult = mainALRM...
%     (TestExample,SBM,SurrModelPar,ALSMPar,ALSMTimeHis);
% 
% % PDC.PTO = [1,1,1,1,1,1,1];
% PDC.bound = [-8,-8,8,8]
% close all
% ALRM_PostProcessingV2(ALRMResult,PDC);

%%

NofRun = 10;
tic
for jj = 1:NofRun
disp(['--------------------','RUN-',num2str(jj),'--------------------------'])
    ALRMResult = mainALRM...
        (TestExample,SBM,SurrModelPar,ALSMPar,ALSMTimeHis);  
    % Save the data of interest
    Result{jj} = ALRMResult;   
    NN(jj) = ALRMResult.NofDoE;
    errorW_y{jj} = ALRMResult.ALSMTimeHis.errorW_y;
    DoE{jj} = ALRMResult.SurrModelPar.DoE;
    NofDoE(jj,1) = ALRMResult.NofDoE;
end
toc
% save data_Example17_10runs NN errorW_y



%% the Covergence trend
% load data_Example18_10runs.mat
figure
ave_errorW_y = zeros(NofRun,max(NN)-14);
for jj = 1:NofRun
    nn = 15:NN(jj);
    plot(nn,errorW_y{jj},'-.k')
    ave_errorW_y(jj,1:NN(jj)-14) = errorW_y{jj};
    hold on
end
Ave_errorW_y = [];
for ii = 1:max(NN)-14
    temp = ave_errorW_y(:,ii);
    Ave_errorW_y(ii) = mean(mean(temp(temp>0)));
end
p(1) = plot(15:max(NN),Ave_errorW_y,'--r','linewidth',2);
set(gca,'YScale','log')
xlabel('Number of DoEs')
ylabel('${{\varepsilon}_W}$')
grid on

%% the DoE distribution
% figure
% DoEY =[];
% for jj = 1:NofRun
%     DoEY = [DoEY;DoE{jj}.Y];
%     plot(16:numel(DoE{jj}.Y),DoE{jj}.Y(16:end),'-bo')
%     hold on
%     plot([16,max(NN)],[ALRMResult.ProSys.ymin,ALRMResult.ProSys.ymin],'r')
%     hold on
%     plot([16,max(NN)],[ALRMResult.ProSys.ymax,ALRMResult.ProSys.ymax],'r')
% end
figure
DoEY = [];DoEX=[];
for jj = 1:size(Result,2)
    DoEY = [DoEY;Result{1,jj}.SurrModelPar.DoE.Y];
    DoEX = [DoEX;Result{1,jj}.SurrModelPar.DoE.X];  
end
figure
histogram(DoEY,100);



%% plot the relative distance 
figure
for  jj = 1:size(Result,2)
    DoE = Result{1,jj}.SurrModelPar.DoE;
    MeanC = ProSys.muX(:)';
    NormDoEX = DoE.X./MeanC;
    D = pdist(NormDoEX, 'euclidean');
    Z = squareform(D);
    Z = Z+eye(size(Z))*max(max(Z));
    minZ = min(Z);
    p(1) = plot(1:numel(DoE.Y),minZ,'k--')
    hold on
end

%%
figure('OuterPosition',[0,0,1000,450])
ALRMResult = Result{8};
Y_true = ALRMResult.SBM.G_ture;
nmax = size(ALRMResult.ALSMTimeHis.Gpredict_his,2);
nn = ceil(linspace(1,nmax,10));
for jj = 1:numel(nn)
    subplot(2,5,jj)
    plot(Y_true,ALRMResult.ALSMTimeHis.Gpredict_his(:,nn(jj)),'.')
    alpha(0.1)
    y_max = max(Y_true);
    y_min = min(Y_true);
    hold on
    plot([y_min,y_max],[y_min,y_max],'-','Color',[178,34,34]/255,'LineWidth',2)
    %     hold on
    %     plot(DoE.Y,DoE.Y,'oc')
    ylabel('$\hat\mathcal{M}$')
    xlabel('$\mathcal{M}$')
    axis([y_min,y_max,y_min,y_max])
    title(['$N_M=',num2str(nn(jj)+15),'$'])
    grid on
end

%%
% vm =5;
% bound = [ProSys.muX-vm*ProSys.sigmaX;ProSys.muX+vm*ProSys.sigmaX]'
% density = 10;
% Myfun = @(X) ProSys.fname(X,ProSys.fun_par)
% data = PlotHighDimSurface(Myfun,bound,density)
% colormap(gca,'hot')