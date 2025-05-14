clc,clear,close all
TestExample = 'eg17_UQ';                             % select the example from LSF_library, LimtStateFunction_select.m;
SBM.method = 'MCS_UQ';                            % numerical simulation method£¬= 'MCS'¡¢'IS'¡¢'SDMCS' for Pf; = 'MCS_UQ' for full distribution
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
ALSMPar.LF_type = 'TwoStepLF';                      % lerarning function£¬ = 'MoV', 'TwoStepLF','TwoStepLF_modified'
ALSMPar.LF_Par.Kernel = 'TwoStepLF_GaussianKernel'; %  'TwoStepLF_GaussianKernel','TwoStepLF_DiracKernel'
ALSMPar.LF_Par.LearnFun = 'U';                      % LearnFunction: 'U','P'
% ALSMPar.LF_Par.c = 1e2;  
ALSMPar.Stopcon.type='SC_FPD';                      % Stop condition for full probability distribution (FPD), ='SC_FPD_Stability','SC_FPD'
ALSMPar.Stopcon.etol= 0.2;

% Self-input Candidate Pool
ALSMTimeHis.File_iniState = 'eg17_iniState';
ALSMTimeHis.false_ture_statistics = 1;
%% single run

% ALRMResult = mainALRM...
%     (TestExample,SBM,SurrModelPar,ALSMPar,,ALSMTimeHis);
% 
% PDC.PTO = [1,1,1,1,1,1,1];
% PDC.bound = [-8,-8,8,8]
% close all
% ALRM_PostProcessingV2(ALRMResult,PDC);

%%

NofRun = 15;
tic
for jj = 1:NofRun
disp(['--------------------','RUN-',num2str(jj),'--------------------------'])
    ALRMResult = mainALRM...
        (TestExample,SBM,SurrModelPar,ALSMPar,ALSMTimeHis);  
     % Save the data of interest
    Result{jj} = ALRMResult;
    DoE{jj} = ALRMResult.SurrModelPar.DoE;
    NofDoE(jj,1) = ALRMResult.NofDoE;
    errorW_y{jj} = ALRMResult.ALSMTimeHis.errorW_y;
   
    % Save the data of interest in brief
    Result_Simple.DoE{jj,1} = Result{1,jj}.SurrModelPar.DoE;
    Result_Simple.NofDoE(jj,1) = Result{1,jj}.NofDoE;
    Result_Simple.errorW_y{jj} = Result{1,jj}.ALSMTimeHis.errorW_y;
    Result_Simple.w_y{jj,1} = Result{1,jj}.ALSMTimeHis.w_y;
    Result_Simple.W_y{jj,1} =  Result{1,jj}.ALSMTimeHis.w_y;
    Result_Simple.CDF{jj,1} =  Result{1,jj}.ALSMTimeHis.CDF;
    Result_Simple.CDF1{jj,1} =  Result{1,jj}.ALSMTimeHis.CDF1;
    Result_Simple.CDF3{jj,1} =  Result{1,jj}.ALSMTimeHis.CDF3;
    Result_Simple.CDF_ture{jj,1} =  Result{1,jj}.ALSMTimeHis.CDF_ture;
end
toc

% save AL_eg18_Matern_MoV Result_Simple ALSMPar SurrModelPar
% save AL_eg18_Matern_P(C=0.01) Result_Simple ALSMPar SurrModelPar
% save AL_eg18_Matern_new_P(C=100) Result_Simple ALSMPar SurrModelPar
% save AL_eg18_Gaussian_new_P(C=100) Result_Simple ALSMPar SurrModelPar
% save AL_eg18_PCE_U Result_Simple ALSMPar SurrModelPar
% save AL_eg18_PCE_P(C=100) Result_Simple ALSMPar SurrModelPar
% save AL_eg18_PCE_MoV Result_Simple ALSMPar SurrModelPar
% save data_Example17_10runs NN errorW_y
%% the Covergence trend
% load data_Example17_15runs.mat
figure
% NN = cell2mat(NN);
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
plot(15:max(NN),Ave_errorW_y,'--r','linewidth',2)
set(gca,'YScale','log')
xlabel('Number of DoEs')
ylabel('${{\varepsilon}_W}$')
grid on

%% the DoE distribution
figure
DoEY =[];
for jj = 1:NofRun
    DoEY = [DoEY;DoE{jj}.Y];
    plot(16:numel(DoE{jj}.Y),DoE{jj}.Y(16:end),'-bo')
    hold on
    plot([16,max(NN)],[ALRMResult.ProSys.ymin,ALRMResult.ProSys.ymin],'r')
    hold on
    plot([16,max(NN)],[ALRMResult.ProSys.ymax,ALRMResult.ProSys.ymax],'r')
end
%%
% vm =5;
% bound = [ProSys.muX-vm*ProSys.sigmaX,ProSys.muX+vm*ProSys.sigmaX]
% density = 10;
% Myfun = @(X) ProSys.fname(X,ProSys.fun_par)
% data = PlotHighDimSurface(Myfun,bound,density)
% % colormap(gca,'hot')