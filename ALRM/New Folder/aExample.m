clc,clear,close all
TestExample = 'eg1';                               % select the example from LSF_library, LimtStateFunction_select.m;
SBM.method = 'MCS_UQ';                             % numerical simulation method£¬= 'MCS'¡¢'IS'¡¢'SDMCS' for Pf; = 'MCS_UQ' for full distribution
SBM.iniNoS = 1e5;                                  % the size of the population of initial sampling-based-method 
SBM.NofInterval = 100;
SurrModelPar.Type = 'ooKriging';                   % surrogate model£¬= 'Kriging','ooKriging','SVM','SVR','PCE'(Bootstap),'PCE_Jackknife','PCKriging'
% SurrModelPar.opts.corrFunc = @corrmatern52;
% SurrModelPar.MetaOpts.Degree = 100;
% SurrModelPar.MetaOpts.Degree = 3:15;
ALSMPar.IniDoE.GenType = 'iniRandom';              % = 'iniLHS','iniRandom' 
ALSMPar.IniDoE.N0 = 15;                     
ALSMPar.LF_type = 'MoV';                 % lerarning function£¬ = 'MoV', 'TwoStepL_ZiQi','TwoStepL_modified'
ALSMPar.LF_Par.Kernel = 'TwoStepL_GaussianKernel'; %  'TwoStepL_GaussianKernel','TwoStepL_DiracKernel'
ALSMPar.Stopcon.type='SC_FPD';                     % Stop condition for full probability distribution (FPD), ='SC_FPD_Stability','SC_FPD'
%
% ALSMTimeHis.File_iniState = 'eg26_iniState';
ALSMTimeHis.false_ture_statistics = 1;

NofRun = 1;
tic
for jj = 1:NofRun
disp(['--------------------','RUN-',num2str(jj),'--------------------------'])
    ALRMResult = mainALRM...
        (TestExample,SBM,SurrModelPar,ALSMPar,ALSMTimeHis);  
    % Save the data of interest
    DoE{jj} = ALRMResult.SurrModelPar.DoE;
    NofDoE(jj,1) = ALRMResult.NofDoE; 
    Moment(jj,:) = ALRMResult.SBM.Moment;
    Wy_ture{jj} = ALRMResult.ALSMTimeHis.Wy_ture;
    index1{jj} = ALRMResult.ALSMTimeHis.errorCDF;
    Result{jj} = ALRMResult;    
%     index2{jj} = ALRMResult.ALSMTimeHis.absDetaG_iter;
%     Wy{jj} = ALRMResult.ALSMTimeHis.W_y;
%     DetaG_ture{jj} = ALRMResult.ALSMTimeHis.absDetaG_ture;     
end
toc

close all
 PTO = [1,1,1,1,1,1];
    ALRM_PostProcessingV2(ALRMResult,PTO);

% % average the result
% for jj = 1:NofRun
%    temp_Wy(1,jj) =  Wy_ture{jj}(end);
% end
% 
% figure
% for jj = 1:NofRun
%     plot(DoE{jj}.X(13:end,1),DoE{jj}.X(13:end,2),'k.');
%     hold on
% end
% 
% temp_Wy(isnan(temp_Wy))=[];
% [mean(temp_Wy),mean(Moment),mean(NofDoE);
%     std(temp_Wy'),std(Moment),std(NofDoE)]
%% 
% clc,clear,close all
% TestExample = 'eg18';              % select the example from LSF_library, LimtStateFunction_select.m;
% SBM.method = 'MCS';             % numerical simulation method£¬= 'MCS'¡¢'IS'¡¢'SDMCS'; = 'MCS_UQ' for full pro-distribution
% SBM.iniNoS = 1e5;
% SurrModelPar.Type = 'Kriging';        % surrogate model£¬= 'Kriging'¡¢'SVM'¡¢'SVR'¡¢'PCE'
% ALSMPar.LF_type = 'P';            % lerarning function£¬= 'U'¡¢'D'¡¢'EFF'¡¢'P'
% ALSMPar.Stopcon.type='SC_UQLab';  % termination criteria
% mainProcedure                     % Run....


