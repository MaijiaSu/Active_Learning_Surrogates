%% 
clc,clear,close all
TestExample = 'eg1';              % select the example from LSF_library, LimtStateFunction_select.m;
SBM.method = 'MCS';             % numerical simulation method£¬= 'MCS'¡¢'IS'¡¢'SDMCS'; = 'MCS_UQ' for full pro-distribution
SBM.iniNoS = 1e5;
SurrModelPar.Type = 'Kriging';        % surrogate model£¬= 'Kriging'¡¢'SVM'¡¢'SVR'¡¢'PCE'
ALSMPar.LF_type = 'P';            % lerarning function£¬= 'U'¡¢'D'¡¢'EFF'¡¢'P'
ALSMPar.Stopcon.type='SC_UQLab';  % termination criteria
% ALSMTimeHis.File_iniState = 'eg26_iniState';
ALSMTimeHis.false_ture_statistics = 1;
ALRMResult = mainALRM...
        (TestExample,SBM,SurrModelPar,ALSMPar,ALSMTimeHis); 
 
close all
 PTO = [1     1     1     1];
 ALRM_PostProcessingV1(ALRMResult,PTO)   