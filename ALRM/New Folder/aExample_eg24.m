clc,clear,close all
TestExample = 'eg24';                              % select the example from LSF_library, LimtStateFunction_select.m;
SBM.method = 'MCS_UQ';                    % numerical simulation method£¬= 'MCS'¡¢'IS'¡¢'SDMCS' for Pf; = 'MCS_UQ' for full distribution
LimtStateFunction_select
SBM.iniNoS = 1e5;                                 % the size of the population of initial sampling-based-method 
SBM.NofInterval = 100;

% Metamodel
SurrModelPar.Type = 'PCE';

% AL strategy
ALSMPar.IniDoE.GenType = 'iniMDS';              % = 'iniLHS','iniRandom','iniMDS'
ALSMPar.IniDoE.N0 = max(12,ProSys.Ndim);                     
ALSMPar.LF_type = 'MoV-distance';                 % lerarning function£¬ = 'MoV', 'TwoStepLF','TwoStepL_modified','MoV-distance'
ALSMPar.LF_Par.Kernel = 'TwoStepLF_GaussianKernel'; %  'TwoStepLF_GaussianKernel','TwoStepLF_DiracKernel'
ALSMPar.Stopcon.type='SC_FPD';                     % Stop condition for full probability distribution (FPD), ='SC_FPD_Stability','SC_FPD'
ALSMPar.Stopcon.etol= 0.2;

% Other Setting
ALSMTimeHis.MaxNoDoE = 100+ProSys.Ndim*50;

%% singal run
ALRMResult = mainALRM...
    (TestExample,SBM,SurrModelPar,ALSMPar);

% PTO = [1,1,1,1,1,1];
PDC = [];
ALRM_PostProcessingV2(ALRMResult,PDC);

