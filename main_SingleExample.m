addpath('ALRM')
addpath(genpath('Toolboxes'))
clc,clear,close all

%% Define the UQ Problem 
% #1: select the example from "LimtStateFunction_select.m";
TestExample =  'BatchProcessEg15'
BenchmarkProblemForFullDistribution

%% Select the simulation method
SBM.method = 'MCS_UQ';                            
% #1 numerical simulation methods = 'MCS','IS','SDMCS' only for computing Pf; 
% #2 numerical simulation methods = 'MCS_UQ' for computing full distribution
SBM.iniNoS = 1e5;   % the size of the population of initial sampling-based-method 
SBM.NofInterval = 100;

%% Metamodel
SurrModelPar.Type = 'Kriging'; % requries the dace toolbox
% SurrModelPar.MetaOpts.MetaType = 'Kriging';
% SurrModelPar.MetaOpts.Trend.Type = 'ordinary';     % 'ordinary'(default), 'linear', and 'quadratic'; 'polynomial'
% SurrModelPar.MetaOpts.Corr.Type = 'separable';     % 'separable', 'Ellipsoidal'(default)
% SurrModelPar.MetaOpts.Corr.Family = 'Gaussian';    % Linear','exponential','Gaussian', 'Matern-3_2', 'Matern-5_2'(default)
% SurrModelPar.MetaOpts.Corr.Isotropic = false;
% SurrModelPar.MetaOpts.EstimMethod = 'ML';

%% AL strategy
ALSMPar.IniDoE.GenType = 'iniMDS';                  % = 'iniLHS','iniRandom','iniMDS'
ALSMPar.IniDoE.N0 = 15;                             % initial size of DoE
ALSMPar.LF_type = 'TwoStepLF';                      % learning functions = 'MoV', 'TwoStepLF','TwoStepLF_modified'
ALSMPar.LF_Par.Kernel = 'TwoStepLF_GaussianKernel'; %  'TwoStepLF_GaussianKernel','TwoStepLF_DiracKernel'
ALSMPar.LF_Par.LearnFun = 'U';                      % LearnFunction: 'U','P'  
ALSMPar.Stopcon.type='SC_FPD';                      % Stop condition for full probability distribution (FPD), ='SC_FPD_Stability','SC_FPD'
ALSMPar.Stopcon.etol= 0.2;                          % stopping threshold
ALSMTimeHis.MaxNoDoE = 300; % Maximum size of DoE


%% single run
ALRMResult = mainALRM...
    ([],SBM,SurrModelPar,ALSMPar,ALSMTimeHis,ProSys);


%% Post-processing
close all
PDC.PTO = [1,1,1,1,1,1,1]; % Figure plotting control, =1,'on'
ALRM_PostProcessingV2(ALRMResult,PDC);



