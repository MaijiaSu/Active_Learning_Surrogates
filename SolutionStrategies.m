% This file defines the 12 strategies used in the paper, inlcuding 3 
% surrogates models and 4 enrichment methods

SurrModelPar = [];
ALSMPar = [];
%% Metamodel
if MetamodelStype == 1
    % GP
    SurrModelPar.Type = 'UQLab_Kriging';
    SurrModelPar.MetaOpts.Trend.Type = 'ordinary';       % 'ordinary'(default), 'linear', and 'quadratic'; 'polynomial'
    SurrModelPar.MetaOpts.Corr.Type = 'separable';     % 'separable', 'Ellipsoidal'(default)
    SurrModelPar.MetaOpts.Corr.Family = 'Gaussian';    % Linear','exponential','Gaussian', 'Matern-3_2', 'Matern-5_2'(default)
    SurrModelPar.MetaOpts.Corr.Isotropic = false;
    SurrModelPar.MetaOpts.EstimMethod = 'ML';
elseif MetamodelStype == 2
    % PCE
    SurrModelPar.Type = 'PCE';
    SurrModelPar.MetaOpts.TruncOptions.qNorm = 0.80;
    SurrModelPar.MetaOpts.Bootstrap.Replications = 100;
    SurrModelPar.MetaOpts.Degree = 1:20;
    SurrModelPar.MetaOpts.Method = 'LARS';
elseif MetamodelStype == 3
    % PC-Kriging
    SurrModelPar.Type = 'PCKriging';
    SurrModelPar.MetaOpts.MetaType = 'PCK';
    SurrModelPar.MetaOpts.Kriging.Corr.Family = 'Gaussian';
    SurrModelPar.MetaOpts.Kriging.Corr.Family = 'separable';
    SurrModelPar.MetaOpts.Kriging.Corr.Isotropic = false;
    SurrModelPar.MetaOpts.Kriging.EstimMethod = 'ML';
    SurrModelPar.MetaOpts.Kriging.Trend.Type = 'ordinary';
    SurrModelPar.MetaOpts.Mode = 'sequential';               % ='optimal','sequential'
    SurrModelPar.MetaOpts.PCE.Method = 'LARS';
    SurrModelPar.MetaOpts.PCE.TruncOptions.qNorm = 0.80;
    SurrModelPar.MetaOpts.PCE.Degree = 1:3;
end
%% AS1: AL strategy-experimental design
% Initial experimental design
if AS1 == 1
    ALSMPar.IniDoE.GenType = 'iniRandom';                  % = 'iniLHS','iniRandom','iniMDS'

elseif  AS2 == 2
    ALSMPar.IniDoE.GenType = 'iniMDS';   
end
 

%% AS2: Learning function
% Learning function
if AS2 == 1
    ALSMPar.LF_type = 'MoV';
elseif AS2 == 2
    ALSMPar.LF_type = 'TwoStepLF';                      % lerarning function£¬ = 'MoV', 'TwoStepLF','TwoStepLF_modified'
    ALSMPar.LF_Par.Kernel = 'TwoStepLF_GaussianKernel'; %  'TwoStepLF_GaussianKernel','TwoStepLF_DiracKernel'
    ALSMPar.LF_Par.LearnFun = 'U';                      % LearnFunction: 'U','P'
    % ALSMPar.LF_Par.c = 1e2; 
elseif AS2 == 3
    ALSMPar.LF_type = 'MoV-Gradient';     
elseif AS2 == 4
    ALSMPar.LF_type = 'MoV-distance'
end

%% Stop criteria
% ALSMPar.Stopcon.type='SC_FD_composite';                      % Stop condition for full probability distribution (FPD), ='SC_FPD_Stability','SC_FPD','SC_FD_composite'
% ALSMPar.Stopcon.etol= 0.0001;
ALSMPar.Stopcon.type='SC_FD_composite';                      % Stop condition for full probability distribution (FPD), ='SC_FPD_Stability','SC_FPD','SC_FD_composite'
ALSMPar.Stopcon.etol= 0.2;

%% Other;
ALSMTimeHis.MaxNoDoE = min(100+ProSys.Ndim*20,300);
ALSMPar.IniDoE.N0 = max(12,ProSys.Ndim);


%% ------------------------------------------------------------------------
% -------------------------------------------------------------------------
% LS1 = [1,2,3];
% LS2 = [1,2];
% LS3 = [1,2,3];
% tt = 0;
% ID = []
% for ii = 1:numel(LS1)
%     for jj = 1:numel(LS2)
%         for kk = 1:numel(LS3)
%             tt = tt+1;
%             temp = [tt,ii,jj,kk];
%             ID = [ID;temp]
%         end
%     end
% end





