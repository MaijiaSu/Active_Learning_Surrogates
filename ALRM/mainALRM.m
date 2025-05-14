function ALRMResult = mainALRM(TestExample,SBM,SurrModelPar,ALSMPar,ALSMTimeHis,ProSys)
%% 
% active learning reliability method combining ...
% Sampling-based method and Active learing surrogate model
% Completed on 2019.12.12 in Guangzhou university, vesion 1.0
% Completed on 2020.07.13 in Guangzhou university, vesion 2.0
% Completed on 2021.03.21 in Guangzhou university, vesion 3.0
% Completed on 2021.08.02 in Guangzhou university, vesion 4.0
% --------------------- Parmeters list ------------------------------------
% 1. Probabilitic System Par.         
% 2. sampling-based method Par.
% 3. Surrogate model Par.
% 4. Active learning surrogate model Par.
% 4.1. IniDoE£ºInitial DoE
% 4.2. LearningFun£ºLearning Function
% 4.3. StopCon: Stoping condition
% 5. other Par.
% -------------------------------------------------------------------------

%% Default Parameter
% default Probabilitic system: TestExample = 'eg1', FBSS (Four Branch series system)
% default sampling-baed method: MCS
% default surrogate model:Kriging
% default ALSM Paramete  
% 1. iniDoE(LHS,N0=12)
% 2. LearnFun of Kring: U; LearnFun of SVR : D
% 3. stopcon: U¡úminU,EFF¡úmaxEFF,D¡úSC_UQLab

%% 1. Probabilitic system Par
if ~exist('ProSys')
    if ~exist('TestExample')
        TestExample = 'eg1';
    end
    try
       LimtStateFunction_select
    catch
       addpath('LSF_Library')
       error('Please make sure the LSF_Library is placed on the current path')
    end
end

%% 2. Sampling-baed method Parameter
% Paremeter of MCS
% the default simulation method is set to 'MCS'

if ~exist('SBM')
    SBM = [];
end

if ~isfield(SBM,'method')
    SBM.method = 'MCS';  
end

% Set the default value and check the the input Par if exist
if strcmp(SBM.method,'MCS')||strcmp(SBM.method,'MCS_UQ')
    DefaultPar = inputParser;
    DefaultPar.KeepUnmatched = 1;
    DefaultPar.addParameter('method','MCS');
    DefaultPar.addOptional('cov_pf_tol',0.05,@isnumeric);
    if strcmp(SBM.method,'MCS_UQ')
        DefaultPar.addOptional('NofInterval',100,@isnumeric);
    end
    if ProSys.Ndim <= 2
        DefaultPar.addOptional('iniNoS',1e4,@isnumeric);
        DefaultPar.addOptional('IncSize',1e4,@isnumeric);
    else  
        DefaultPar.addOptional('iniNoS',1e5,@isnumeric);
        DefaultPar.addOptional('IncSize',1e5,@isnumeric);
    end
    DefaultPar.addOptional('MaxPoolSize',1e6,@isnumeric);
    parse(DefaultPar,SBM);
    SBM = DefaultPar.Results;

elseif strcmp(SBM.method,'IS')
    DefaultPar = inputParser;
    DefaultPar.KeepUnmatched = 1; % allow the SBM with unuseful input
    DefaultPar.addParameter('method','IS');
    DefaultPar.addOptional('cov_pf_tol',0.05,@isnumeric);
    DefaultPar.addOptional('iniNoS',1*1e4,@isnumeric);
    DefaultPar.addOptional('IncSize',1e4,@isnumeric);
    DefaultPar.addOptional('MaxPoolSize',1e6,@isnumeric);
    parse(DefaultPar,SBM,'cov_pf_tol',300);
    SBM = DefaultPar.Results;
    
elseif strcmp(SBM.method,'SDMCS')
    DefaultPar = inputParser;
    DefaultPar.KeepUnmatched = 1;
    DefaultPar.addParameter('method','SDMCS');
    DefaultPar.addOptional('cov_pf_tol',0.05,@isnumeric);
    DefaultPar.addOptional('MaxPoolSize',1e6,@isnumeric);
    DefaultPar.addOptional('detaNs_max',1e4,@isnumeric);
    DefaultPar.addOptional('detaNs_min',1e4,@isnumeric);
    DefaultPar.addOptional('P_Divde_CDF',[0,1 - 10.^-(1:18)],@isnumeric);   % Spherical ring decompostion 
    parse(DefaultPar,SBM);
    SBM = DefaultPar.Results;    
end
 
%% 3. Surrogate model 
if ~exist('SurrModelPar')
    SurrModelPar = [];
end
if ~isfield(SurrModelPar,'Type')
    SurrModelPar.Type = 'Kriging';
end

% Transform the required data format of Uqlab
if strcmp(SurrModelPar.Type,'PCE_UQLAB')||...
        strcmp(SurrModelPar.Type,'PCE')||...
        strcmp(SurrModelPar.Type,'PCKriging')
    for ii = 1:  ProSys.Ndim
        InputOpts.Marginals(ii).Name = ['x',num2str(ii)]; % resistance variable
        if  ProSys.Distri(ii)== 1|| ProSys.Distri(ii)==6
            InputOpts.Marginals(ii).Type = 'Gaussian';
        elseif strcmp(ProSys.Distri,'mvnrnd')
            InputOpts.Marginals(ii).Type = 'Gaussian';
        elseif ProSys.Distri(ii)== 2
            InputOpts.Marginals(ii).Type = 'Lognormal';
        elseif ProSys.Distri(ii)== 3
            InputOpts.Marginals(ii).Type = 'Gumbel';
        elseif ProSys.Distri(ii)== 4
            InputOpts.Marginals(ii).Type = 'Gumbelmin';
        elseif ProSys.Distri(ii)== 5
            InputOpts.Marginals(ii).Type = 'Uniform';
        elseif  ProSys.Distri(ii)== 1|| ProSys.Distri(ii)==7
            InputOpts.Marginals(ii).Type = 'Weibull';
        elseif strcmp(ProSys.Distri,'mvnrnd')||strcmp(ProSys.Distri,'MCMC')
            InputOpts.Marginals(ii).Type = 'Gaussian';
        else
            error('Need to update the code')
        end
        InputOpts.Marginals(ii).Moments = [ProSys.muX(ii),ProSys.sigmaX(ii)];
    end
    myInput = uq_createInput(InputOpts);
end

% 3.1 Surrogate model Parameter
% DefaultPar_Surrogate_model _setting

% 3.2 Surrogate model calling function
%  Unify the calling method of surrogate model
%  SurrogateModel = MoedlTrain(DoE,G) 
%  [G_predict,Gmse] = ModelPredictor(SurrogateModel,PredictPool)
[SurrModelPar.MoedlTrain,SurrModelPar.ModelPredictor] = SurrogateModelSlect(SurrModelPar,ProSys);

%% 4. Active learning surrogate model Parameter
% 4.1 Initial DoE
% ALSMPar.IniDoE.GenType = 2;  % Method of generating initial Doe
                          % #1£ºiniRandom -- select from candidate Pool randomly
                          % #2£ºiniLHS -- Latin hypercubic sampling £¨LHS£©
                          % #3£ºiniPMC -- Pseudo Markov Chains
% Initialization
if ~exist('ALSMPar')
  ALSMPar = struct();
end
if ~isfield(ALSMPar,'IniDoE')
    ALSMPar.IniDoE = [];
end
% Default method of generating initial DoE                  
if  ~isfield(ALSMPar.IniDoE,'GenType')
    ALSMPar.IniDoE.GenType = 'iniMDS';
end
% Default parameter of different method                           
if strcmp(ALSMPar.IniDoE.GenType, 'iniRandom')
    if ~isfield(ALSMPar.IniDoE,'N0')
        ALSMPar.IniDoE.N0 = max(12,ProSys.Ndim+2);
    end    
elseif strcmp(ALSMPar.IniDoE.GenType, 'iniLHS')
    DefaultPar = inputParser;
    DefaultPar.addParameter('GenType','iniLHS');
    % Number of initial DoE 
    DefaultPar.addOptional('N0',max(12,ProSys.Ndim),@isnumeric);
    % Variance magnification of LHS 
    DefaultPar.addOptional('VM',3,@isnumeric);                  
    parse(DefaultPar,ALSMPar.IniDoE);
    ALSMPar.IniDoE = DefaultPar.Results;      
elseif strcmp(ALSMPar.IniDoE.GenType, 'iniPMC')
    DefaultPar = inputParser;
    DefaultPar.addParameter('GenType','iniPMC')
    DefaultPar.addOptional('VM',5,@isnumeric);                              % Generate the first state of each chain
    DefaultPar.addOptional('PMC_N0',15,@isnumeric);                         % Number of chains
    DefaultPar.addParameter('PMC_LF_type','D',@(x)any(validatestring(x,...  % Type of learing function to find out local best point 
		{'U','EFF','D'})));
    DefaultPar.addOptional('N0',15,@isnumeric);
    DefaultPar.addOptional('PMC_MaxL',5,@isnumeric);                        % Maximum length of sigle chain
    DefaultPar.addOptional('PMC_VM2',0.2,@isnumeric);                       % Variance magnification of local seaching region
    parse(DefaultPar,ALSMPar.IniDoE);
    ALSMPar.IniDoE = DefaultPar.Results;   
end

% 4.2 Learning Function
% = 'U'(default); = 'D'; = 'EFF';
if ~isfield(ALSMPar,'LF_type')
    if strcmp(SurrModelPar.Type,'Kriging')
        ALSMPar.LF_type = 'U';
    elseif strcmp(SurrModelPar.Type,'SVM')
        ALSMPar.LF_type = 'D';
    elseif strcmp(SurrModelPar.Type,'SVR')
        ALSMPar.LF_type = 'D';
    elseif strcmp(SBM.method,'MCS_UQ')
        ALSMPar.LF_type = 'TwoStepLF';
    end
else 
    if ~any(validatestring(ALSMPar.LF_type,{'U','EFF','D','UD','P','MoV',...
             'TwoStepLF','TwoStepLF_modified','MoV-distance','MoV-composite','MoV-Gradient'}))
        error('Please input correct learning function')
    end
end
% Parameters of learning function 
if ~isfield(ALSMPar,'LF_Par')
    ALSMPar.LF_Par = struct();
end
if strcmp(ALSMPar.LF_type,'P')
    DefaultPar = inputParser;
    DefaultPar.KeepUnmatched = 1;
    if ProSys.Ndim <= 2
        DefaultPar.addOptional('c',1e2,@isnumeric);
    else
        DefaultPar.addOptional('c',1e3,@isnumeric);
    end    
    DefaultPar.addOptional('alpha',0.5,@isnumeric);
    parse(DefaultPar,ALSMPar.LF_Par);
    ALSMPar.LF_Par = DefaultPar.Results; 
end
if strcmp(ALSMPar.LF_type,'TwoStepLF')
    if ~isfield(ALSMPar.LF_Par,'Kernel')
        ALSMPar.LF_Par.Kernel = 'TwoStepLF_GaussianKernel';
        %     ALSMPar.LF_Par.Kernel = 'TwoStepLF_DiracKernel';
    end
    if ~isfield(ALSMPar.LF_Par,'LearnFun')
        ALSMPar.LF_Par.LearnFun = 'P';
    end
    % c only useful for P Learning function
    if ~isfield(ALSMPar.LF_Par,'c')
        ALSMPar.LF_Par.c = 1e2;
    end
end


% 3.3 Stopping condition
% Default Stopping condition type 
% #1: minU : minU > 2
% #2: maxEFF : maxEFF < 0.001
% #3: SC_UQLab: SC1<0.001&SC2<0.001
% #4: CompositeSC :(minU > 2)||(SC1<0.001&SC2<0.001) !!{Composite condition of #1 and #3}

% Stopping threshold
% Initialization
if ~isfield(ALSMPar,'Stopcon')
    ALSMPar.Stopcon = [];
end
if ~isfield(ALSMPar.Stopcon,'type')
    if (strcmp(ALSMPar.LF_type,'U'))
        ALSMPar.Stopcon.type = 'minU';
    elseif (strcmp(ALSMPar.LF_type,'EFF'))
        ALSMPar.Stopcon.type = 'maxEFF';
    elseif (strcmp(ALSMPar.LF_type,'D'))
        ALSMPar.Stopcon.type = 'SC_UQLab';
    end
end

% Threshold of U
if (strcmp(ALSMPar.Stopcon.type,'minU'))
    ALSMPar.Stopcon.minU = 2;
% Threshold of EFF
elseif (strcmp(ALSMPar.Stopcon.type,'maxEFF'))
    ALSMPar.Stopcon.maxEFF = 0.001;
% Threshold of stopping condition from UQLab
elseif strcmp(ALSMPar.Stopcon.type,'SC_UQLab')
    ALSMPar.Stopcon.SC1 = 1e-3;    % threshold#1
    ALSMPar.Stopcon.SC2 = 1e-3;    % threshold#2
% Threshold of CompositeSC
elseif strcmp(ALSMPar.Stopcon.type,'CompositeSC')
    ALSMPar.Stopcon.minU = 2;        % default:minU = 2
    ALSMPar.Stopcon.SC1 = 0.001;    % threshold#1
    ALSMPar.Stopcon.SC2 = 0.001;    % threshold#2
elseif strcmp(ALSMPar.Stopcon.type,'Modifie_SC_UQLab')
    ALSMPar.Stopcon.SC1 = 0.001;    % threshold#1
    ALSMPar.Stopcon.SC2 = 0.01;     % threshold#2
elseif strcmp(ALSMPar.Stopcon.type,'SC_ASVM-MCS')
    ALSMPar.Stopcon.SC1 = 1e-5;    % threshold#1
    ALSMPar.Stopcon.SC2 = 1e-4;     % threshold#2
elseif strcmp(ALSMPar.Stopcon.type,'SC_FPD')  % stopcon for full probability distribution
    if ~isfield(ALSMPar.Stopcon,'etol')
        ALSMPar.Stopcon.etol = 0.02;    %
    end
elseif strcmp(ALSMPar.Stopcon.type,'SC_FPD_Stability')
    if ~isfield(ALSMPar.Stopcon,'etol')
        ALSMPar.Stopcon.etol = 0.001;    %
    end
end

%% 5. Other Parameter
if ~exist('ALSMTimeHis')
    ALSMTimeHis = struct();
end
DefaultPar = inputParser;
DefaultPar.KeepUnmatched = 1; % allow the struct with unuseful input
DefaultPar.addParameter('File_iniState','lack');                % control whether intialize the DoE or SBM
if strcmp(TestExample,'eg14')||...
        strcmp(TestExample,'eg17')||...
        strcmp(TestExample,'eg18')||...
         strcmp(TestExample,'eg19')||...
         strcmp(TestExample,'eg20')||...
        strcmp(TestExample,'eg21')||...
        strcmp(TestExample,'eg22')||...
        strcmp(TestExample,'eg23')
    DefaultPar.addOptional('false_ture_statistics',0,@isnumeric);
else
    DefaultPar.addOptional('false_ture_statistics',1,@isnumeric);
end
DefaultPar.addOptional('ProcessDisplay',1,@isnumeric);           % =1,display; otherwise, not;
DefaultPar.addOptional('Maximum_runing_time',10*3600,@isnumeric);   % 
DefaultPar.addOptional('MaxNoDoE',500,@isnumeric);               % Maxium number of DoE;
parse(DefaultPar,ALSMTimeHis);
ALSMTimeHis = DefaultPar.Results;
     
%% Parameter setting check

% ......................................

%% ---------------------------------------------------------------
%                 Active learning reliability method  
% -------------------------------------------------------------------------
text_ALRM = [];
text_ALRM{1} = ['AL','-',SurrModelPar.Type,'-',SBM.method];
text_ALRM{2} = ['Method of generating initial Doe:',ALSMPar.IniDoE.GenType];
text_ALRM{3} = [ 'Learning Function:',ALSMPar.LF_type];
text_ALRM{4} = ['StopCondition:',ALSMPar.Stopcon.type];
disp(['---------------------------',text_ALRM{1},'----------------------------'])
    fprintf([text_ALRM{2},'\n',text_ALRM{3},'\n',text_ALRM{4},'\n'])
disp('---------------------------------------------------------------------')
% -------------------------------------------------------------------------
if strcmp(SBM.method,'MCS')    
    AL_Metamodel_MCS
elseif strcmp(SBM.method,'SDMCS')
    AL_Metamodel_SDMCS
elseif strcmp(SBM.method,'IS')
    AL_Metamodel_IS
elseif strcmp(SBM.method,'MCS_UQ')
    AL_Metamodel_FullDistribution    
end
     SimTime = toc(ALSMTimeHis.t_start);
% -------------------------------------------------------------------------
%                                 END  
% -------------------------------------------------------------------------
%% Storage of Result
ALRMResult.SBM = SBM;
ALRMResult.SimTime = SimTime;
ALRMResult.NofDoE = size(DoE.X,1);
ALRMResult.ALSMTimeHis = ALSMTimeHis;
ALRMResult.SurrModelPar = SurrModelPar;
ALRMResult.ProSys = ProSys;
ALRMResult.ALSMPar = ALSMPar;
ALRMResult.text_ALRM = text_ALRM;

%%  PostProcessing of Pf estimation   
if strcmp(SBM.method,'MCS')||strcmp(SBM.method,'IS')||strcmp(SBM.method,'SDMCS')
% PostProcessing1
     disp(['Pf=',num2str(Pf),' ',...
           'SBM.NofSamples= ',num2str(SBM.NofSamples),' ',...
           'NofDoE= ',num2str(size(DoE.X,1))])
end

%% PostProcessing of CDF/CCDF estiamation 
if strcmp(SBM.method,'MCS_UQ')
% PostProcessingV2
    disp(['SBM.NofSamples=',num2str(SBM.NofSamples),' ',...
        'NofDoE=',num2str(size(DoE.X,1)),' ',...
        'Error=',num2str(error_Wy)])
end
return