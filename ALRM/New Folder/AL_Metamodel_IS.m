%  MPP-IS
%% Step 1. FORM approximation to find the MPP
f=2;
error_tol=1e-3;      
X0 = ProSys.muX;
[Xi_DP,beta,DoE_IS] = RsUsingJC(ProSys.muX,ProSys.sigmaX,ProSys.Distri,...
                               f,error_tol,X0,ProSys.fname,ProSys.fun_par);
SBM.MPP = Xi_DP';
SBM.MPP_DoE = DoE_IS;

%% Step 2. Pool-based adaptive metamodel

% 0. Initialzation of SBM-ALSM
    ALSMTimeHis.times = 0;
    ALSMTimeHis.t_start = tic;
    DoE = struct('X',[],'Y',[]);
    SBM.CPI = [];
    SBM.SPI = [];
    
% 1. Genetate initial samples
    [SBM.SamplePool,SBM.Ppdf] = RNgeneratorV2...
                      (SBM.MPP,ProSys.sigmaX,ProSys.Distri,SBM.iniNoS);
    SBM.Fpdf = MyPoolJointPdf(ProSys.muX,ProSys.sigmaX,ProSys.Distri,SBM.SamplePool);                 
    SBM.NofSamples = SBM.iniNoS;  
    SBM.CPI = 1:SBM.iniNoS;         % Candidate Pool Index in SamplePool

% 2. Train the suurogate model acoording sequentially
% -------------------------------------------------------------------------
%                     Active learning surrogate model 
% -------------------------------------------------------------------------    
    SBM.CPI = 1:SBM.NofSamples;
    SBM.CPI(unique(SBM.SPI)) = [];
    [SBM,DoE,SurrModelPar,ALSMPar,ALSMTimeHis] = ...
            AL_Metamodel(ProSys,SBM,ALSMPar,DoE,SurrModelPar,ALSMTimeHis);
% -------------------------------------------------------------------------
%                                end
% -------------------------------------------------------------------------
    
% 3. Calculate Pf and Cov_Pf      
    Fail_index = find(SBM.G_predict<=0);
    Pf = 1/SBM.NofSamples *sum(SBM.Fpdf(Fail_index)./SBM.Ppdf(Fail_index));
    State = zeros(SBM.NofSamples,1); State(Fail_index,1) = 1;
    temp1 = sum(State .*(SBM.Fpdf./SBM.Ppdf).^2);
    Var_Pf = 1/SBM.NofSamples*(1/SBM.NofSamples*temp1-Pf.^2);
    Cov_Pf = sqrt(Var_Pf)/Pf;
        
% 4. Contronling the coeifficient of variation
     while Cov_Pf > SBM.cov_pf_tol

         [tempPool,~] = RNgeneratorV2...
             (SBM.MPP,ProSys.sigmaX,ProSys.Distri,SBM.IncSize);
              
         SBM.SamplePool = [SBM.SamplePool;tempPool];
         SBM.CPI = [SBM.CPI,SBM.NofSamples+1:SBM.NofSamples+SBM.IncSize];
         SBM.NofSamples = SBM.NofSamples + SBM.IncSize;
         
         str1 = ['Increasing Pool population size, IncSize=',num2str(SBM.IncSize)];
         str2 = ['NofSamples=',num2str(SBM.NofSamples)];
         disp([str1,' ',str2])
         
% -------------------------------------------------------------------------
%                     Active learning surrogate model 
% -------------------------------------------------------------------------    
    SBM.CPI = 1:SBM.NofSamples;
    SBM.CPI(unique(SBM.SPI)) = [];
    [SBM,DoE,SurrModelPar,ALSMPar,ALSMTimeHis] = ...
            AL_Metamodel(ProSys,SBM,ALSMPar,DoE,SurrModelPar,ALSMTimeHis);
% -------------------------------------------------------------------------
%                                end
% -------------------------------------------------------------------------
        
         Fail_index = find(SBM.G_predict<=0);
         Pf = 1/SBM.NofSamples *sum(SBM.Fpdf(Fail_index)./SBM.Ppdf(Fail_index));
         State = zeros(SBM.NofSamples,1); State(Fail_index,1) = 1;
         temp1 = sum(State .*(SBM.Fpdf./SBM.Ppdf).^2);
         Var_Pf = 1/SBM.NofSamples*(1/SBM.NofSamples*temp1-Pf.^2);
         Cov_Pf = sqrt(Var_Pf)/Pf;
         
         % whether exceed the maxium Pool population size
         if SBM.NofSamples >= SBM.MaxPoolSize
            break
         end
     end