%  MCS
% 0. Initialzation of SBM-ALSM
    ALSMTimeHis.times = 0;
    ALSMTimeHis.t_start = tic;
    DoE = struct('X',[],'Y',[]);
    SBM.CPI = [];
    SBM.SPI = [];
    
% 1. Genetate initial samples        
    [SBM.SamplePool,~] = RNgeneratorV2...
                      (ProSys.muX,ProSys.sigmaX,ProSys.Distri,SBM.iniNoS);
    SBM.NofSamples = SBM.iniNoS;  
    SBM.CPI = 1:SBM.iniNoS;         % Candidate Pool Index in SamplePool

    
% 2. Train the suurogate model acoording sequentially   
% -------------------------------------------------------------------------
%                     Active learning surrogate model 
% -------------------------------------------------------------------------    
        [SBM,DoE,SurrModelPar,ALSMPar,ALSMTimeHis] = ...
            AL_Metamodel(ProSys,SBM,ALSMPar,DoE,SurrModelPar,ALSMTimeHis);
% -------------------------------------------------------------------------
%                                end
% -------------------------------------------------------------------------    

% 3. Calculate Pf and Cov_Pf         
    Pf = SBM.Pf_predict;
    Cov_Pf = sqrt((1-Pf)/Pf/SBM.NofSamples);
 
% 4. Contronling the coeifficient of variation
     while Cov_Pf > SBM.cov_pf_tol

         [tempPool,~] = RNgeneratorV2...
             (ProSys.muX,ProSys.sigmaX,ProSys.Distri,SBM.IncSize);
     
         SBM.SamplePool = [SBM.SamplePool;tempPool];
         SBM.CPI = [SBM.CPI,SBM.NofSamples+1:SBM.NofSamples+SBM.IncSize];
         SBM.NofSamples = SBM.NofSamples + SBM.IncSize;

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
               
         str1 = ['Increasing Pool population size, IncSize=',num2str(SBM.IncSize)];
         str2 = ['NofSamples=',num2str(SBM.NofSamples)];
         disp(['MCS:',str1,' ',str2])
         
         Pf = SBM.Pf_predict;
         Cov_Pf = sqrt((1-Pf)/Pf/SBM.NofSamples);
         
         % whether exceed the maxium Pool population size
         if SBM.NofSamples >= SBM.MaxPoolSize
            break
         end
     end
     
%% Sensitivity Analysis
% Indicator function
tempID = find(SBM.G_predict<=0);
Ix_pre = zeros(SBM.NofSamples,1);
Ix_pre(tempID) = 1;
SBM.Pf_sensi = zeros(2,ProSys.Ndim);

for nn = 1:ProSys.Ndim
    SBM.Pf_sensi(1,nn) = sum(Ix_pre.*(SBM.SamplePool(:,nn)-ProSys.muX(nn))...
                           /ProSys.sigmaX(nn)^2)/SBM.NofSamples;
end

for nn = 1:ProSys.Ndim
    temp = 1/ProSys.sigmaX(nn)*...
           (((SBM.SamplePool(:,nn)-ProSys.muX(nn))/ProSys.sigmaX(nn)).^2-1);
    SBM.Pf_sensi(2,nn) = sum(Ix_pre.*temp)/SBM.NofSamples;
end