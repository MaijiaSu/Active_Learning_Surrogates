% SDMCS
% 0. Initialzation of SBM-ALmetamodel
    ALSMTimeHis.times = 0;
    ALSMTimeHis.t_start = tic;
    DoE = struct('X',[],'Y',[]);
    SBM.CPI = [];
    SBM.SPI = [];
    
% 1. Initializtion of SDMCS   
    SBM.NoR = 1;                    % Novber of Ring
    n = 1;                          % the current ring which decides to samplig in     
    SBM.SubPoolIndex{1} = [];       % Index of Sub-Sample-Pool in the whole SamplePool
    SBM.NoS_ring(1) = 0;            % Nov of samples in each ring
    SBM.SamplePool = [];
    SBM.NofSamples = 0;             % the whole samples pool size
    SBM.G_ture = [];
    % JoinOption = 1;               % A option that control whether joining the next subspace
    detaNs = SBM.detaNs_max;  
    
while 1   
% 2. Spherical decomposition of space and sampling in subregion Dm-1
    % Divide the outer spherical ring into two subregions
    SBM.NoR = SBM.NoR+1;            % Novber of Ring
    n = n +1;
    p1 = SBM.P_Divde_CDF(n-1);p2 = SBM.P_Divde_CDF(n);
    SBM.thtai(n-1) = p2 - p1;       % Probability volume of subregion Dm-1
    SBM.thtai(n) = 1-p2;    
    % Sampling in subregion D(m-1)
    tempPool = SDMCS_RVGenerator(p1,p2,ProSys.Ndim,detaNs);
    SBM.SamplePool = [SBM.SamplePool;tempPool];
    SBM.SubPoolIndex{n-1} = SBM.NofSamples+1:SBM.NofSamples+detaNs;
    SBM.NofSamples = SBM.NofSamples+detaNs;
    SBM.NoS_ring(n-1) = detaNs;
    SBM.CPI = [SBM.CPI,SBM.SubPoolIndex{n-1}];
    str1 = ['Generate ',num2str(detaNs),' samples in subregion D',num2str(n-1)];
    str2 = [',Current NofSamples=',num2str(SBM.NofSamples)];
    disp(['SDMCS:',str1,' ',str2])
    
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

% 3. Estimation of the current failure probability
    % Calculate conditional probability  
    % Pf equals to the sum of subfailure probabilities in the first m-1 subregions
    for ni = 1:SBM.NoR-1
        SubPool_response = SBM.G_predict(SBM.SubPoolIndex{ni},1);
        nf_SubPool = length(find(SubPool_response<=0));
        Pi(ni) = nf_SubPool/ SBM.NoS_ring(ni);
        PFi(ni) = Pi(ni) * SBM.thtai(ni);
    end
    % Current fialure probability
    Pf =  sum(PFi);
    
% 4. Stop condition on spherical decomposition of space.   
    N_MCS = ceil((1-Pf)/Pf/SBM.cov_pf_tol^2 * SBM.thtai(n)); 
    if   N_MCS > SBM.detaNs_min
        detaNs = SBM.detaNs_max;      
    else
    % the N_MCS samples are sufficient to represent the outermost spherical
    % ring, the ring number m is fixed to its current value and the method
    % goes to Step 5.
        SBM.P_Divde = SBM.P_Divde_CDF(1:SBM.NoR);
        SBM.P_Divde(SBM.NoR+1) = 1;
        break
    end
end   

% 5. Sampling in the last spherical ring Dm
    p1 = SBM.P_Divde_CDF(n);
    p2 = SBM.P_Divde_CDF(n+1);
    SBM.SubPoolIndex{n} = [];
    SBM.NoS_ring(n) = 0;
    detaNs = N_MCS;
     
while 1  
    % Sampling in subregion
    tempPool = SDMCS_RVGenerator(p1,p2,ProSys.Ndim,detaNs);
    SBM.SamplePool = [SBM.SamplePool;tempPool];
    SBM.SubPoolIndex{n} = [SBM.SubPoolIndex{n},SBM.NofSamples+1:SBM.NofSamples+detaNs];
    SBM.NofSamples = SBM.NofSamples+detaNs;
    SBM.NoS_ring(n) = SBM.NoS_ring(n)+detaNs;
    SBM.CPI = [SBM.CPI,SBM.SubPoolIndex{n}];
    str1 = ['Generate ',num2str(detaNs),' samples in subregion D',num2str(n)];
    str2 = ['Current NofSamples=',num2str(SBM.NofSamples)];
    disp(['SDMCS:',str1,' ',str2])
    if ALSMTimeHis.false_ture_statistics == 1
        temp_G_ture = ProSys.fname(tempPool,ProSys.fun_par);
        SBM.G_ture = [SBM.G_ture;temp_G_ture];
    end
 
    
% -------------------------------------------------------------------------
%                     Active learning surrogate model
% -------------------------------------------------------------------------
        SBM.CPI = 1:SBM.NofSamples;
        SPI2 = find(ALSMTimeHis.NofSignUnchange>30);  % 符号保持不变n次就不再进行主动学习      
        SBM.CPI(unique([SBM.SPI';SPI2]))=[];   
        if length(SBM.CPI)==1
            SBM.CPI = [SBM.CPI,1:100];
        end
       [SBM,DoE,SurrModelPar,ALSMPar,ALSMTimeHis] = ...
            AL_Metamodel(ProSys,SBM,ALSMPar,DoE,SurrModelPar,ALSMTimeHis);    
% -------------------------------------------------------------------------
%                                end
% -------------------------------------------------------------------------
    
% 6. Computation of Cov(Pf). 
    % Calculation of Pf
    for ni = 1:SBM.NoR
        SubPool_response = SBM.G_predict(SBM.SubPoolIndex{ni},1);
        nf_SubPool = length(find(SubPool_response<=0));
        Pi(ni) = nf_SubPool/ SBM.NoS_ring(ni);
        PFi(ni) = Pi(ni) * SBM.thtai(ni);
        VarPi(ni) =  1/SBM.NoS_ring(ni) * Pi(ni)*(1-Pi(ni));
    end
    Pf =  sum(PFi);
    % Coffcient of variacne and fialure probability   
    VarPFi = SBM.thtai.^2.*VarPi;
    Cov_pf = sum(sqrt(VarPFi))/Pf;
    
% 7. Enrich additional samples to the subregion with the largest estimator variance.     
    if Cov_pf > SBM.cov_pf_tol        
%         [~,n] = max(VarPFi);                % 按照对总变异系数的贡献大小进行选择
          [~,n] = max(VarPFi./SBM.NoS_ring);  % 按照对变异系数降低的效率进行选择

        p2 = SBM.P_Divde_CDF(n+1); p1 = SBM.P_Divde_CDF(n);
        detaNs = min(SBM.detaNs_max,SBM.NoS_ring(n));
        
% 8. End of SDMCS.        
    else
        break
    end
end