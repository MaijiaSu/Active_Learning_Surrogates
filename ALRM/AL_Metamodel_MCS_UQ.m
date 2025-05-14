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

% 3. Output of CDF and CCDF       
         ymin = ProSys.ymin;
         ymax = ProSys.ymax;
         BW = (ymax-ymin)/SBM.NofInterval;
         yy = ymin:BW:ymax;
         Edges = [-inf,ymin:BW:ymax,+inf];         
         % CDF
         CDF = histcounts(SBM.G_predict,Edges,'Normalization','CDF');CDF(end)=[];
         CDF_ture = histcounts(SBM.G_ture,Edges,'Normalization','CDF');CDF_ture(end)=[];
         CDF_error = (CDF- CDF_ture)./CDF_ture;
         CDF_CoV = sqrt((1-CDF)./CDF/SBM.NofSamples);     
         % CCDF
         CCDF = 1-CDF;
         CCDF_ture = 1-CDF_ture;
         CCDF_error = (CCDF-CCDF_ture)./CCDF_ture;
         CCDF_CoV = sqrt((1-CCDF)./CCDF/SBM.NofSamples);    
         
% 4. Estimation of statistical moment
        Moment(1) = mean(SBM.G_predict);
        Moment(2) = std(SBM.G_predict);
        Moment(3) = skewness(SBM.G_predict);
        Moment(4) = kurtosis(SBM.G_predict);
        
        Moment_ture(1) = mean(SBM.G_ture);
        Moment_ture(2) = std(SBM.G_ture);
        Moment_ture(3) = skewness(SBM.G_ture);
        Moment_ture(4) = kurtosis(SBM.G_ture);
        
% 5. Error mearsure        
         error_wy = abs(CDF-CDF_ture)./min([CDF_ture;1-CDF_ture]);
         error_Wy = trapz(yy,error_wy)/(ymax-ymin);
% 6. Stroage of data
SBM.CDF = CDF;
SBM.CDF_ture = CDF_ture;
SBM.CCDF = CCDF;
SBM.CCDF_ture = CCDF_ture;
SBM.Moment = Moment;
SBM.Moment_ture = Moment_ture;

%          figure
%          plot(yy,error_wy)
%          figure
%          plot(yy,CDF)
%          hold on
%          plot(yy,CDF_ture)
% 5. Enrichment of samples population?       

