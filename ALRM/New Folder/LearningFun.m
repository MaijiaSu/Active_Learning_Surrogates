function [BPloc,LF_index,ALSMTimeHis]= LearningFun(SurrModelPar,G_predict,Gmse,LF_type,LF_Par,CandidatePool,DoE,ProSys,SBM,ALSMTimeHis,times)
% Call to the learning function to identify the best points
% input: G_predict
%        G_mse
%        LF_type
%        CandidatePool
%        DoE
% output: BPloc      £ºLocation fo BestPoint
%         LF_index   £ºThe index of learnfun value of best point
%% 
% Learnig function U, refferred to [1]
if strcmp(LF_type,'U')   
    Gstd = sqrt(Gmse);
    U = abs(G_predict./Gstd);
    [LF_index,BPloc] = min(U);
    
% Learnig function D, refferred to [2]    
elseif strcmp(LF_type,'D')
    Gmse_d = pdist2(CandidatePool,DoE.X,'euclidean');
    Gmse_d = min(Gmse_d')';
    Gstd = sqrt(Gmse_d);   
    U = abs(G_predict./Gstd);
    U(find(abs(G_predict)>=1)) = inf;
    [LF_index,BPloc] = min(U);
    
% Learnig function EFF, refferred to [3]      
elseif strcmp(LF_type,'EFF')
    Gstd = sqrt(Gmse);
    deta = 2*Gstd;
    G1 = -G_predict./Gstd;
    G2 = (-deta-G_predict)./Gstd;
    G3 = (deta-G_predict)./Gstd; 
    p1 = normcdf(G1,0,1);
    p2 = normcdf(G2,0,1);
    p3 = normcdf(G3,0,1);
    p4 = normpdf(G1,0,1);
    p5 = normpdf(G2,0,1);
    p6 = normpdf(G3,0,1);   
    EFF =  G_predict.*(2*p1-p2-p3)...
        - Gstd.*(2*p4-p5-p6)...
        + deta.*(p3 - p2);
    [LF_index,BPloc] = max(EFF);

% Learnig function UD, refferred to [4]    
elseif strcmp(LF_type,'UD')
    % Gmse_d
    Gmse_d = pdist2(CandidatePool,DoE.X,'euclidean');
    Gmse_d = min(Gmse_d')';
    Gmse_d = sqrt(Gmse_d);       
    Gmse_d = Gmse_d/max(abs(Gmse_d));
    % Gmse_e
    Gmse = Gmse/max(abs(Gmse)); 
    % learning function
    UD = abs(G_predict./(0.5*Gmse_d+0.5*Gmse));
    [LF_index,BPloc] = min(UD);

% Weighting penalty learning function P, proposed
elseif strcmp(LF_type,'P')
    % variance 
    % distance-variance    
    Gstd_d = pdist2(CandidatePool,DoE.X,'euclidean');
    Gstd_d = min(Gstd_d')';
    Gstd_d = sqrt(Gstd_d);
    % predictive variance    
    Gstd = sqrt(Gmse);
    Gstd_d = Gstd_d/max(abs(Gstd_d));
    Gstd = Gstd/max(abs(Gstd));
    % Learning function P
    c = LF_Par.c;
    alpha1 = LF_Par.alpha;
    alpha2 = 1-alpha1;
    P = abs(1./(alpha1*Gstd_d+alpha2*Gstd))+c*abs(G_predict)/max(abs(G_predict));
    [LF_index,BPloc] = min(P);

elseif strcmp(LF_type,'MoV')
    ymin = ProSys.ymin;
    ymax = ProSys.ymax;
    temp_index = 1:length(G_predict);
    Gstd = sqrt(Gmse);
%     indexL = temp_index((G_predict-2*Gstd>ymin)&(G_predict+2*Gstd<ymax));
     indexL = temp_index;
    [LF_index,index2] = max(Gstd(indexL));
    BPloc = indexL(index2);
    % Storage of data
    [w_y,CDF1,CDF2,CDF3] = wy_Calulate(ProSys,SBM,G_predict,Gstd);
    ALSMTimeHis.w_y(times,:) =  w_y;
    ALSMTimeHis.CDF3(times,:) = CDF3;
    ALSMTimeHis.CDF1(times,:) = CDF1;
    ALSMTimeHis.CDF(times,:) = CDF2;

elseif strcmp(LF_type,'MoV-distance')
    MeanC = mean(CandidatePool);
    stdC = std(CandidatePool);
    CandidatePool_norm = (CandidatePool-MeanC)./stdC;
    DoEX_norm = (DoE.X-MeanC)./stdC;
    tempDist = pdist2(CandidatePool_norm,DoEX_norm,'euclidean');
    DistVar = min(tempDist')';
    [LF_index,BPloc] = max(DistVar);
    % Storage of data
    Gstd = sqrt(Gmse);
    [w_y,CDF1,CDF2,CDF3] = wy_Calulate(ProSys,SBM,G_predict,Gstd);
    ALSMTimeHis.w_y(times,:) =  w_y;
    ALSMTimeHis.CDF3(times,:) = CDF3;
    ALSMTimeHis.CDF1(times,:) = CDF1;
    ALSMTimeHis.CDF(times,:) = CDF2;
   
elseif strcmp(LF_type,'MoV-composite')
    % variance 
    % distance-variance    
    Gstd_d = pdist2(CandidatePool,DoE.X,'euclidean');
    Gstd_d = min(Gstd_d')';
    Gstd_d = sqrt(Gstd_d);
    % predictive variance    
    Gstd = sqrt(Gmse);
    Gstd_d = Gstd_d/max(abs(Gstd_d));
    Gstd = Gstd/max(abs(Gstd));
    % composite variance
    alpha1 = 0.5;
    alpha2 = 0.5;
    std_c = (alpha1*Gstd_d+alpha2*Gstd);
    [LF_index,BPloc] = max(std_c);
    % Storage of data
    [w_y,CDF1,CDF2,CDF3] = wy_Calulate(ProSys,SBM,G_predict,Gstd);
    ALSMTimeHis.w_y(times,:) =  w_y;
    ALSMTimeHis.CDF3(times,:) = CDF3;
    ALSMTimeHis.CDF1(times,:) = CDF1;
    ALSMTimeHis.CDF(times,:) = CDF2;

elseif strcmp(LF_type,'MoV-Gradient')
    % variance
    % distance-variance
    MeanC = mean(CandidatePool);
    stdC = std(CandidatePool);
    CandidatePool_norm = (CandidatePool-MeanC)./stdC;
    DoEX_norm = (DoE.X-MeanC)./stdC;
    d2DoEs = pdist2(CandidatePool_norm,DoEX_norm,'euclidean');
    [Dmin,ID] = min(d2DoEs');  % Dmin is the distance to the nearest points
    Dmin = Dmin';
%     d2DoEs = pdist2(CandidatePool,DoE.X,'euclidean');

    % Gradient Residual 
    Resi = zeros(size(CandidatePool,1),1);
    Ndim = size(DoE.X,2);
    NofDoE = size(DoE.X,1); 
    for ii = 1:NofDoE
        id2 = find(ID==ii);
%         %     theta = SurrModelPar.SurrogateModel.theta;
%         dd = DoE.X(ii,:)-DoE.X; % differences between given data points
%         r = corrgauss(theta, dd);
        deta_r = zeros(Ndim,1);
        for kk = 1:Ndim                   
                    detaX = ProSys.sigmaX(kk)/1000;
                    xi = DoE.X(ii,:); xj = xi; xj(kk) = xj(kk)+detaX;
                    tempX = [xi;xj];
                    tempG = SurrModelPar.ModelPredictor(SurrModelPar,tempX);
                    deta_r(kk) = (tempG(1)-tempG(2))/detaX;
%             for jj = 1:NofDoE
% %                 deta_r(kk) = deta_r(kk)+2*theta(kk)*r(jj).*dd(jj);
%             end
        end
        dd = CandidatePool(id2,:)-DoE.X(ii,:); % differences between given data points
        G_L = DoE.Y(ii)+dd*deta_r;
        resi = abs(G_predict(id2)-G_L);
        Resi(id2) = resi;       
    end
    
    % Learnig function
    alpha1 = 1;
    Dmax_DoE = max(max(pdist2(DoE.X,DoE.X,'euclidean')));
    alpha2 = Dmin/Dmax_DoE;
%     alpha1 = 0.5;alpha2=0.5;
    LF_score = (alpha1*Dmin/max(Dmin)+alpha2.*Resi/max(Resi));
    [LF_index,BPloc] = max(LF_score);
    % Storage of data
    Gstd = sqrt(Gmse);
    [w_y,CDF1,CDF2,CDF3] = wy_Calulate(ProSys,SBM,G_predict,Gstd);
    ALSMTimeHis.w_y(times,:) =  w_y;
    ALSMTimeHis.CDF3(times,:) = CDF3;
    ALSMTimeHis.CDF1(times,:) = CDF1;
    ALSMTimeHis.CDF(times,:) = CDF2;
%     MyCmap
%     figure
%     for ii = 1:12
%         id2 = find(ID==ii);
%         plot(CandidatePool(id2,1),CandidatePool(id2,2),'.','Color',mycmap2(ii,:));
%         hold on
%     end

elseif strcmp(LF_type,'TwoStepLF')||strcmp(LF_type,'TwoStepLF_modified')   
    % Step 1. indentify candidate threshold
    ymin = ProSys.ymin;
    ymax = ProSys.ymax;
    BW = (ymax-ymin)/SBM.NofInterval;
    yy = ymin:BW:ymax;
    Edges = [-inf,ymin:BW:ymax,+inf];
    % Step 1.1 Calculation of symmetric mearsure w_y
        % Method#1
        if strcmp(LF_type,'TwoStepLF')
            Gstd = sqrt(Gmse);
            [w_y,CDF1,CDF2,CDF3] = wy_Calulate(ProSys,SBM,G_predict,Gstd);
            % Storage of data
            ALSMTimeHis.w_y(times,:) =  w_y;
            ALSMTimeHis.CDF3(times,:) = CDF3;
            ALSMTimeHis.CDF1(times,:) = CDF1;
            ALSMTimeHis.CDF(times,:) = CDF2;
        % Method#2
        elseif strcmp(LF_type,'TwoStepLF_modified')
            error('Updating')
        end 
    % Step 1.2 Calculation of WL_y 
        % % Case 1: Dirac kernel, WL_y = w_y
        if strcmp(LF_Par.Kernel,'TwoStepLF_DiracKernel')            
            WL_y = w_y;
        end
        % % Case 2: Gauusain kernel,WL_y = int(w_y,GaussianKernel)
        if strcmp(LF_Par.Kernel,'TwoStepLF_GaussianKernel')
            % i. determination of kernel width sigma_y of y'
            for ii = 1:length(yy)
                [~,temp_index] = min(abs(G_predict-yy(ii)));
                sigma_y(ii) = Gstd(temp_index);
            end
            % ii. Calculation of WL_y
            for ii = 1:length(yy)
                if sigma_y(ii)<BW
                    WL_y(ii) = interp1(yy,w_y,yy(ii));
                else
                    Z = sqrt(2*pi)*sigma_y(ii)*(normcdf((ymax-yy(ii))/sigma_y(ii))-...
                        normcdf((ymin-yy(ii))/sigma_y(ii)));
                    G_kernel = exp(-(yy-yy(ii)).^2/2/sigma_y(ii)^2);
                    WL_y(ii) = 1/Z*trapz(ymin:BW:ymax,w_y.*G_kernel);
                end
            end
        end
        [~,index] = max(WL_y);
        y_opt = yy(index);     

    % Step 2, indentify the optimal traning sample
    if strcmp(LF_Par.LearnFun,'U')
        U = abs(G_predict-y_opt)./Gstd;
        [LF_index,BPloc] = min(U);
        aa11 = 1;
    elseif strcmp(LF_Par.LearnFun,'P')
        c = LF_Par.c;   %  weighting coefficient
        P =   1./(Gstd/max(Gstd))+...
                      c*abs(G_predict-y_opt)/max(abs(G_predict-y_opt));
%         P =   1./(Gstd)+...
%                       c*abs(G_predict);
        [LF_index,BPloc] = min(P);
         
    end
end
   if numel(LF_index) == 0
       pause
   end
    
%% 
% [1] B. Echard, N. Gayton, M. Lemaire, AK-MCS: An active learning reliability
% method combining Kriging and Monte Carlo Simulation, STRUCT SAF, 33 (2011) 145-154.
% [2] Q. Pan, D. Dias, An efficient reliability method combining adaptive Support
% Vector Machine and Monte Carlo Simulation, 24 STRUCT SAF, 67 (2017) 85-95.
% [3] B.J. Bichon, M.S. Eldred, L.P. Swiler, S. Mahadevan, J.M. McFarland,
% Efficient Global Reliability Analysis for Nonlinear 22 Implicit Performance
% Functions, AIAA J, 46 (2008) 2459-2468.
% [4] N. Xiao, M.J. Zuo, C. Zhou, A new adaptive sequential sampling method to 
% construct surrogate models for efficient 36 reliability analysis, 
% RELIAB ENG SYST SAFE, 169 (2018) 330-33