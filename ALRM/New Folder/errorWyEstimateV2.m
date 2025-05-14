function [FullDistr,errorWy,MetamodelPredict,yy,error_w_y,error_G,std_G] = errorWyEstimateV2(SurrModelPar,DoEX,DoEY,...
                                            CandidatePool,Y_ture,ProSys)

% input: SurrModelPar.Type, DoE, Candidate_Pool,Y_ture,ProSys
% output：CDF,CCDF,errorWy

ymin = ProSys.ymin;
ymax = ProSys.ymax;

DoE.X = DoEX;
DoE.Y = DoEY;
%% Surrogate Model
% Transfrom the required data format of Uqlab
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

% Metamodels !!
if strcmp(SurrModelPar.Type,'Kriging')
    % Kriging Par.
    SurrModelPar.lob = 1e-6*ones(1,ProSys.Ndim);
    SurrModelPar.upb =1e6*ones(1,ProSys.Ndim);       
    SurrModelPar.theta=1e0*ones(1,ProSys.Ndim);   
    if  ~isfield(SurrModelPar,'RegressionModels')
        SurrModelPar.RegressionModels = @regpoly1;
    end
    if  ~isfield(SurrModelPar,'CorrelationModels')
        SurrModelPar.CorrelationModels = @corrgauss;
    end

elseif strcmp(SurrModelPar.Type,'ooKriging')  
    if  ~isfield(SurrModelPar,'opts')
        SurrModelPar.opts = struct();
    end
    if  ~isfield(SurrModelPar.opts,'type')
        SurrModelPar.opts.type = 'BlindKriging'; 
    end
    if  ~isfield(SurrModelPar.opts,'corrFunc')
       SurrModelPar.opts.corrFunc = @corrgauss; 
    end
%     SurrModelPar.opts.type = 'BlindKriging'; 
%     SurrModelPar.opts.corrFunc = @corrgauss;   
%     SurrModelPar.opts.corrFunc = @corrmatern52;
elseif strcmp(SurrModelPar.Type,'SVR')
    % SVR Par.
    SurrModelPar.KernelScale = 'auto';
    SurrModelPar.BoxConstraint = 1000;
    SurrModelPar.Epsilon = 0.001 ;
    SurrModelPar.KernelFunction =  'gaussian';       
    % KernelFunction = 'polynomial';PolynomialOrder = 2;
    SurrModelPar.ks = 5;           
elseif strcmp(SurrModelPar.Type,'SVM')
     % SVM Par.
     SurrModelPar.KernelScale = 'auto';
     SurrModelPar.KernelFunction = 'rbf';
     SurrModelPar.BoxConstraint = inf;
     SurrModelPar.ks = 5;           % 交叉验证分组数
elseif strcmp(SurrModelPar.Type,'PCE_Jackknife')    
    for ii = 1:  ProSys.Ndim
        InputOpts.Marginals(ii).Name = ['x',num2str(ii)]; % resistance variable
        if  ProSys.Distri(ii)== 1
            InputOpts.Marginals(ii).Type = 'Gaussian';
        elseif ProSys.Distri(ii)== 5
            InputOpts.Marginals(ii).Type = 'Uniform';
        else
            error('Need to update the code')
        end
        InputOpts.Marginals(ii).Moments = [ProSys.muX(ii),ProSys.sigmaX(ii)];
    end
    myInput = uq_createInput(InputOpts);   
    SurrModelPar.MetaOpts.Type = 'Metamodel';
    SurrModelPar.MetaOpts.MetaType = 'PCE';
    SurrModelPar.MetaOpts.TruncOptions.qNorm = 0.80;
    if  ~isfield(SurrModelPar.MetaOpts,'Method')
       SurrModelPar.MetaOpts.Method = 'LARS';
    end
    if  ~isfield(SurrModelPar.MetaOpts,'Degree')
         SurrModelPar.MetaOpts.Degree = 3:18;
    end 
    SurrModelPar.ks = 5;
    
elseif strcmp(SurrModelPar.Type,'PCE_UQLAB')||strcmp(SurrModelPar.Type,'PCE') 
    SurrModelPar.MetaOpts.Type = 'Metamodel';
    SurrModelPar.MetaOpts.MetaType = 'PCE';
    SurrModelPar.MetaOpts.TruncOptions.qNorm = 0.80;
    SurrModelPar.MetaOpts.Bootstrap.Replications = 100;  
    if  ~isfield(SurrModelPar.MetaOpts,'Method')
        SurrModelPar.MetaOpts.Method = 'LARS';
    end
    if  ~isfield(SurrModelPar.MetaOpts,'Degree')
        SurrModelPar.MetaOpts.Degree = 3:18;
    end
    
elseif strcmp(SurrModelPar.Type,'PCKriging')
    SurrModelPar.MetaOpts.Type = 'Metamodel';
    SurrModelPar.MetaOpts.MetaType = 'PCK';
    SurrModelPar.MetaOpts.Kriging.Corr.Family = 'Gaussian';
    SurrModelPar.MetaOpts.Mode = 'sequential'; % ='optimal','sequential'
    SurrModelPar.MetaOpts.PCE.Method = 'LARS';
     if  ~isfield(SurrModelPar.MetaOpts,'Degree')
        SurrModelPar.MetaOpts.PCE.Degree = 1:3;
     end

elseif  strcmp(SurrModelPar.Type,'UQLab_Kriging')
    SurrModelPar.MetaOpts.Type = 'Metamodel';
    SurrModelPar.MetaOpts.MetaType = 'Kriging';  
    % Options
    if ~isfield(SurrModelPar.MetaOpts,'Trend')
        SurrModelPar.MetaOpts.Trend.Type = 'ordinary';    % 'ordinary'(default), 'linear', and 'quadratic'; 'polynomial'(MetaOpts.Trend.Degree = q)
    end
    % Correlation function
    if ~isfield(SurrModelPar.MetaOpts,'Corr')
        SurrModelPar.MetaOpts.Corr = struct();
    end
    if ~isfield(SurrModelPar.MetaOpts.Corr,'Type')
        SurrModelPar.MetaOpts.Corr.Type = 'ellipsoidal';  % 'separable', 'Ellipsoidal'(default)
    end
    if ~isfield(SurrModelPar.MetaOpts.Corr,'Family')
        SurrModelPar.MetaOpts.Corr.Family = 'matern-5_2';  % Linear','exponential','Gaussian', 'Matern-3_2', 'Matern-5_2'(default)
    end
    if ~isfield(SurrModelPar.MetaOpts.Corr,'Isotropic')
        SurrModelPar.MetaOpts.Corr.Isotropic = false;
    end
    % Estimation method options
    SurrModelPar.MetaOpts.EstimMethod = 'CV';   % 'CV'(default, Metaopts.CV.LeaveKOut = 1), 'ML'
    SurrModelPar.MetaOpts.Optim.Method = 'HGA';
end

% 3.2 Surrogate model calling function
%  Unify the calling method of surrogate model
%  SurrogateModel = MoedlTrain(DoE,G) 
%  [G_predict,Gmse] = ModelPredictor(SurrogateModel,PredictPool)
[SurrModelPar.MoedlTrain,SurrModelPar.ModelPredictor] = SurrogateModelSlect(SurrModelPar);

%% Metemodels Training and Prediction
 %  Model Train
 SurrModelPar.SurrogateModel = SurrModelPar.MoedlTrain(DoE.X,DoE.Y);

 % Prediction
 if ~exist('Y_ture')
     Y_ture = ProSys.fname(CandidatePool,ProSys.fun_par);
 end

%  nn = 10000;
 NofSamples = size(CandidatePool,1);
%  Y_predict = zeros(NofSamples,1);
%  Ymse = zeros(NofSamples,1);

 [Y_predict,Ymse] =  SurrModelPar.ModelPredictor(SurrModelPar,CandidatePool);
 Ystd = sqrt(Ymse);

MetamodelPredict.Y_predict = Y_predict;
MetamodelPredict.Ymse = Ymse;
MetamodelPredict.ModelPredictor = @(X) SurrModelPar.ModelPredictor(SurrModelPar,X);

%% Compute relative absolute error
temperrorG = abs(Y_predict-Y_ture);
error_G = mean(temperrorG);
std_G = std(temperrorG);

% error_G = (abs(Y_predict-Y_ture));
% figure
% plot3(CandidatePool(:,1),CandidatePool(:,2),error_G,'.')
% hold on
% plot3(DoE.X(:,1),DoE.X(:,2),DoE.Y,'ro')
% 
% figure
% plot3(CandidatePool(:,1),CandidatePool(:,2),Y_predict,'.')
% hold on
% plot3(DoE.X(:,1),DoE.X(:,2),DoE.Y,'ro')
% 
% figure
% plot3(CandidatePool(:,1),CandidatePool(:,2),Y_ture,'.')
% hold on
% plot3(DoE.X(:,1),DoE.X(:,2),DoE.Y,'ro')


%% CDF,CCDF,errorWy 
% ymin = min(Y_ture);
% ymax = max(Y_ture);
BW = (ymax-ymin)/100;
yy = ymin:BW:ymax;
G1 = Y_predict-2*Ystd;
G = Y_predict;
G3 = Y_predict+2*Ystd;
% 
CDF1 = histcounts(G1,[-inf,yy],'Normalization','cumcount')/NofSamples;
CDF = histcounts(G,[-inf,yy],'Normalization','cumcount')/NofSamples;
CDF_exact = histcounts(Y_ture,[-inf,yy],'Normalization','cumcount')/NofSamples;
CDF3 = histcounts(G3,[-inf,yy],'Normalization','cumcount')/NofSamples;
% 
CCDF1 = 1-CDF1;
CCDF = 1-CDF;
CCDF_exact = 1-CDF_exact;
CCDF3 = 1-CDF3;
error_w_y = abs(CDF_exact-CDF)./min([CDF_exact;1-CDF_exact]);
error_w_y(isnan(error_w_y))=0;error_w_y(isinf(error_w_y))=0;
errorWy = trapz(yy,error_w_y)/(ymax-ymin);

FullDistr.CDF1 = CDF1;
FullDistr.CDF = CDF;
FullDistr.CDF_exact = CDF_exact;
FullDistr.CDF3 = CDF3;

FullDistr.CCDF1 = CCDF1;
FullDistr.CCDF = CCDF;
FullDistr.CCDF_exact = CCDF_exact;
FullDistr.CCDF3 = CCDF3;

return