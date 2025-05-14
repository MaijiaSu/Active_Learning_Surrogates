function [MoedlTrain,ModelPredictor] = SurrogateModelSlect(SurrModelPar,ProSys)
% 统一代理模型训练函数及预测函数的调用格式
%  SurrogateModel = MoedlTrain(DoE,G) 
%  [G_predict,Gmse] = ModelPredictor(SurrogateModel,PredictPool)
 
%% Default parameters of metamodels
if strcmp(SurrModelPar.Type,'Kriging')
    % Kriging Par.
    SurrModelPar.lob = 1e-6*ones(1,ProSys.Ndim);
    SurrModelPar.upb =1e6*ones(1,ProSys.Ndim);       
    SurrModelPar.theta=1e0*ones(1,ProSys.Ndim);   
    if  ~isfield(SurrModelPar,'RegressionModels')
        SurrModelPar.RegressionModels = @regpoly0;
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
    if  ~isfield(SurrModelPar.MetaOpts.PCE,'Degree')
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
        SurrModelPar.MetaOpts.Corr.Family = 'matern-5_2'  % Linear','exponential','Gaussian', 'Matern-3_2', 'Matern-5_2'(default)
    end
    if ~isfield(SurrModelPar.MetaOpts.Corr,'Isotropic')
        SurrModelPar.MetaOpts.Corr.Isotropic = false;
    end
    % Estimation method options
    SurrModelPar.MetaOpts.EstimMethod = 'CV';   % 'CV'(default, Metaopts.CV.LeaveKOut = 1), 'ML'
    SurrModelPar.MetaOpts.Optim.Method = 'HGA';
end



%%  Kriging   
if strcmp(SurrModelPar.Type,'Kriging')   
    MoedlTrain = @(DoE,G) dacefit...
        (DoE, G, SurrModelPar.RegressionModels,SurrModelPar.CorrelationModels,...
        SurrModelPar.theta,SurrModelPar.lob,SurrModelPar.upb);
    
    ModelPredictor = @(SurrModelPar,PredictPool)...
                          predictor_modified(PredictPool,SurrModelPar.SurrogateModel);  
%         ModelPredictor = @(SurrModelPar,PredictPool)...
%                           predictor(PredictPool, SurrogateModel);  

end
%% ooKriging

if strcmp(SurrModelPar.Type,'ooKriging')   
    MoedlTrain = @(DoE,G) oodacefit(DoE,G,SurrModelPar.opts);
    
    ModelPredictor = @(SurrModelPar,PredictPool)...
                          oopredictor(PredictPool,SurrModelPar.SurrogateModel);  
                      
end

%% ooKriging

if strcmp(SurrModelPar.Type,'ooKriging')   
    MoedlTrain = @(DoE,G) oodacefit(DoE,G,SurrModelPar.opts);
    
    ModelPredictor = @(SurrModelPar,PredictPool)...
                          oopredictor(PredictPool,SurrModelPar.SurrogateModel);  
                      
end

%% UQLab_Kriging
if strcmp(SurrModelPar.Type,'UQLab_Kriging')
        MoedlTrain = @(DoE,G) UQLabMetaModelTrain(DoE,G,SurrModelPar.MetaOpts);
    ModelPredictor = @(SurrModelPar,PredictPool) ...
        uq_evalModel(SurrModelPar.SurrogateModel,PredictPool);
end

%% SVR
if strcmp(SurrModelPar.Type,'SVR')
    
    if strcmp(SurrModelPar.BoxConstraint,'default') 
        MoedlTrain = @(Doe,G) fitrsvm(Doe,G,'KernelFunction',SurrModelPar.KernelFunction,...
            'Standardize',true,'Epsilon',iqr(G)/1.349,'KernelScale',SurrModelPar.KernelScale,...
             'BoxConstraint',iqr(G)/1.349,'Solver','SMO');
    else
        MoedlTrain = @(Doe,G) fitrsvm(Doe,G,'KernelFunction',SurrModelPar.KernelFunction,...
            'Standardize',true,'Epsilon',SurrModelPar.Epsilon,'KernelScale',SurrModelPar.KernelScale,...
            'BoxConstraint',SurrModelPar.BoxConstraint,'Solver','SMO');
    end
    
    ModelPredictor = @(SurrModelPar,PredictPool)...
         CsModelPredictor(SurrModelPar,PredictPool);
                         
end
%% SVM
if strcmp(SurrModelPar.Type,'SVM')
    
    MoedlTrain = @(DoE,G) fitcsvm(DoE,sign(G),...
        'KernelFunction',SurrModelPar.KernelFunction,'KernelScale',SurrModelPar.KernelScale,...
        'BoxConstraint',SurrModelPar.BoxConstraint,'Solver','SMO');
                                
    ModelPredictor = @(SurrModelPar,PredictPool) ...
         CsModelPredictor(SurrModelPar,PredictPool);
    
end
%% PCE — Jackknife Variance
 if strcmp(SurrModelPar.Type,'PCE_Jackknife')
    

    MoedlTrain = @(DoE,G) UQLabMetaModelTrain(DoE,G,SurrModelPar.MetaOpts)       
    ModelPredictor = @(SurrModelPar,PredictPool) ...
        CsModelPredictor(SurrModelPar,PredictPool);
                          
 end
%% PCE — Bootstrap
if strcmp(SurrModelPar.Type,'PCE')||strcmp(SurrModelPar.Type,'PCE_UQLAB')||strcmp(SurrModelPar.Type,'PCKriging')
    
    MoedlTrain = @(DoE,G) UQLabMetaModelTrain(DoE,G,SurrModelPar.MetaOpts)
    ModelPredictor = @(SurrModelPar,PredictPool) ...
        uq_evalModel(SurrModelPar.SurrogateModel,PredictPool);
    
end

return