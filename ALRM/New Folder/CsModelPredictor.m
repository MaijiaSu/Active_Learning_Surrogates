function [G_predict,Gmse] = CsModelPredictor(SurrModelPar,PredictPool)
% This function is used to calculate the Predicted Response and Jackknife
% variance
% innput£º  SurrModelPar, which is structural variables, with the filed: 
%  SurrogateModel, DoEs, and ks (the number of subsets of Jackknife method)
%           PredictPool, a sample pool
% output £º G_predict, predicted response at given sample pool points
%           Gmse, Jackknife variance at given sample pool points 

%%
   %% Calculate the predicted response
    N = size(PredictPool,1);
    if strcmp(SurrModelPar.Type,'SVM')
        [~,temp_G,~] = predict(SurrModelPar.SurrogateModel,PredictPool);
        temp_G(:,1)=[]; 
        G_predict = temp_G;
    elseif strcmp(SurrModelPar.Type,'SVR')
        G_predict = predict(SurrModelPar.SurrogateModel,PredictPool);
    elseif strcmp(SurrModelPar.Type,'PCE')
        G_predict = uq_evalModel(SurrModelPar.SurrogateModel,PredictPool);
    end
    Gmse = 0;

    %% Calculate the Jackknife variance
    if nargout == 2  
       
       % 1.Bulid the cross-validation models
         CSModel = MyCrossval(SurrModelPar);                
       
       % 2.Call to the cross-validation models
       ks = SurrModelPar.ks;
       g = zeros(N,ks);
       if strcmp(SurrModelPar.Type,'SVR') 
           for i = 1:ks
               g(:,i) = predict(CSModel{i},PredictPool);
           end
       elseif strcmp(SurrModelPar.Type,'SVM') 
           for i = 1:ks
               [~,temp_G,~] = predict(CSModel{i},PredictPool);
               temp_G(:,1)=[];
               g(:,i) = temp_G;
           end
       elseif strcmp(SurrModelPar.Type,'PCE')
           for i = 1:ks
               g(:,i) = uq_evalModel(CSModel{i},PredictPool);
           end    
       end
       
       % 3.Calculate the Jacknife value
       Jackknife = zeros(N,ks);
       Jackknife = ks*G_predict*ones(1,ks)-(ks-1)*g;
       
       % 4. Estimate the Jackknife variance
       Jsum = sum(Jackknife,2)/ks;
       temp = (Jackknife - Jsum *ones(1,ks)).^2;
       Gmse = 1/ks/(ks-1)*sum(temp,2);
       
   end
end