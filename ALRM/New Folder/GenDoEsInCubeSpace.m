function [DoE,errorWy,error_G,std_G] = GenDoEsInCubeSpace(ProSys,SurrModelPar,NN)

%% Example Definaiton
% clc,clear,close all
% TestExample = 'BatchProcessEg5'
% BenchmarkProblemForFullDistribution
% From my data-format to Uqlab
OriginInput = MyProSys2Uqlab(ProSys);


%%
% use different Experimental design method in constructing Metamodels
% NN = [32,64,128,256];
% NN = [128];
SampleMethod = {'MC','Sobol','Halton','LHS','MD'};

% 0. Generate the MCS samples and determine the simulation region
MCSPool = uq_getSample(OriginInput,1e5,'MC');
Y_ture = ProSys.fname(MCSPool,ProSys.fun_par);
bound = [min(MCSPool)',max(MCSPool)'];
for ii = 1:ProSys.Ndim
    InputUniformOpts.Marginals(ii).Name = ['U',num2str(ii)];
    InputUniformOpts.Marginals(ii).Type = 'Uniform';
    InputUniformOpts.Marginals(ii).Parameters = [bound(ii,:)];
end

UnifromInput = uq_createInput(InputUniformOpts);
% MCSPool_U = uq_getSample(UnifromInput,1e5,'MC');

for kk = 1:numel(NN)
    NofDoE = NN(kk);   
% 1. generate the samples in the Hypercube/normal space
    for nn = 1:numel(SampleMethod)   
        disp([num2str(kk),'/',num2str(numel(NN)),'---',num2str(nn),'/',num2str(numel(SampleMethod))])
        if strcmp(SampleMethod{nn},'MC')||strcmp(SampleMethod{nn},'Sobol')||...
                strcmp(SampleMethod{nn},'Halton')||strcmp(SampleMethod{nn},'LHS')
            DoE{kk,nn}.U = uq_getSample(UnifromInput,NofDoE,SampleMethod{nn});
        end
        if strcmp(SampleMethod{nn},'MD')
            DoE{kk,nn}.U = MaximunDistanceSelection(MCSPool,NofDoE);     
        end
    
% 2. Isoprobabilistic transform the points to original space   
DoE{kk,nn}.X  = DoE{kk,nn}.U; % do not need to transform

% 3. Evaluation on the perfomance function
DoE{kk,nn}.Y = ProSys.fname(DoE{kk,nn}.X,ProSys.fun_par);

% 4. Construct the metamodel and compute the error measure
% metamodel#1
% SurrModelPar.Type = 'UQLab_Kriging';
% SurrModelPar.MetaOpts.Trend.Type = 'ordinary';     % 'ordinary'(default), 'linear', and 'quadratic'; 'polynomial'
% SurrModelPar.MetaOpts.Corr.Type =  'Ellipsoidal';   % 'separable', 'Ellipsoidal'(default)
% SurrModelPar.MetaOpts.Corr.Family = 'Matern-5_2';    % Linear','exponential','Gaussian', 'Matern-3_2', 'Matern-5_2'(default)
% SurrModelPar.MetaOpts.EstimMethod = 'CV';    % 'CV','ML'
% metamodel#2
%  SurrModelPar = []
% SurrModelPar.Type = 'PCE';
% SurrModelPar.MetaOpts.Degree = 3:100;  % 3:18, 100
% metamodel#3
% SurrModelPar.Type = 'PCK';

[FullDistr{kk,nn},errorWy(kk,nn),MetamodelPredict{kk,nn},yy{kk,nn},error_w_y{kk,nn},error_G(kk,nn),std_G(kk,nn)] = ...
    errorWyEstimateV2(SurrModelPar,DoE{kk,nn}.U,DoE{kk,nn}.Y,MCSPool,Y_ture,ProSys);

    end
end

%% Compare the result
% figure('OuterPosition',[300,100,800,500])
% tempstr1 = {'$N_\mathcal{M}=32$',...
%            '$N_\mathcal{M}=64$','$N_\mathcal{M}=128$','$N_\mathcal{M}=256$',...
%            '$N_\mathcal{M}=512$'};
% tempstr = tempstr1(1:numel(NN))
% Xl = categorical(tempstr)
% Xl = reordercats(Xl,tempstr);
% bar(Xl,errorWy)
% legend(SampleMethod)
% xlabel('Number of DoE')
% ylabel('error measure')
% title(TestExample)
% tightfig
% saveas(gcf,['PCK-',TestExample,'.png'])

end