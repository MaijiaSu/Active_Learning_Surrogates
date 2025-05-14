clc,clear
% ProSys.muX = [0,0];
% ProSys.sigmaX = [1,1];
% ProSys.Distri = [1,1];
% ProSys.fname=@FBSS;
% ProSys.fun_par.type=1;ProSys.fun_par.a = 3;ProSys.fun_par.b = 7;
% % ProSys.fun_par.type=1;ProSys.fun_par.a = 5.5;ProSys.fun_par.b = 11;
% ProSys.fun_par.f0 = 0;
% ProSys.Ndim = 2;

TestExample = 'eg1';
LimtStateFunction_select
       
%%
% ! 注意
% 1.采用Subset Simulation (SS)作为生成样本候选池的方法，中间过程都进行主动学习，而
%   非AK-SS只进行一次主动学习
% 2.该方法仅适用于简单算例，计算耗时的算例不建议采用（考虑算法简洁性，实际样本点
%   会被调用两次功能函数）
% 3.暂时不支持计算数值模拟法SS的估计方差
% 4.暂时放弃采用SS，因为SS的子样本是由第i次迭代代理模型生成；生成后需要进行主动学习，
%   模型更新后会产生新的
%%
% 0. Initializtion of Subset Simulation 
     SBM.iniNoS = 1e4;
     SBM.SubSize = 1e4;
     SBM.P0 = 0.1;          
     % 随机变量联合概率密度函数
     RVPdf = @(x) MyPoolJointPdf(ProSys.muX,ProSys.sigmaX,ProSys.Distri,x); 
     % 建议分布概率密度函数
     VM = 1;         %  建议分布参数   
     PropPdf = @(x,y) mvnpdf(x,y,diag((VM*ProSys.sigmaX).^2));
     % 建议分布生成器
     PropRnd = @(x) mvnrnd(x,diag((VM*ProSys.sigmaX).^2));
     % 其它（未实现）
     %      nchain = 1;     %  马氏链数，正整数
     %      burnin = 100;   %  马氏链忽略初始状态数
     %      thin = 1;       %  切片步长，正整数
         
% 1. Crude MCS
     SBM.NoSubset = 1;   % November of subset
     [SBM.SamplePool,~] = RNgeneratorV2...
                      (ProSys.muX,ProSys.sigmaX,ProSys.Distri,SBM.iniNoS);
     SBM.NofSamples = SBM.iniNoS;
     SBM.SubPoolIndex{1} = 1:SBM.NofSamples;
     SBM.NoS_SubSet{1} = SBM.iniNoS;
          
% 2. Automatic generation of sub-samples
      bk = inf;

while  1
    
    % Calculate the response, and determine the next threshold 
    % ---------------------------------------------------------------------
    G = ProSys.fname(SBM.SamplePool(SBM.SubPoolIndex{SBM.NoSubset},:),ProSys.fun_par);
    % ---------------------------------------------------------------------    
    [Gsort,temp_index] = sort(G,'descend');
    bk = Gsort(ceil(SBM.NoS_SubSet{SBM.NoSubset}*(1-SBM.P0)));
    
    if bk<=0
        break
    end
    
    % Geneate the subset samples obeying the new Sub-TargetPdf
    index = temp_index(ceil(SBM.NoS_SubSet{SBM.NoSubset}*(1-SBM.P0)));   
    Xstart = SBM.SamplePool(SBM.SubPoolIndex{SBM.NoSubset}(index),:);
    SBM.bk(SBM.NoSubset) = bk;
      
    SBM.NoSubset = SBM.NoSubset+1;    
    Pfk = SBM.P0;
    TargetPdf = @(x) IndicatorFun(x,ProSys.fname,ProSys.fun_par,bk)*RVPdf(x)/Pfk;                        
    tempPool = [];
        
    tic
    tempPool =  MCMCgenerator(Xstart,SBM.SubSize,TargetPdf,PropPdf,PropRnd);
    toc
       
    % Save the data
    SBM.SamplePool = [SBM.SamplePool;tempPool];
    SBM.SubPoolIndex{SBM.NoSubset} = SBM.NofSamples+1:SBM.NofSamples+SBM.SubSize;
    SBM.NofSamples = SBM.NofSamples+SBM.SubSize;
    SBM.NoS_SubSet{SBM.NoSubset} = SBM.SubSize;
                  
end

% 3. calculate the Pf
P_end = length(find(G<0))/SBM.NoS_SubSet{end};
Pf = SBM.P0^(SBM.NoSubset-1) * P_end;

%% 绘图
colo={'bo'  'go' 'ro' 'co' 'mo'  'yo' 'ko' ...
    'b*'  'g*' 'r*' 'c*' 'm*'  'y*' 'k*'};
colo = repmat(colo,1,10);
for ii = 1:SBM. NoSubset
    plot(SBM.SamplePool(SBM.SubPoolIndex{ii},1),SBM.SamplePool(SBM.SubPoolIndex{ii},2),colo{ii});
    hold on
end
Pf

%%
function value = IndicatorFun(X,fname,fun_par,bk)
        Gk = fname(X,fun_par)-bk;
        value =1; 
       if Gk >= 0
           value =0;
       end
end