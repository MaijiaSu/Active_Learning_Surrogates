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
% ! ע��
% 1.����Subset Simulation (SS)��Ϊ����������ѡ�صķ������м���̶���������ѧϰ����
%   ��AK-SSֻ����һ������ѧϰ
% 2.�÷����������ڼ������������ʱ��������������ã������㷨����ԣ�ʵ��������
%   �ᱻ�������ι��ܺ�����
% 3.��ʱ��֧�ּ�����ֵģ�ⷨSS�Ĺ��Ʒ���
% 4.��ʱ��������SS����ΪSS�����������ɵ�i�ε�������ģ�����ɣ����ɺ���Ҫ��������ѧϰ��
%   ģ�͸��º������µ�
%%
% 0. Initializtion of Subset Simulation 
     SBM.iniNoS = 1e4;
     SBM.SubSize = 1e4;
     SBM.P0 = 0.1;          
     % ����������ϸ����ܶȺ���
     RVPdf = @(x) MyPoolJointPdf(ProSys.muX,ProSys.sigmaX,ProSys.Distri,x); 
     % ����ֲ������ܶȺ���
     VM = 1;         %  ����ֲ�����   
     PropPdf = @(x,y) mvnpdf(x,y,diag((VM*ProSys.sigmaX).^2));
     % ����ֲ�������
     PropRnd = @(x) mvnrnd(x,diag((VM*ProSys.sigmaX).^2));
     % ������δʵ�֣�
     %      nchain = 1;     %  ����������������
     %      burnin = 100;   %  ���������Գ�ʼ״̬��
     %      thin = 1;       %  ��Ƭ������������
         
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

%% ��ͼ
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