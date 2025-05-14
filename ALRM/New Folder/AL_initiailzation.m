% Initialization of ALSModel

%%  
SBM.G_ture = [];
if ~isfield(ALSMTimeHis,'Gpredict_his')
    ALSMTimeHis.Gpredict_his = [];
end
if ~isfield(ALSMTimeHis,'StopLearningTimes')
    ALSMTimeHis.StopLearningTimes = [];
end

%% load DoEs from file
% iniState.iniDoE = DoE;
% iniState.iniSBM = SBM;
if ~strcmp(ALSMTimeHis.File_iniState,'lack')
    load(ALSMTimeHis.File_iniState);    
    if isfield(iniState,'iniDoE')
        DoE = StrucPass(DoE,iniState.iniDoE);
        disp('!!! Sucessfully load the initial DoE')
    end
    if isfield(iniState,'iniSBM')  
        SBM = StrucPass(SBM,iniState.iniSBM);
        SBM.CPI = 1:size(SBM.SamplePool,1);
        disp('!!! Sucessfully load the initial SBM')
    end
end

%%  Generate the initial DoEs (DoE.X��DoE.Y)
   if size(DoE.X,1) == 0
       if strcmp(ALSMPar.IniDoE.GenType,'iniRandom')
           % '=1', select from candidate Pool randomly
           CandidatePool = SBM.SamplePool(SBM.CPI,:);
           temp_index = randperm(size(CandidatePool,1),ALSMPar.IniDoE.N0);
           DoE.X = CandidatePool(temp_index,:);
           DoE.Y = ProSys.fname(DoE.X,ProSys.fun_par);
           SPI = temp_index;
       elseif strcmp(ALSMPar.IniDoE.GenType,'iniLHS')
           % '=2', Latin hypercubic sampling ��LHS��
           if sum(ProSys.Distri==5)>0
               DoE.X = RNgeneratorV2(ProSys.muX,ProSys.sigmaX,...
                   'LHS_uniform',ALSMPar.IniDoE.N0);
          elseif strcmp(ProSys.Distri,'mvnrnd')
              DoE.X = RNgeneratorV2(ProSys.muX,ProSys.sigmaX,...
                  'LHS',ALSMPar.IniDoE.N0);
           else 
               DoE.X = RNgeneratorV2(ProSys.muX,ALSMPar.IniDoE.VM *ProSys.sigmaX,...
                   'LHS_uniform',ALSMPar.IniDoE.N0);
           end
           DoE.Y = ProSys.fname(DoE.X,ProSys.fun_par);
       elseif strcmp(ALSMPar.IniDoE.GenType,'iniPMC')
           % '=3', Pseudo Markov Chains
           CandidatePool = SBM.SamplePool(SBM.CPI,:);
           [DoE,PMC_DoE] = PMC_dsmj(ProSys,ALSMPar.IniDoE,...
               CandidatePool,SurrModelPar);
           DoE.PMC_DoE = PMC_DoE;
           ALSMPar.IniDoE.N0 = length(DoE.Y);
           disp(['PMC_DOE=' num2str(ALSMPar.IniDoE.N0)]);
       elseif strcmp(ALSMPar.IniDoE.GenType,'iniMDS')
           CandidatePool = SBM.SamplePool(SBM.CPI,:);
           DoE.X = MaximunDistanceSelection(CandidatePool,ALSMPar.IniDoE.N0);
           DoE.Y =  ProSys.fname(DoE.X,ProSys.fun_par);
       end
   end
   
   if strcmp(SBM.method,'IS')
       DoE.X = [DoE.X;SBM.MPP_DoE.X];
       DoE.Y = [DoE.Y;SBM.MPP_DoE.Y];
   end
%% �޳��쳣�ĳ�ʼ���ݵ�  -
   % ��������Ԫ������Ϊ����Ŵ󣬿��ܻ���һЩ�쳣�����ݵ�
   while 1
       tempidnex = find(abs(DoE.Y)>1e20*abs(iqr(DoE.Y)));
       if isempty(tempidnex)
           break
       end
       DoE.Y(tempidnex,:) = [];
       DoE.X(tempidnex,:) = [];
       tempX = RNgeneratorV2(ProSys.muX,ALSMPar.IniDoE.VM *ProSys.sigmaX,'LHS_uniform',length(tempidnex));
       tempY = ProSys.fname(tempX ,ProSys.fun_par);
       DoE.X = [DoE.X;tempX];
       DoE.Y = [DoE.Y;tempY];
       add_DoE = 1;
   end
   
%% ѵ����ʼ����ģ�� 
% ѵ������ģ��
   SurrModelPar.SurrogateModel = SurrModelPar.MoedlTrain(DoE.X,DoE.Y);
   SurrModelPar.DoE = DoE;
   ALSMTimeHis.old_SurrogateModel{1} = SurrModelPar;
    
% ---------------------------����ģ��ѵ������-------------------------------
    % ���´���ģ��-ѵ��ģ��MoedlTrain����SVR/SVM��ģ�Ͳ����ǲ������Զ���̶��ķ�ʽ
    if strcmp(SurrModelPar.Type ,'SVR')||strcmp(SurrModelPar.Type ,'SVM')
       SurrModelPar.KernelScale =  SurrModelPar.SurrogateModel.KernelParameters.Scale;
       [SurrModelPar.MoedlTrain,SurrModelPar.ModelPredictor] = ...
                                         SurrogateModelSlect(SurrModelPar);
    end
% ------------------------------------------------------------------------- 
%% sbufunction
function struc1 = StrucPass(struc1,struc2)
% function: pass the same field of struc2 to struct1
    field_names =   fieldnames(struc2);
    for ii = 1:size(field_names,1)
        temp1 = strcmp(fieldnames(struc1),field_names(ii));
        if sum(temp1)>0
            jj = find(temp1==1);
            struc1 = setfield(struc1,field_names{ii},getfield(struc2,field_names{ii}));
        end
    end
end