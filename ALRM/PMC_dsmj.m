function  [DoE,PMC_DoE] = PMC_dsmj4(ProSys,PMC_Par,Pool_AL,SurrModel)
                                                
% initial_doe_creat
% 采用伪马尔科夫链生成初始DoE.X
% input: 随机变量统计参数  niuX sigmaX 
%        LHS生成初始个数 即多少条马尔科夫链 n0 及LHS方差放大系数VM
%        代理模型训练函数 MoedlTrain， 代理模型预测函数 ModelPredictor 
%        每一步生成备选状态个数nc；建议分布标准差调整系数VM2
%        功能函数名 fname;及响应参数fun_par 
% output 伪马尔卡夫链生成的DoE.X 及响应功能函数值G PMC_DoE

%% 初始化 
muX =  ProSys.muX;
sigmaX = ProSys.sigmaX;

Chain_N0 = PMC_Par.PMC_N0;                % Number of chains
VM =  PMC_Par.VM;                         % Generate the first state of each chain
VM2 = PMC_Par.PMC_VM2;                    % Variance magnification of local seaching region
PMC_MaxL = PMC_Par.PMC_MaxL;              % Maximum length of sigle chain
PMC_LF_type = PMC_Par.PMC_LF_type;        % Type of learing function to find out local best point 

fname = ProSys.fname;
fun_par = ProSys.fun_par;

MoedlTrain = SurrModel.MoedlTrain;
ModelPredictor = SurrModel.ModelPredictor;

t=0;                         % 第t个状态
PMC_developing = 1:Chain_N0; % 还未收敛的马氏链编号

for i = 1:Chain_N0
    PMC_DoE{i} = i;
end
%%
% 1.将数据Pool_AL进行归一化
    Ndim = size(Pool_AL,2);
    for i=1:Ndim
        Norm_coeff(1,i) = max(abs(Pool_AL(:,i)));  % Normalization_coefficient
        Pool_AL(:,i) = Pool_AL(:,i)/ Norm_coeff(1,i);
    end
    
% 2.LHS生成初始DoE.X
%     DoE.X = RNgeneratorV2(muX,VM*sigmaX,'LHS',N0);
      DoE.X = RNgeneratorV2(muX,VM*sigmaX,'LHS_uniform',Chain_N0);
      DoE.Y = fname(DoE.X,fun_par);
      
% 3. 循环生成
 while Chain_N0 > 0
     
     t = t + 1;   %马氏链第t状态 
     PMC_STOP_NO = []; 
     disp('---------------------------')
     for i = 1:length(PMC_developing)
         
         disp(['t=' num2str(t) '：链数迭代进度' num2str(i) '/' num2str(length(PMC_developing))])
         % 3.1基于当前doe构造代理模型
         
         SurrogateModel = MoedlTrain(DoE.X,DoE.Y);

         % 3.2 生成第i条伪马氏链的备选状态     
         
         N_iPMC_t = PMC_DoE{PMC_developing(i)}(end);   % 第i条链当前状态的DoE.X编号
         Dpool2DoE.X = pdist2(Pool_AL,DoE.X(N_iPMC_t,:)./Norm_coeff,'euclidean');  % 备选池中点到当前马氏链状态点的距离
         doe_candidate = [];vm = 0; 
         
         while size(doe_candidate,1) == 0  %直至备选状态集合不为空
             vm = vm+VM2;
             doe_candidate_loc = find(Dpool2DoE.X<vm*2*sqrt(Ndim));
             doe_candidate = Pool_AL(doe_candidate_loc,:);
         end
     
         % 3.3 将备选点归一化还原
         
         doe_candidate = doe_candidate.*...
             (ones(size(doe_candidate,1),1)*Norm_coeff);
         
         % 3.4 调用学习函数选取最佳的实验点
         
         if  strcmp(PMC_LF_type,'U')|| strcmp(PMC_LF_type,'EFF')
             [G_predict,mse] = ModelPredictor(SurrogateModel,doe_candidate);
             [BPloc,LF_index] = LearningFun(G_predict,Gmse,...
                 PMC_LF_type,doe_candidate,DoE);
         elseif  strcmp(PMC_LF_type,'D')
              G_predict = ModelPredictor(SurrogateModel,doe_candidate);
              Gmse = [];
              [BPloc,LF_index] = LearningFun(G_predict,Gmse,...
                  PMC_LF_type,doe_candidate,DoE);
         else
             error('PMC_LF_type is out of selected range')
         end
         
%          G_predict = ModelPredictor(SurrogateModel,doe_candidate);
%          Gmse = pdist2(doe_candidate,DoE.X,'euclidean');
%          Gmse = min(Gmse')'; Gmse = sqrt(Gmse);        
%          [~,index] = min(abs(G_predict)./Gmse);
         
         DoE.X = [DoE.X;doe_candidate(BPloc,:)];
         DoE.Y = [DoE.Y;fname(DoE.X(end,:),fun_par)];
         
         % 3.5 判断是否接受该状态,及是否将该点加入到马氏链中
                  
         if abs(DoE.Y(end)) <= abs(DoE.Y(N_iPMC_t))  
             PMC_DoE{PMC_developing(i)} = [PMC_DoE{PMC_developing(i)},length(DoE.Y)];   
         end
         
         % 3.6 判断第i条马氏链是否满足收敛条件
         
         s1 = DoE.Y(PMC_DoE{PMC_developing(i)}(1));
         s2 = DoE.Y(PMC_DoE{PMC_developing(i)}(end));
         
         if s1*s2 <= 0  % 记录停止迭代的马氏链编号，状态t结束将该链删除
             PMC_STOP_NO = [PMC_STOP_NO,i];            
         end
         
     end 
     
     PMC_developing(PMC_STOP_NO) = [];  % 删除已收敛的马氏链编号
     
     if t >= PMC_MaxL
         break
     end
         
 end
%%
% plot(Pool_AL(:,1),Pool_AL(:,2),'ko')
% hold on
% plot(doe_candidate(:,1),doe_candidate(:,2),'yo')
% hold on
% plot(0.8317,0.9763,'ro')