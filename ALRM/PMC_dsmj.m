function  [DoE,PMC_DoE] = PMC_dsmj4(ProSys,PMC_Par,Pool_AL,SurrModel)
                                                
% initial_doe_creat
% ����α����Ʒ������ɳ�ʼDoE.X
% input: �������ͳ�Ʋ���  niuX sigmaX 
%        LHS���ɳ�ʼ���� ������������Ʒ��� n0 ��LHS����Ŵ�ϵ��VM
%        ����ģ��ѵ������ MoedlTrain�� ����ģ��Ԥ�⺯�� ModelPredictor 
%        ÿһ�����ɱ�ѡ״̬����nc������ֲ���׼�����ϵ��VM2
%        ���ܺ����� fname;����Ӧ����fun_par 
% output α������������ɵ�DoE.X ����Ӧ���ܺ���ֵG PMC_DoE

%% ��ʼ�� 
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

t=0;                         % ��t��״̬
PMC_developing = 1:Chain_N0; % ��δ���������������

for i = 1:Chain_N0
    PMC_DoE{i} = i;
end
%%
% 1.������Pool_AL���й�һ��
    Ndim = size(Pool_AL,2);
    for i=1:Ndim
        Norm_coeff(1,i) = max(abs(Pool_AL(:,i)));  % Normalization_coefficient
        Pool_AL(:,i) = Pool_AL(:,i)/ Norm_coeff(1,i);
    end
    
% 2.LHS���ɳ�ʼDoE.X
%     DoE.X = RNgeneratorV2(muX,VM*sigmaX,'LHS',N0);
      DoE.X = RNgeneratorV2(muX,VM*sigmaX,'LHS_uniform',Chain_N0);
      DoE.Y = fname(DoE.X,fun_par);
      
% 3. ѭ������
 while Chain_N0 > 0
     
     t = t + 1;   %��������t״̬ 
     PMC_STOP_NO = []; 
     disp('---------------------------')
     for i = 1:length(PMC_developing)
         
         disp(['t=' num2str(t) '��������������' num2str(i) '/' num2str(length(PMC_developing))])
         % 3.1���ڵ�ǰdoe�������ģ��
         
         SurrogateModel = MoedlTrain(DoE.X,DoE.Y);

         % 3.2 ���ɵ�i��α�������ı�ѡ״̬     
         
         N_iPMC_t = PMC_DoE{PMC_developing(i)}(end);   % ��i������ǰ״̬��DoE.X���
         Dpool2DoE.X = pdist2(Pool_AL,DoE.X(N_iPMC_t,:)./Norm_coeff,'euclidean');  % ��ѡ���е㵽��ǰ������״̬��ľ���
         doe_candidate = [];vm = 0; 
         
         while size(doe_candidate,1) == 0  %ֱ����ѡ״̬���ϲ�Ϊ��
             vm = vm+VM2;
             doe_candidate_loc = find(Dpool2DoE.X<vm*2*sqrt(Ndim));
             doe_candidate = Pool_AL(doe_candidate_loc,:);
         end
     
         % 3.3 ����ѡ���һ����ԭ
         
         doe_candidate = doe_candidate.*...
             (ones(size(doe_candidate,1),1)*Norm_coeff);
         
         % 3.4 ����ѧϰ����ѡȡ��ѵ�ʵ���
         
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
         
         % 3.5 �ж��Ƿ���ܸ�״̬,���Ƿ񽫸õ���뵽��������
                  
         if abs(DoE.Y(end)) <= abs(DoE.Y(N_iPMC_t))  
             PMC_DoE{PMC_developing(i)} = [PMC_DoE{PMC_developing(i)},length(DoE.Y)];   
         end
         
         % 3.6 �жϵ�i���������Ƿ�������������
         
         s1 = DoE.Y(PMC_DoE{PMC_developing(i)}(1));
         s2 = DoE.Y(PMC_DoE{PMC_developing(i)}(end));
         
         if s1*s2 <= 0  % ��¼ֹͣ��������������ţ�״̬t����������ɾ��
             PMC_STOP_NO = [PMC_STOP_NO,i];            
         end
         
     end 
     
     PMC_developing(PMC_STOP_NO) = [];  % ɾ�������������������
     
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