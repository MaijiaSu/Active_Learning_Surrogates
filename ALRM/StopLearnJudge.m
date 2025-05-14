function [LFSC,ALSMTimeHis] = StopLearnJudge(ProSys,ALSMTimeHis,times,ALSMPar,G_predict,Gmse,SurrModel,SBM,CandidatePool) 
% 计算终止学习函数指标
% 根据选定的终止条件，判断是否终止学习

LFSC = 0;
%%  Storage of the predicted results in successive n_save iterations
    n_save = 500; % Candidate pool keep uncharged
    if length(SBM.G_predict) == size(ALSMTimeHis.Gpredict_his,1) 
        ALSMTimeHis.Gpredict_his = [ALSMTimeHis.Gpredict_his,SBM.G_predict];
        if size(ALSMTimeHis.Gpredict_his,2)>n_save
            ALSMTimeHis.Gpredict_his(:,1) = [];
        end
    else % Candidate pool is updateed              
        nn1 = length(SBM.G_predict); nn2 = size(ALSMTimeHis.Gpredict_his,1);
        detaN = nn1-nn2;
        n_temp = min(size(ALSMTimeHis.old_SurrogateModel,2),n_save);        
        temp_Gpredict_his = zeros(detaN,n_temp);
        for n = 1:n_temp
            temp_Gpredict_his(:,n) = SurrModel.ModelPredictor...
                (ALSMTimeHis.old_SurrogateModel{end-n_temp+n},SBM.SamplePool(nn2+1:nn1,:));
        end
        ALSMTimeHis.Gpredict_his = [ALSMTimeHis.Gpredict_his;temp_Gpredict_his];    
    end
    
%% Record of the times of signs' changed 
if ~strcmp(SBM.method,'MCS_UQ')
    if size(ALSMTimeHis.Gpredict_his,2) == 1
        ALSMTimeHis.NofSignUnchange = zeros(size(SBM.G_predict,1),1);
    elseif size(ALSMTimeHis.Gpredict_his,2)>= 2
        nn1 = length(SBM.G_predict); nn2 = size(ALSMTimeHis.NofSignUnchange,1);
        if nn1 > nn2    % 样本池更新
            ALSMTimeHis.NofSignUnchange = [ALSMTimeHis.NofSignUnchange;zeros(nn1-nn2,1)];
        end
         temp_state = sign(ALSMTimeHis.Gpredict_his(:,end))-sign(ALSMTimeHis.Gpredict_his(:,end-1));
         temp_index1 = find(temp_state==0);
         temp_index2 = find(temp_state==1);
         ALSMTimeHis.NofSignUnchange(temp_index1) = ALSMTimeHis.NofSignUnchange(temp_index1)+1; % 保持不变次数加1
         ALSMTimeHis.NofSignUnchange(temp_index2) = 0;  % 保持不变次数清零
    end
end 

% if strcmp(SBM.method,'MCS_UQ')
%     % ------------------------------------------------------------------
%     
%     ymin = ProSys.ymin;
%     ymax = ProSys.ymax;
%     BW = (ymax-ymin)/SBM.NofInterval;
%     yy = ymin:BW:ymax;
%     if times == 1
%         ALSMTimeHis.errorCDF(times) = NaN;
%     else
%         CDFi = ALSMTimeHis.CDF(end-1,:);
%         CDFj = ALSMTimeHis.CDF(end,:);
%         CCDFi = 1-CDFi;CCDFi(CCDFi<0)=0;
%         deta = abs(CDFi-CDFj)./min(CDFi,CCDFi);
%         ALSMTimeHis.errorCDF(times) = trapz(yy,deta)/(ymax-ymin);
%     end
%     % ------------------------------------------------------------------
% end
%%  Calculating index of SC_UQLib  
if strcmp(ALSMPar.Stopcon.type,'SC_UQLab')...
       ||strcmp(ALSMPar.Stopcon.type,'CompositeSC')...
       ||strcmp(ALSMPar.Stopcon.type,'Modifie_SC_UQLab')

    % 储存连续3次的预测结果
    SC_iter1 = 0;
    SC_iter2 = 0;
    if size(ALSMTimeHis.Gpredict_his,2)>= 3
        % index 1
        nf1 = length(find(ALSMTimeHis.Gpredict_his(:,end-2)<0));
        nf2 = length(find(ALSMTimeHis.Gpredict_his(:,end-1)<0));
        nf3 = length(find(ALSMTimeHis.Gpredict_his(:,end)<0));
        ALSMTimeHis.SC1(times-1) = abs(nf2-nf1)/(nf1+1e-20);
        ALSMTimeHis.SC1(times) = abs(nf3-nf2)/(nf2+1e-20);
        
        % index 2
        temp1 = sign(ALSMTimeHis.Gpredict_his(:,end-1))-sign(ALSMTimeHis.Gpredict_his(:,end-2));
        temp2 = sign(ALSMTimeHis.Gpredict_his(:,end))-sign(ALSMTimeHis.Gpredict_his(:,end-1));
        NoSignChange1 = length(find(temp1~=0));
        NoSignChange2 = length(find(temp2~=0));
        ALSMTimeHis.SC2(times-1) = NoSignChange1/length(G_predict);
        ALSMTimeHis.SC2(times) = NoSignChange2/length(G_predict);
        
        % Judgement
        SC_iter1 = (ALSMTimeHis.SC1(times-1)<ALSMPar.Stopcon.SC1)&&...
            (ALSMTimeHis.SC2(times-1)<ALSMPar.Stopcon.SC2);
        SC_iter2 = (ALSMTimeHis.SC1(times)<ALSMPar.Stopcon.SC1)&&...
            (ALSMTimeHis.SC2(times)<ALSMPar.Stopcon.SC2);
    end
end

%% Judgement of stopping condition
if strcmp(ALSMPar.Stopcon.type,'minU')
    StopConditonStyle = 1;
elseif strcmp(ALSMPar.Stopcon.type,'maxEFF')
    StopConditonStyle = 2;
elseif strcmp(ALSMPar.Stopcon.type,'SC_UQLab')
    StopConditonStyle = 3;
elseif  strcmp(ALSMPar.Stopcon.type,'CompositeSC')
    StopConditonStyle = 4;
elseif strcmp(ALSMPar.Stopcon.type,'Modifie_SC_UQLib')
    StopConditonStyle = 5;
elseif strcmp(ALSMPar.Stopcon.type,'SC_ASVM-MCS')
    StopConditonStyle = 6;
elseif strcmp(ALSMPar.Stopcon.type,'SC_FPD')
    StopConditonStyle = 7;
elseif strcmp(ALSMPar.Stopcon.type,'SC_FPD_Stability')
    StopConditonStyle = 8;
elseif strcmp(ALSMPar.Stopcon.type,'SC_FD_composite')
     StopConditonStyle = 9;
end

switch StopConditonStyle
% #1: minU > 2
% #2: maxEFF < 0.001
% #3: SC_UQLib: SC1<0.001&SC2<0.001
% #4: minU > 2 || (SC1<0.001&SC2<0.001) !!{Compound condition of #1 and #3}
    case 1 % #1: minU > 2
        % Judgement
        if ALSMTimeHis.LFindex(times) > ALSMPar.Stopcon.minU
            LFSC = 1;
            disp('Stop learning：minU > 2')
        end
        
    case 2 % #2: maxEFF < 0.001 
        % Judgement
        if ALSMTimeHis.LFindex(times) < ALSMPar.Stopcon.maxEFF 
            LFSC = 1;
            disp('Stop learning：maxEFF < 0.001')
        end
        
    case 3 % #3: SC_UQLab: SC1<0.001&SC2<0.0             
        if SC_iter1&&SC_iter2
            LFSC = 1;
            disp('Stop learning：SC_UQLib: SC1<0.001&SC2<0.0001')
        end
                
    case 4  % minU > 2 || (SC1<0.001&SC2<0.001) 
        
        % Judgement1
        if (ALSMTimeHis.LFindex(times) > ALSMPar.Stopcon.minU)
            LFSC = 1;
            disp('Stop learning：minU > 2 (√)|| (SC1<0.001&SC2<0.001)')
        end     
        
        % Judgement2
        if SC_iter1&&SC_iter2
            LFSC = 1;
            disp('stop learning：minU > 2|| (SC1<0.001&SC2<0.001) (√)')
        end     
        
    case 5
        SC_iter1 = 0;
        SC_iter2 = 0;
        if size(ALSMTimeHis.Gpredict_his,2) >= 3
            % index 1
            nf1 = length(find(ALSMTimeHis.Gpredict_his(:,end-2)<0));
            nf2 = length(find(ALSMTimeHis.Gpredict_his(:,end-1)<0));
            nf3 = length(find(ALSMTimeHis.Gpredict_his(:,end)<0));
            ALSMTimeHis.SC1(times-1) = abs(nf2-nf1)/(nf1+1e-20);
            ALSMTimeHis.SC1(times) = abs(nf3-nf2)/(nf2+1e-20);          
            % index 2
            temp1 = sign(ALSMTimeHis.Gpredict_his(:,end-1))-sign(ALSMTimeHis.Gpredict_his(:,end-2));
            temp2 = sign(ALSMTimeHis.Gpredict_his(:,end))-sign(ALSMTimeHis.Gpredict_his(:,end-1));
            NoSignChange1 = length(find(temp1~=0));
            NoSignChange2 = length(find(temp2~=0));
            ALSMTimeHis.SC2(times-1) = NoSignChange1/(length(find(G_predict<0))+1e-20);
            ALSMTimeHis.SC2(times) = NoSignChange2/(length(find(G_predict<0))+1e-20);          
            % Judgement
            SC_iter1 = (ALSMTimeHis.SC1(times-1)<ALSMPar.Stopcon.SC1)&&...
                (ALSMTimeHis.SC2(times-1)<ALSMPar.Stopcon.SC2);
            SC_iter2 = (ALSMTimeHis.SC1(times)<ALSMPar.Stopcon.SC1)&&...
                (ALSMTimeHis.SC2(times)<ALSMPar.Stopcon.SC2);
        end
        
        if SC_iter1&&SC_iter2
            LFSC = 1;
            disp('Stop learning：SC_UQLib: SC1<0.001&SC2<0.0001')
        end
        
    case 6      
        MP=find(abs(G_predict)<1);   %magin point set,MP是相对在S中的位置
        % 终止评判标准1：超平面间隙中点的比例
        deta = length(MP)/length(G_predict);
        ALSMTimeHis.ASVM_deta(times) = deta;       
        %最小二乘收敛曲线
        k = length(ALSMTimeHis.ASVM_deta);
        p = polyfit(1:k,log(ALSMTimeHis.ASVM_deta),1);
        deta_predit = exp(polyval(p,k));
        e1 = ALSMPar.Stopcon.SC1; e2 = ALSMPar.Stopcon.SC2;      
        % 终止评判标准1：
        pd1 =max(deta,deta_predit)<e1;
        % 终止评判标准2：即-B*A*exp(B*K)
        pd2 = -p(2)*deta_predit<e2 && -p(2)*deta_predit>0; 
        if pd1&&pd2
            LFSC=1;  %如此终止学习
            disp('Stop learning：Judged by ASVM-MCS')
        end
        
    case 7
         % Calculation of three kinds of models
         ymin = ProSys.ymin;
         ymax = ProSys.ymax;
         BW = (ymax-ymin)/SBM.NofInterval;
         yy = ymin:BW:ymax;
         w_y = ALSMTimeHis.w_y(end,:);
         W_y = trapz(yy,w_y)/(ymax-ymin);    
         ALSMTimeHis.W_y(times) = W_y;
         % Judgement
        if W_y<ALSMPar.Stopcon.etol
            LFSC=1;  % stop learning
            disp('Stop learning：W_y<=e_tol')
        end

    case 8
        ymin = ProSys.ymin;
        ymax = ProSys.ymax;
        BW = (ymax-ymin)/SBM.NofInterval;
        yy = ymin:BW:ymax;
        w_y = ALSMTimeHis.w_y(end,:);
        W_y = trapz(yy,w_y)/(ymax-ymin);
        ALSMTimeHis.W_y(times) = W_y;
        if times >= 3
            Wy_final3 = ALSMTimeHis.W_y(end-2:end);
            error1 = abs(Wy_final3(2)-Wy_final3(1))/Wy_final3(1);
            error2 = abs(Wy_final3(3)-Wy_final3(2))/Wy_final3(2);     
            error_tol = ALSMPar.Stopcon.etol;
            if error1<error_tol&&error2<error_tol
                LFSC=1;  
                disp('Stop learning: reach to stability of error measure')
            end
        end         

    case 9
        % Calculation of three kinds of models
        ymin = ProSys.ymin;
        ymax = ProSys.ymax;
        BW = (ymax-ymin)/SBM.NofInterval;
        yy = ymin:BW:ymax;
        w_y = ALSMTimeHis.w_y(end,:);
        W_y = trapz(yy,w_y)/(ymax-ymin);
        ALSMTimeHis.W_y(times) = W_y;
        if times >= 3
            Wy_final3 = ALSMTimeHis.W_y(end-2:end);
            error1 = abs(Wy_final3(2)-Wy_final3(1))/Wy_final3(1);
            error2 = abs(Wy_final3(3)-Wy_final3(2))/Wy_final3(2);
            error_tol = ALSMPar.Stopcon.etol;
            if error1<1e-4&&error2<1e-4&&W_y<0.001
                LFSC=1;
                disp('Stop learning: reach to stability of error measure')
            end
        end

end
end