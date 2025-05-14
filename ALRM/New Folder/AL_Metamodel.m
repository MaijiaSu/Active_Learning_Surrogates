function    [SBM,DoE,SurrModelPar,ALSMPar,ALSMTimeHis] = ...
           AL_Metamodel(ProSys,SBM,ALSMPar,DoE,SurrModelPar,ALSMTimeHis)
%% 
% ALSModel is preformed to train a surrogate model using straegy of adaptively
% sqential sampling. The approximate sapce is perferred to be as minimal as the
% distribution region of candidate pool. 
% input: 
%     candidatePool and involved ALSM's Par.
%     ProSys Par. is involved in initial DoE and updation of DoE
%     DoE, SurrModel, ALSMTimeHis Par. wiil be updated in procedure
%  output:
%     G_predict 
%     DoE
%     SurrModel,ALSMTimeHis

%% 
%# AL_Step-1. Initialization of ALSModel
times = ALSMTimeHis.times;
LFSC = 0;
if ALSMTimeHis.times == 0
    AL_initiailzation
end
CandidatePool = SBM.SamplePool(SBM.CPI,:);

% determine to if it is necessary to calculate Pf_ture
% -------------------------------------------------------------------------     
 if ALSMTimeHis.false_ture_statistics == 1
     SBM.G_ture = SBM.G_ture(:);
     if size(SBM.G_ture,1)<SBM.NofSamples
         tempPool = SBM.SamplePool(length(SBM.G_ture)+1:SBM.NofSamples,:);
         temp_G_ture = ProSys.fname(tempPool,ProSys.fun_par);
         SBM.G_ture = [SBM.G_ture;temp_G_ture];
     end
 end
% -------------------------------------------------------------------------

% ------------------- start the Active learnig process --------------------
while LFSC == 0   

    times = times + 1;
    
%# AL_Step-2. Predict the respose of Candidate Pool   
    if  strcmp(ALSMPar.LF_type,'D')   % don't need to output the Gmse
        G_predict = SurrModelPar.ModelPredictor(SurrModelPar,CandidatePool);
        Gmse = [];
    else 
        [G_predict,Gmse] = SurrModelPar.ModelPredictor(SurrModelPar,CandidatePool);
    end

%# AL_Step-3. Calculate the Pf_predict or Pf_ture 
    SBM.G_predict(SBM.CPI,1) = G_predict;
    SBM.Gmse(SBM.CPI,1) = Gmse;
    if strcmp(SBM.method,'MCS_UQ')
        ymin = ProSys.ymin; ymax = ProSys.ymax;
        BW = (ymax-ymin)/SBM.NofInterval;
        Edges = [-inf,ymin:BW:ymax,+inf];   
        if ALSMTimeHis.false_ture_statistics == 1
            CDF_ture = histcounts(SBM.G_ture,Edges,'Normalization','cumcount')/SBM.NofSamples;CDF_ture(end)=[];
            ALSMTimeHis.CDF_ture = CDF_ture;
        end
    else
        SBM = AL_Estimation_Pf(SBM);
        % Storage of iterative data
        ALSMTimeHis.Pf_predict(times) = SBM.Pf_predict;
        ALSMTimeHis.Cov_predict(times) = SBM.Cov_predict;
        if ALSMTimeHis.false_ture_statistics == 1
            ALSMTimeHis.Pf_ture(times) = SBM.Pf_ture;
            ALSMTimeHis.Cov_ture(times) = SBM.Cov_ture;
        end
    end
    
%# AL_Step-4. Empolyment of learning function 
    [BPloc,LF_index,ALSMTimeHis] = LearningFun(SurrModelPar,G_predict,Gmse,...
                             ALSMPar.LF_type,ALSMPar.LF_Par,CandidatePool,DoE,ProSys,SBM,ALSMTimeHis,times);                                                    
    ALSMTimeHis.LFindex(times) = LF_index(1);    
    % Check the learning effect 
    ALSMTimeHis.NBP_dis(times) = min(pdist2(CandidatePool(BPloc,:),DoE.X,'euclidean'));
    if ALSMTimeHis.NBP_dis(times)<1e-20
       warning('The selected optimal point is too close to the existent DoEs')
    end
    
%# AL_Step-5. Judgement of teminataing the active learning process       
    [LFSC,ALSMTimeHis] = StopLearnJudge(ProSys,ALSMTimeHis,times,ALSMPar,G_predict,Gmse,SurrModelPar,SBM,CandidatePool);

%# AL_Step-6. Select the optimal points and update the metamodel
    if LFSC == 0      
        % Enrichment of DoE
        DoE.X = [DoE.X; CandidatePool(BPloc,:)];   
        DoE.Y = [DoE.Y; ProSys.fname(DoE.X(end,:),ProSys.fun_par)];     
        SBM.SPI = [SBM.SPI,BPloc]; % record the index of selected point       
              
        % Trian the surrogate model
        SurrModelPar.SurrogateModel = SurrModelPar.MoedlTrain(DoE.X,DoE.Y);
        SurrModelPar.DoE = DoE;
        % Storage of the meatamode
        N = size(ALSMTimeHis.old_SurrogateModel,2);
        ALSMTimeHis.old_SurrogateModel{N+1} = SurrModelPar;                                                      
    end

%# AL_Step-7. Terminate the active learning process
    % LFSC == 1, break the program loop

%# AL_Step-8.1 ProcessDisplay of active learning (Estiamtion of Pf )
% ------------------------------------------------------------------------- 
if strcmp(SBM.method,'MCS')||strcmp(SBM.method,'IS')||strcmp(SBM.method,'SDMCS')
% Display Text 
 str{1} = ['Iter#',num2str(times),':',...
     '(',num2str(toc(ALSMTimeHis.t_start),'%.1f'),'s)',...
     ' Pf=',num2str(SBM.Pf_predict),...
     ' Cov_pf=',num2str(SBM.Cov_predict)];
 str{2} = [];
 % Truth and False statistics
if ALSMTimeHis.false_ture_statistics == 1
    T1=sign(SBM.G_ture);T2=sign(SBM.G_predict);
    nf_ture = length(find(T1<=0));
    nf_predict = length(find(T2<=0));
    state = T1-T2;
    %  FPCP: FALSE positive classified points, '-' ¡ú '+' (È¡Î±,flase positive)
    %  FNCP: FALSE negative classified points, '+' ¡ú '-' (ÆúÕæ,flase negative)
    S1 = find(state==0);             % set of correct indentification 
    ALSMTimeHis.S2 = find(state>0);  % set of FNCP
    ALSMTimeHis.S3 = find(state<0);  % set of FPCP
    L1 = length(S1); L2 = length(ALSMTimeHis.S2); L3 = length(ALSMTimeHis.S3);
    ALSMTimeHis.MCP_state(times,:) = [L1,L2,L3]/ nf_ture;
    ALSMTimeHis.Pf_Error(times) = (SBM.Pf_predict-SBM.Pf_ture)/SBM.Pf_ture;        
    str{2}=[' Pf_ture=',num2str(ALSMTimeHis.Pf_ture(times)),...
        ' NoFPCP(-¡ú+)=',num2str(L3),' NoFNCP(+¡ú-)=',num2str(L2),...
        ' Pf_Error=',num2str(ALSMTimeHis.Pf_Error(times)*100),'%'];
end  
if ALSMTimeHis.ProcessDisplay == 1
    disp([str{1},str{2}])
end
end
% -------------------------------------------------------------------------

%# AL_Step-8.2 ProcessDisplay of active learning (Estimation of Full proability distribution )
% ------------------------------------------------------------------------- 
if strcmp(SBM.method,'MCS_UQ')
    str{1} = ['Iter#',num2str(times),':',...
        '(',num2str(toc(ALSMTimeHis.t_start),'%.1f'),'s)'];   
    if ALSMTimeHis.false_ture_statistics == 1
        CDF_ture = ALSMTimeHis.CDF_ture;
        error_wy = abs(CDF_ture-ALSMTimeHis.CDF(times,:))./min([CDF_ture;1-CDF_ture]);
        error_wy(isnan(error_wy))=0;error_wy(isinf(error_wy))=0;
        errorW_y = trapz(ymin:BW:ymax,error_wy)/(ymax-ymin);
        ALSMTimeHis.errorW_y(times) = errorW_y;
        str{2} = ['  W_y','=',num2str(ALSMTimeHis.W_y(times)),' ',...
                  'errorW_y = ',num2str(errorW_y)];
    else
        str{2} = ['  W_y','=',num2str(ALSMTimeHis.W_y(times))];
    end
    % Display Text
    disp([str{1},str{2}])
end

if ProSys.Ndim == 1
    % Display the distribution of DoE
    if times == 1
        p = plot(DoE.X(:,1),DoE.X(:,2),'k.');
        p.XDataSource = 'DoEX1t';
        p.YDataSource = 'DoEX2t';
        hold on
        vm = 5;
        for ii = 1:ProSys.Ndim
            if ProSys.Distri(ii) == 5
                bound(ii,:) = [ProSys.RVPar.P1(ii),ProSys.RVPar.P2(ii)];
            elseif ProSys.Distri(ii) == 6 || ProSys.Distri(ii) == 2
                bound(ii,:) = [0,ProSys.muX(ii)+vm*ProSys.sigmaX(ii)];
            else
                bound(ii,:) = [ProSys.muX(ii)-vm*ProSys.sigmaX(ii),ProSys.muX(ii)+vm*ProSys.sigmaX(ii)];
            end
        end
        gap = 100;
        MyFun = @(X) ProSys.fname(X,ProSys.fun_par);
        contour2D(bound,gap,MyFun);
        xlabel('$X_1$')
        ylabel('$X_2$')
    else
        DoEX1t = DoE.X(:,1);
        DoEX2t = DoE.X(:,2);
        refreshdata(p,'caller');
        drawnow
    end
end
% -------------------------------------------------------------------------

% Storage of data
ALSMTimeHis.times = times;
SurrModelPar.DoE = DoE;
% record the time of stop learning
ALSMTimeHis.StopLearningTimes = [ALSMTimeHis.StopLearningTimes,times];

%# AL_Step-9. Forced to stop?
    if size(DoE.X,1) >=ALSMTimeHis.MaxNoDoE        
        SBM.cov_pf_tol = inf;
        disp('!!! Maximum number of DoE is reached')
        break
    end
    if toc(ALSMTimeHis.t_start)>=ALSMTimeHis.Maximum_runing_time
        SBM.cov_pf_tol = inf;
        disp('!!! Maximum runing time is reached')
        break
    end
end
    
return   