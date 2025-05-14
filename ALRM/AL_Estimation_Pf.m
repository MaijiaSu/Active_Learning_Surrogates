function SBM = AL_Estimation_Pf(SBM)


GG = SBM.G_predict(:); 
if isfield(SBM,'G_ture')
   GG = [GG,SBM.G_ture];
end

for nn = 1:size(GG,2)
    
    G = GG(:,nn);
    
    if strcmp(SBM.method,'MCS')||strcmp(SBM.method,'MCS_UQ')
        
        nf = length(find(G<=0));
        Pf = nf/SBM.NofSamples;
        Cov_Pf = sqrt((1-Pf)/Pf/SBM.NofSamples);

    elseif strcmp(SBM.method,'IS')
        
        Fail_index = find(G<=0);
        Pf = 1/SBM.NofSamples *sum(SBM.Fpdf(Fail_index)./SBM.Ppdf(Fail_index));
        State = zeros(SBM.NofSamples,1); State(Fail_index,1) = 1;     
        temp1 = sum(State .*(SBM.Fpdf./SBM.Ppdf).^2);
        Var_Pf = 1/SBM.NofSamples*(1/SBM.NofSamples*temp1-Pf.^2);
        Cov_Pf = sqrt(Var_Pf)/Pf;
        
    elseif strcmp(SBM.method,'SDMCS')
        
        for ni = 1:size(SBM.SubPoolIndex,2)
            SubPool_response = G(SBM.SubPoolIndex{ni},1);
            nf_SubPool = length(find(SubPool_response<=0));
            Pi(ni) = nf_SubPool/ SBM.NoS_ring(ni);
            PFi(ni) = Pi(ni) * SBM.thtai(ni);
            VarPi(ni) =  1/SBM.NoS_ring(ni) * Pi(ni)*(1-Pi(ni));
        end
        % Current fialure probability
        Pf = sum(PFi); 
        VarPFi = SBM.thtai(1:ni).^2.*VarPi;
        Cov_Pf = sum(sqrt(VarPFi))/Pf;
        
    elseif strcmp(SBM.method,'SS')
        
        P_end = length(find(G<0))/SBM.NoS_SubSet{end};
        Pf = SBM.P0^(SBM.NoSubset-1) * P_end;
        Cov_Pf = nan;  %暂时不支持计算
        
    end
    
    if nn == 1
        SBM.Pf_predict = Pf;
        SBM.Cov_predict = Cov_Pf;
    elseif nn == 2
        SBM.Pf_ture = Pf;
        SBM.Cov_ture = Cov_Pf;
    end
    
end


%     ALSMTimeHis.nf(times) = length(find(G_predict<=0));
%     if strcmp(SBM.method,'MCS')
%         tempIndex = find(G_predict<=0); 
%         SBM.State(SBM.CPI) = 0;
%         SBM.State(SBM.CPI(tempIndex),1) = 1;
%         nf = sum(SBM.State);
%         ALSMTimeHis.Pf(times) = nf/SBM.NofSamples;
%         if ALSMTimeHis.false_ture_statistics == 1
%             nf_ture = length(find(SBM.G_ture<=0));
%             ALSMTimeHis.Pf_ture(times) = nf_ture/SBM.NofSamples;
%         end
%         
%     elseif strcmp(SBM.method,'IS')
%          tempIndex = find(G_predict<=0); 
%          SBM.State(SBM.CPI) = 0;
%          SBM.State(SBM.CPI(tempIndex),1) = 1;
%          ALSMTimeHis.Pf(times) = 1/SBM.NofSamples *sum(SBM.State.*SBM.Fpdf./SBM.Ppdf);
%          if ALSMTimeHis.false_ture_statistics == 1
%              Ture_state = zeros(SBM.NofSamples,1);
%              Ture_state(find(SBM.G_ture<=0)) = 1;
%              ALSMTimeHis.Pf_ture(times) = 1/SBM.NofSamples *sum(Ture_state.*SBM.Fpdf./SBM.Ppdf);
%          end
%     elseif strcmp(SBM.method,'SDMCS')        
% %              tempindex = find( ProSys.fname(CandidatePool,ProSys.fun_par) <=0 );
% %              SBM.State(tempindex) = 1;
%         tempIndex = find(G_predict<=0);
%         SBM.State(SBM.CPI(:)) = 0;
%         SBM.State(SBM.CPI(tempIndex),1) = 1;
%         for ni = 1:size(SBM.SubPoolIndex,2)
%             SubPoolState_predict = SBM.State(SBM.SubPoolIndex{ni},1);
%             Pi(ni) = length( find( SubPoolState_predict == 1) ) / SBM.NoS_ring(ni);
%             PFi(ni) = Pi(ni) * SBM.thtai(ni);
%             VarPi(ni) =  1/SBM.NoS_ring(ni) * Pi(ni)*(1-Pi(ni));
%         end
%         % Current fialure probability
%         ALSMTimeHis.Pf(times) =  sum(PFi); 
%         %------------------------------------------------------------------
%         if ALSMTimeHis.false_ture_statistics == 1
%             State_ture = zeros(SBM.NofSamples,1);
%             State_ture(find(SBM.G_ture<=0)) = 1;
%             SubPoolState_ture = State_ture(SBM.SubPoolIndex{ni},1);
%             Pi(ni) = length( find( SubPoolState_ture == 1) ) / SBM.NoS_ring(ni);
%             PFi(ni) = Pi(ni) * SBM.thtai(ni);       
%             ALSMTimeHis.Pf_ture(times) = sum(PFi);
%         end
%         %------------------------------------------------------------------
%     end
% end