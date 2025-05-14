% clc,clear,close all
% TestExample = 'eg1';                           % select the example from LSF_library, LimtStateFunction_select.m;
% SBM.method = 'MCS_UQ';
% LimtStateFunction_select
% -------------------------------------------------------------------------
%     for nn = [15,50,100,300,1000]
%             ALSMPar.IniDoE.VM = 1.5;
%             ALSMPar.IniDoE.N0 = nn;
% %            DoE.X = RNgeneratorV2(ProSys.muX,ALSMPar.IniDoE.VM*ProSys.sigmaX,...
% %                 'LHS_uniform',ALSMPar.IniDoE.N0);
% %             DoE.X = uq_getSample(myInput,ALSMPar.IniDoE.N0,'LHS');
%      % method #3
%             a = min(CandidatePool);
%             b = max(CandidatePool);
%             DoE.X = lhsdesign(ALSMPar.IniDoE.N0,ProSys.Ndim);           
%             for i = 1:ProSys.Ndim
%                 DoE.X(:,i) = a(i)+(b(i)-a(i))*DoE.X(:,i);
%             end          
%             DoE.Y = ProSys.fname(DoE.X,ProSys.fun_par);    
%             filename = [TestExample,'_','DoE',num2str(nn)]
%             save(filename,'DoE')
%     end 
% -------------------------------------------------------------------------