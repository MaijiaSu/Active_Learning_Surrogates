function UqlabInput = MyProSys2Uqlab(ProSys)
% Parameters
if isfield(ProSys,'RVPar')
    IsInputPar = ProSys.RVPar.isInputPar;
else
    IsInputPar = zeros(1,ProSys.Ndim);
end

    for ii = 1:ProSys.Ndim
        InputOpts.Marginals(ii).Name = ['X',num2str(ii)]; 
        % Random variables type
        if  ProSys.Distri(ii)== 1|| ProSys.Distri(ii)==6
            InputOpts.Marginals(ii).Type = 'Gaussian';
        elseif ProSys.Distri(ii)== 2
            InputOpts.Marginals(ii).Type = 'Lognormal';
        elseif ProSys.Distri(ii)== 3
            InputOpts.Marginals(ii).Type = 'Gumbel';
        elseif ProSys.Distri(ii)== 4
            InputOpts.Marginals(ii).Type = 'Gumbelmin';
        elseif ProSys.Distri(ii)== 5
            InputOpts.Marginals(ii).Type = 'Uniform';
        elseif  ProSys.Distri(ii)== 1|| ProSys.Distri(ii)==7
            InputOpts.Marginals(ii).Type = 'Weibull';
        elseif strcmp(ProSys.Distri,'mvnrnd')||strcmp(ProSys.Distri,'MCMC')
            InputOpts.Marginals(ii).Type = 'Gaussian';
        else
            error('Need to update the code')
        end


         if IsInputPar(ii)==1
             InputOpts.Marginals(ii).Parameters = [ProSys.RVPar.P1(ii),ProSys.RVPar.P2(ii)];
         else
            InputOpts.Marginals(ii).Moments = [ProSys.muX(ii),ProSys.sigmaX(ii)];
         end
       
    end
    % Dependency
%     InputOpts.Copula.Type = 'Independent';

    UqlabInput = uq_createInput(InputOpts);
end