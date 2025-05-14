function [Pool,PoolPdf] = RNgeneratorV4(ProSys,NEP)
% RNgeneratot：random number generator
% input  A stru varibales ProSys

% output Sampling Pool (Size: NEP*Ndim)
%        The PDF of Pool PoolPdf (Size: NEP*1)

% ———————————————————Type————————————————————————
% TYPE = 1 : Normal distribution
% TYPE = 2 : Log normal
% TYPE = 3 : Gumbel-I Max
% TYPE = 4 : Gumbel-I Min
% TYPE = 5 : Uniform
% TYPE = 6 : Left truancated normal
% TYPE = 7 : Weibull 
% TYPE =  'LHS_normal'    ：LHS -normal
% TYPE =  'LHS_uniform'   ：LHS - unifrom
% TYPE =  'mvnrnd'        ：multi-normal diribution
% TYPE =  'MCMC'        ：Sampling from MCMC
% ———————————————————————————————————————————————

%% Initializaiton
    Ndim = ProSys.Ndim;
    Pool = zeros(NEP,Ndim);    
    TYPE = ProSys.Distri;
    PoolPdf = [];
    if ~ischar(TYPE)&&~(strcmp(TYPE,'MCMC'))
        sigma = ProSys.sigmaX;
        mu = ProSys.muX;
    end

%% Compute the Parametrs of Random model
if ~ischar(TYPE)     % non char input

    if isfield(ProSys,'RVPar')
        IsInputPar = ProSys.RVPar.isInputPar;
    else
        IsInputPar = zeros(1,Ndim);
    end

    for i = 1:Ndim
        type = TYPE(i);

        if IsInputPar(i)==1
            PAR{i} = [ProSys.RVPar.P1(i),ProSys.RVPar.P2(i)];
        else
            switch type
                case 1
                    PAR{i} = [mu(i),sigma(i)];
                case 2
                    deta = sigma(i)/mu(i);
                    sLn = sqrt(log(1+deta^2));
                    mLn = log(mu(i))-sLn^2/2;
                    PAR{i} = [mLn sLn];
                case 3
                    aEv = pi/sqrt(6)/sigma(i);
                    uEv = psi(1)/aEv+mu(i); %-psi(1)为欧拉常数
                    PAR{i} = [uEv aEv];
                case 4
                    aEv = sqrt(6)*sigma(i)/pi;
                    uEv = -psi(1)*aEv+mu(i); %-psi(1)为欧拉常数
                    PAR{i} = [uEv aEv];
                case 5
                    a = mu(i)-sqrt(3)*sigma(i);
                    b = mu(i)+sqrt(3)*sigma(i);
                    PAR{i} = [a,b];
                case 6  % left truncated normal
                    PAR{i} = [mu(i),sigma(i)];
                case 7                   
                    error('how to calculate the paramters')
                    PAR{i} = [uEv aEv];

                otherwise
                    error('Input type is out of range ')
            end
        end
        
    end
end

%%  TYPE = num    

if ~ischar(TYPE)     % non char input
    for i = 1:Ndim
        type = TYPE(i);
        switch type
            case 1   
                par = PAR{i};
                Pool(:,i) = normrnd(par(1),par(2),NEP,1);
                                
            case 2                
                par = PAR{i};
                Pool(:,i) = lognrnd(par(1),par(2),NEP,1);               
                
            case 3   
                par = PAR{i};
                Pool(:,i) = -1*evrnd(-par(1),1/par(2),NEP,1);   

            case 4    
                par = PAR{i};
                Pool(:,i) = evrnd(par(1),par(2),NEP,1);
                               
            case 5  
                par = PAR{i};
                a = par(1);
                b = par(2);
                Pool(:,i) = a+(b-a).*rand(NEP,1);
                               
            case 6  % left truncated normal
                PoolXi = [];
                par = PAR{i};
                while numel(PoolXi>0)<NEP
                    NofTempX = NEP-numel(PoolXi);
                    tempX = normrnd(par(1),par(2),NofTempX,1);
                    ID = find(tempX<=0);
                    tempX(ID) = [];
                    PoolXi = [PoolXi;tempX];
                end
                Pool(:,i) = PoolXi;

            case 7  % left truncated normal
                par = PAR{i};
                Pool(:,i) = wblrnd(par(1),par(2),NEP,1);
                
            otherwise
                error('input type is out of range')
        end
    end
    
    % Jointed PDF
end

%% TYPE = char

if ischar(TYPE)  % 若TYPE数据类型为字符
    
    if strcmp(TYPE,'LHS')  % lhsnorm
        Pool = lhsnorm(mu(:),sigma,NEP); 
%         PoolPdf = normpdf(Pool,mu,sigma);
         
     
    elseif strcmp(TYPE,'LHS_uniform')
        Pool = lhsdesign(NEP,Ndim);      
        for i = 1:Ndim
            a(i) = mu(i)-sqrt(3)*sigma(i);
            b(i) = mu(i)+sqrt(3)*sigma(i);
            Pool(:,i) = a(i)+(b(i)-a(i))*Pool(:,i);
        end        
        

    elseif strcmp(TYPE,'mvnrnd')
        covX = ProSys.sigmaX;           % 直接输入协方差
        mu = ProSys.muX;
        Pool = mvnrnd(mu,covX,NEP);
        
    
    elseif strcmp(TYPE,'MCMC')
        
        ExpextedNofSamples = NEP;
        ThinChain = ProSys.MCMC.ThinChain;
        mccount = ExpextedNofSamples*ThinChain*1.25;  % 1.25 = 1/(1-20%), 20% as burn-in
        NofChain = ProSys.MCMC.NofChain;  %number of Chain
        logPfuns = @(x) log(ProSys.MCMC.TargetPdf(x));

        [Pool,logP]=gwmcmc(randn(Ndim,NofChain),logPfuns,mccount,'ThinChain',ThinChain);
        Pool(:,:,1:floor(size(Pool,3)*0.2))=[];  %remove 20% as burn-in
        Pool=Pool(:,:)';
        logP(:,:,1:floor(size(logP,3)*0.2))=[];  %remove 20% as burn-in
        logP=logP(:,:)';
    else
        error('input type is out of range')
    end
end

return