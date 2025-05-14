function [Pool,PoolPdf] = RNgeneratorV3(ProSys,NEP)
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
% TYPE =  'LHS_normal'    ：LHS -normal
% TYPE =  'LHS_uniform'   ：LHS - unifrom
% TYPE =  'mvnrnd'        ：multi-normal diribution
% TYPE =  'MCMC'        ：Sampling from MCMC
% ———————————————————————————————————————————————


% Initializaiton
    Ndim = ProSys.Ndim;
    Pool = zeros(NEP,Ndim);    
    PoolPdf = zeros(NEP,Ndim);
    TYPE = ProSys.Distri;
    if ~ischar(TYPE)&&~(strcmp(TYPE,'MCMC'))
        sigma = ProSys.sigmaX;
        mu = ProSys.muX;
    end

%%  TYPE = num    

if ~ischar(TYPE)     % non char input
    
    for i = 1:Ndim
        type = TYPE(i);
        switch type
            case 1   
                par = [mu(i),sigma(i)];
                Pool(:,i) = normrnd(par(1),par(2),NEP,1);
                PoolPdf(:,i) = normpdf(Pool(:,i),par(1),par(2));
                
            case 2                
                deta = sigma(i)/mu(i);
                sLn = sqrt(log(1+deta^2));
                mLn = log(mu(i))-sLn^2/2;
                par = [mLn sLn];
                Pool(:,i) = lognrnd(par(1),par(2),NEP,1);
                PoolPdf(:,i) = lognpdf(Pool(:,i),par(1),par(2));
                
            case 3   
                aEv = pi/sqrt(6)/sigma(i);
                uEv = psi(1)/aEv+mu(i); %-psi(1)为欧拉常数
                par = [uEv aEv];
                Pool(:,i) = -1*evrnd(-par(1),1/par(2),NEP,1);
                PoolPdf(:,i) = evpdf(-Pool(:,i),-par(1),1/par(2));
                
            case 4    
                aEv = sqrt(6)*sigma(i)/pi;
                uEv = -psi(1)*aEv+mu(i); %-psi(1)为欧拉常数
                par = [uEv aEv];
                Pool(:,i) = evrnd(par(1),par(2),NEP,1);
                PoolPdf(:,i) = evpdf(Pool(:,i),par(1),par(2));
                
            case 5  
                a = mu(i)-sqrt(3)*sigma(i);
                b = mu(i)+sqrt(3)*sigma(i);
                Pool(:,i) = a+(b-a).*rand(NEP,1);
                PoolPdf(:,i) = 1/(b-a)*ones(NEP,1);
                
            case 6  % left truncated normal
                PoolXi = [];
                par = [mu(i),sigma(i)];
                while numel(PoolXi>0)<NEP
                    NofTempX = NEP-numel(PoolXi);
                    tempX = normrnd(par(1),par(2),NofTempX,1);
                    ID = find(tempX<=0);
                    tempX(ID) = [];
                    PoolXi = [PoolXi;tempX];
                end
                Pool(:,i) = PoolXi;
                PoolPdf(:,i) = normpdf(Pool(:,i),par(1),par(2));

            otherwise
                error('type输入类型超出范围')
        end
    end
    
    % 计算联合概率密度
    PoolPdf = prod(PoolPdf')';  
end

%% TYPE = char

if ischar(TYPE)  % 若TYPE数据类型为字符
    
    if strcmp(TYPE,'LHS')  % lhsnorm
        Pool = lhsnorm(mu(:),sigma,NEP); 
%         PoolPdf = normpdf(Pool,mu,sigma);
         PoolPdf = [];
     
    elseif strcmp(TYPE,'LHS_uniform')
        Pool = lhsdesign(NEP,Ndim);      
        for i = 1:Ndim
            a(i) = mu(i)-sqrt(3)*sigma(i);
            b(i) = mu(i)+sqrt(3)*sigma(i);
            Pool(:,i) = a(i)+(b(i)-a(i))*Pool(:,i);
        end        
        PoolPdf = 1/((b-a)*(b-a)')*ones(NEP,1);

    elseif strcmp(TYPE,'mvnrnd')
        covX = ProSys.sigmaX;           % 直接输入协方差
        mu = ProSys.muX;
        Pool = mvnrnd(mu,covX,NEP);
        PoolPdf = mvnpdf(Pool,mu,covX);
    
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
        PoolPdf = exp(logP);
    else
        error('type输入类型超出范围')
    end
end

return
    

% output 随机变量样本池 Pool (规模为 NEP*Ndim)
%        样本池相应的概率密度 PoolPdf (规模为 NEP*1)

% ——————————————可输入TYPE————————————————————
% TYPE = 1 : 正太分布
% TYPE = 2 : 对数正太分布
% TYPE = 3 : 极大值I型分布
% TYPE = 4 : 极小值I型分布
% TYPE = 5 : 均匀分布
% TYPE =  'LHS_normal'    ：拉丁超立方抽样 - 正态分布
% TYPE =  'LHS_uniform'   ：拉丁超立方抽样 - 均匀分布
% TYPE =  'mvnrnd'        ：正态分布的联合概率分布
% ———————————————————————————————————————

% ———————————————test example——————————————————
% eg1
% niu=[0.36 0.18 20];     
% sigma=[0.036 0.018 5.0];
% TYPE=[1 2 4]; 
% NEP=1000;
% [Pool,PoolPdf] = RNgeneratorV2(mu,sigma,TYPE,NEP)
% eg2
% mu = [0.36 0.18];     
% sigma = [0.036 0.018];
% TYPE = ['LHS_uniform']; 
% NEP=10;
% [Pool,PoolPdf] = RNgeneratorV2(mu,sigma,TYPE,NEP)
% plot(Pool(:,1),Pool(:,2),'ko')
% ————————————————————————————————————————

% 初始化参数
%     Ndim=length(mu);
%     Pool = zeros(NEP,Ndim);    
%     PoolPdf=zeros(NEP,Ndim);
% %%  TYPE = num    
% 
% if ~ischar(TYPE)     %若TYPE数据类型非字符
%     
%     for i = 1:Ndim
%         type = TYPE(i);
%         switch type
%             case 1   %正态随机变量
%                 par = [mu(i),sigma(i)];
%                 Pool(:,i) = normrnd(par(1),par(2),NEP,1);
%                 PoolPdf(:,i) = normpdf(Pool(:,i),par(1),par(2));
%                 
%             case 2   %对数正态随机变量
%                 %参数转化
%                 deta = sigma(i)/mu(i);
%                 sLn = sqrt(log(1+deta^2));
%                 mLn = log(mu(i))-sLn^2/2;
%                 par = [mLn sLn];
%                 Pool(:,i) = lognrnd(par(1),par(2),NEP,1);
%                 PoolPdf(:,i) = lognpdf(Pool(:,i),par(1),par(2));
%                 
%             case 3   %极值1型随机变量(极大值)
%                 %参数转化
%                 aEv = pi/sqrt(6)/sigma(i);
%                 uEv = psi(1)/aEv+mu(i); %-psi(1)为欧拉常数
%                 par = [uEv aEv];
%                 Pool(:,i) = -1*evrnd(-par(1),1/par(2),NEP,1);
%                 PoolPdf(:,i) = evpdf(-Pool(:,i),-par(1),1/par(2));
%                 
%             case 4    %极值1型随机变量(极小值)
%                 %参数转化
%                 aEv = sqrt(6)*sigma(i)/pi;
%                 uEv = -psi(1)*aEv+mu(i); %-psi(1)为欧拉常数
%                 par = [uEv aEv];
%                 Pool(:,i) = evrnd(par(1),par(2),NEP,1);
%                 PoolPdf(:,i) = evpdf(Pool(:,i),par(1),par(2));
%                 
%             case 5  %均匀分布
%                 a = mu(i)-sqrt(3)*sigma(i);
%                 b = mu(i)+sqrt(3)*sigma(i);
%                 Pool(:,i) = a+(b-a).*rand(NEP,1);
%                 PoolPdf(:,i) = 1/(b-a)*ones(NEP,1);
%                 
%             otherwise
%                 error('type输入类型超出范围')
%         end
%     end
%     
%     % 计算联合概率密度
%     PoolPdf = prod(PoolPdf')';
%     
% end
% %% TYPE = char
% 
% if ischar(TYPE)  % 若TYPE数据类型为字符
%     
%     if strcmp(TYPE,'LHS')  % lhsnorm
%         Pool = lhsnorm(mu(:),sigma,NEP); 
% %         PoolPdf = normpdf(Pool,mu,sigma);
%          PoolPdf = [];
%      
%     elseif strcmp(TYPE,'LHS_uniform')
%         Pool = lhsdesign(NEP,Ndim);      
%         for i = 1:Ndim
%             a(i) = mu(i)-sqrt(3)*sigma(i);
%             b(i) = mu(i)+sqrt(3)*sigma(i);
%             Pool(:,i) = a(i)+(b(i)-a(i))*Pool(:,i);
%         end        
%         PoolPdf = 1/((b-a)*(b-a)')*ones(NEP,1);
% 
%     elseif strcmp(TYPE,'mvnrnd')
%         covX = sigma;           % 直接输入协方差
%         Pool = mvnrnd(mu,covX,NEP);
%         PoolPdf = mvnpdf(Pool,mu,covX);
% 
%     elseif strcmp(TYPE,'MCMC')
%         ExpextedNofSamples = 1e5;
%         ThinChain = 5;
%         mccount = ExpextedNofSamples*ThinChain*1.25;
%         NofChain = 100;  %number of Chain
%         logPfuns = @(x) log(TargetPdf(x));
%         [Pool,logP]=gwmcmc(randn(Ndim,NofChain),logPfuns,mccount,'ThinChain',ThinChain);
%         Pool(:,:,1:floor(size(Pool,3)*0.2))=[]; %remove 20% as burn-in
%         Pool=Pool(:,:)';
% 
%         logP(:,:,1:floor(size(logP,3)*0.2))=[];
%         logP=logP(:,:)';
%         PoolPdf = exp(logP);
%     else
%                 error('type输入类型超出范围')
%     end
%     
% 
% end
% 
% return
    