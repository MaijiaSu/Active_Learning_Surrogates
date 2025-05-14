function [Pool,PoolPdf] = RNgeneratorV2(mu,sigma,TYPE,NEP,RvPar)
%% ˵��
% RNgeneratot��random number generator

% input  ����������� mu,sigama 
%        ����������� TYPE
%        ����������� NEPfunction [Pool,PoolPdf] = RNgeneratorV2(mu,sigma,TYPE,NEP,RvPar)
%% ˵��
% RNgeneratot��random number generator

% input  ����������� mu,sigama 
%        ����������� TYPE
%        ����������� NEP

% output ������������� Pool (��ģΪ NEP*Ndim)
%        ��������Ӧ�ĸ����ܶ� PoolPdf (��ģΪ NEP*1)

% ����������������������������������TYPE����������������������������������������
% TYPE = 1 : ��̫�ֲ�
% TYPE = 2 : ������̫�ֲ�
% TYPE = 3 : ����ֵI�ͷֲ�
% TYPE = 4 : ��СֵI�ͷֲ�
% TYPE = 5 : ���ȷֲ�
% TYPE =  'LHS_normal'    ���������������� - ��̬�ֲ�
% TYPE =  'LHS_uniform'   ���������������� - ���ȷֲ�
% TYPE =  'mvnrnd'        ����̬�ֲ������ϸ��ʷֲ�
% ������������������������������������������������������������������������������

% ������������������������������test example������������������������������������
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
% ��������������������������������������������������������������������������������

% ��ʼ������
    Ndim=length(mu);
    Pool = zeros(NEP,Ndim);    
    PoolPdf=zeros(NEP,Ndim);
%%  TYPE = num    

if ~ischar(TYPE)     %��TYPE�������ͷ��ַ�
    
    for i = 1:Ndim
        type = TYPE(i);
        switch type
            case 1   %��̬�������
                par = [mu(i),sigma(i)];
                Pool(:,i) = normrnd(par(1),par(2),NEP,1);
                PoolPdf(:,i) = normpdf(Pool(:,i),par(1),par(2));
                
            case 2   %������̬�������
                %����ת��
                deta = sigma(i)/mu(i);
                sLn = sqrt(log(1+deta^2));
                mLn = log(mu(i))-sLn^2/2;
                par = [mLn sLn];
                Pool(:,i) = lognrnd(par(1),par(2),NEP,1);
                PoolPdf(:,i) = lognpdf(Pool(:,i),par(1),par(2));
                
            case 3   %��ֵ1���������(����ֵ)
                %����ת��
                aEv = pi/sqrt(6)/sigma(i);
                uEv = psi(1)/aEv+mu(i); %-psi(1)Ϊŷ������
                par = [uEv aEv];
                Pool(:,i) = -1*evrnd(-par(1),1/par(2),NEP,1);
                PoolPdf(:,i) = evpdf(-Pool(:,i),-par(1),1/par(2));
                
            case 4    %��ֵ1���������(��Сֵ)
                %����ת��
                aEv = sqrt(6)*sigma(i)/pi;
                uEv = -psi(1)*aEv+mu(i); %-psi(1)Ϊŷ������
                par = [uEv aEv];
                Pool(:,i) = evrnd(par(1),par(2),NEP,1);
                PoolPdf(:,i) = evpdf(Pool(:,i),par(1),par(2));
                
            case 5  %���ȷֲ�
                a = mu(i)-sqrt(3)*sigma(i);
                b = mu(i)+sqrt(3)*sigma(i);
                Pool(:,i) = a+(b-a).*rand(NEP,1);
                PoolPdf(:,i) = 1/(b-a)*ones(NEP,1);
                
            otherwise
                error('type�������ͳ�����Χ')
        end
    end
    
    % �������ϸ����ܶ�
    PoolPdf = prod(PoolPdf')';
    
end
%% TYPE = char

if ischar(TYPE)  % ��TYPE��������Ϊ�ַ�
    
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
        covX = sigma;           % ֱ������Э����
        Pool = mvnrnd(mu,covX,NEP);
        PoolPdf = mvnpdf(Pool,mu,covX);
    
    else
                error('type�������ͳ�����Χ')
    end
    
elseif strcmp(TYPE,'MCMC')
        covX = sigma;           % ֱ������Э����
        Pool = mvnrnd(mu,covX,NEP);
        PoolPdf = mvnpdf(Pool,mu,covX);
        ExpextedNofSamples = 1e5;
        ThinChain = 5;
        mccount = ExpextedNofSamples*ThinChain*1.25;
        Dim = 2;
        NofChain = 100;  %number of Chain
        logPfuns = @(x) log(TargetPdf(x));
        [Pool,logP]=gwmcmc(randn(Dim,NofChain),logPfuns,mccount,'ThinChain',ThinChain);
        Pool(:,:,1:floor(size(Pool,3)*0.2))=[]; %remove 20% as burn-in
        Pool=Pool(:,:)';
        PoolPdf
end

return
    

% output ������������� Pool (��ģΪ NEP*Ndim)
%        ��������Ӧ�ĸ����ܶ� PoolPdf (��ģΪ NEP*1)

% ����������������������������������TYPE����������������������������������������
% TYPE = 1 : ��̫�ֲ�
% TYPE = 2 : ������̫�ֲ�
% TYPE = 3 : ����ֵI�ͷֲ�
% TYPE = 4 : ��СֵI�ͷֲ�
% TYPE = 5 : ���ȷֲ�
% TYPE =  'LHS_normal'    ���������������� - ��̬�ֲ�
% TYPE =  'LHS_uniform'   ���������������� - ���ȷֲ�
% TYPE =  'mvnrnd'        ����̬�ֲ������ϸ��ʷֲ�
% ������������������������������������������������������������������������������

% ������������������������������test example������������������������������������
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
% ��������������������������������������������������������������������������������

% ��ʼ������
    Ndim=length(mu);
    Pool = zeros(NEP,Ndim);    
    PoolPdf=zeros(NEP,Ndim);
%%  TYPE = num    

if ~ischar(TYPE)     %��TYPE�������ͷ��ַ�
    
    for i = 1:Ndim
        type = TYPE(i);
        switch type
            case 1   %��̬�������
                par = [mu(i),sigma(i)];
                Pool(:,i) = normrnd(par(1),par(2),NEP,1);
                PoolPdf(:,i) = normpdf(Pool(:,i),par(1),par(2));
                
            case 2   %������̬�������
                %����ת��
                deta = sigma(i)/mu(i);
                sLn = sqrt(log(1+deta^2));
                mLn = log(mu(i))-sLn^2/2;
                par = [mLn sLn];
                Pool(:,i) = lognrnd(par(1),par(2),NEP,1);
                PoolPdf(:,i) = lognpdf(Pool(:,i),par(1),par(2));
                
            case 3   %��ֵ1���������(����ֵ)
                %����ת��
                aEv = pi/sqrt(6)/sigma(i);
                uEv = psi(1)/aEv+mu(i); %-psi(1)Ϊŷ������
                par = [uEv aEv];
                Pool(:,i) = -1*evrnd(-par(1),1/par(2),NEP,1);
                PoolPdf(:,i) = evpdf(-Pool(:,i),-par(1),1/par(2));
                
            case 4    %��ֵ1���������(��Сֵ)
                %����ת��
                aEv = sqrt(6)*sigma(i)/pi;
                uEv = -psi(1)*aEv+mu(i); %-psi(1)Ϊŷ������
                par = [uEv aEv];
                Pool(:,i) = evrnd(par(1),par(2),NEP,1);
                PoolPdf(:,i) = evpdf(Pool(:,i),par(1),par(2));
                
            case 5  %���ȷֲ�
                a = mu(i)-sqrt(3)*sigma(i);
                b = mu(i)+sqrt(3)*sigma(i);
                Pool(:,i) = a+(b-a).*rand(NEP,1);
                PoolPdf(:,i) = 1/(b-a)*ones(NEP,1);
                
            otherwise
                error('type�������ͳ�����Χ')
        end
    end
    
    % �������ϸ����ܶ�
    PoolPdf = prod(PoolPdf')';
    
end
%% TYPE = char

if ischar(TYPE)  % ��TYPE��������Ϊ�ַ�
    
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
        covX = sigma;           % ֱ������Э����
        Pool = mvnrnd(mu,covX,NEP);
        PoolPdf = mvnpdf(Pool,mu,covX);
    
    else
                error('type�������ͳ�����Χ')
    end
    
elseif strcmp(TYPE,'MCMC')
        covX = sigma;           % ֱ������Э����
        Pool = mvnrnd(mu,covX,NEP);
        PoolPdf = mvnpdf(Pool,mu,covX);
        ExpextedNofSamples = 1e5;
        ThinChain = 5;
        mccount = ExpextedNofSamples*ThinChain*1.25;
        Dim = 2;
        NofChain = 100;  %number of Chain
        logPfuns = @(x) log(TargetPdf(x));
        [Pool,logP]=gwmcmc(randn(Dim,NofChain),logPfuns,mccount,'ThinChain',ThinChain);
        Pool(:,:,1:floor(size(Pool,3)*0.2))=[]; %remove 20% as burn-in
        Pool=Pool(:,:)';
        PoolPdf
end

return
    