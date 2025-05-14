function PoolPdf = MyPoolJointPdf(mu,sigma,TYPE,Pool)
%% 说明
% 2020.04.28 (待完善，目前相当于只编了独立正态分布的情况)

% MyCdf： 得到联合概率密度函数 

% input  随机变量参数 mu,sigama 
%        随机变量类型 TYPE
%        样本池

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
    Ndim=length(mu);   
    NEP = size(Pool,1);
    PoolPdf=zeros(NEP,Ndim);
%%  TYPE = num    

if ~ischar(TYPE)     %若TYPE数据类型非字符
    
    for i = 1:Ndim
        type = TYPE(i);
        switch type
            case 1   %正态随机变量
                par = [mu(i),sigma(i)];
                PoolPdf(:,i) = normpdf(Pool(:,i),par(1),par(2));
                
            case 2   %对数正态随机变量
                %参数转化
                deta = sigma(i)/mu(i);
                sLn = sqrt(log(1+deta^2));
                mLn = log(mu(i))-sLn^2/2;
                par = [mLn sLn];
                Pool(:,i) = lognrnd(par(1),par(2),NEP,1);
                PoolPdf(:,i) = lognpdf(Pool(:,i),par(1),par(2));
                
            case 3   %极值1型随机变量(极大值)
                %参数转化
                aEv = pi/sqrt(6)/sigma(i);
                uEv = psi(1)/aEv+mu(i); %-psi(1)为欧拉常数
                par = [uEv aEv];
                Pool(:,i) = -1*evrnd(-par(1),1/par(2),NEP,1);
                PoolPdf(:,i) = evpdf(-Pool(:,i),-par(1),1/par(2));
                
            case 4    %极值1型随机变量(极小值)
                %参数转化
                aEv = sqrt(6)*sigma(i)/pi;
                uEv = -psi(1)*aEv+mu(i); %-psi(1)为欧拉常数
                par = [uEv aEv];
                Pool(:,i) = evrnd(par(1),par(2),NEP,1);
                PoolPdf(:,i) = evpdf(Pool(:,i),par(1),par(2));
                
            case 5  %均匀分布
                a = mu(i)-sqrt(3)*sigma(i);
                b = mu(i)+sqrt(3)*sigma(i);
                Pool(:,i) = a+(b-a).*rand(NEP,1);
                PoolPdf(:,i) = 1/(b-a)*ones(NEP,1);
                
            otherwise
                error('type输入类型超出范围')
        end
    end
    
    % 计算联合概率密度
    PoolPdf = prod(PoolPdf')';
    
end