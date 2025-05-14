function [G,T_stru] = StruDynamicNF_2RV(u,fun_par)
%% 转化数据
% 将标准正态变量u映射回原始空间
    [Ns,Ndim] = size(u);
if strcmp(fun_par.StadNrom,'Yes')     
    % 对数正态分布参数
     muY = fun_par.muX;
     sigmaY = fun_par.sigmaX; 
    % 相应正态分布参数
    deta = sigmaY./muY;
    sigmaX = sqrt(log(1+deta.^2));
    muX = log(muY)-sigmaX.^2/2;    
    % 标转正态分布转换为变量X
    X = zeros(Ns,Ndim);
    for n = 1:Ndim
        X(:,n) = u(:,n)*sigmaX(n) + muX(n);
    end
    % 正态变量X转化为变量Y
    u = exp(X);
else  % 不需转化  
    u = u;    
end

G = zeros(Ns,1);
T_stru  = zeros(Ns,3);

for iii = 1:Ns
%% variables
beta1 = u(iii,1);
beta2 = u(iii,2);

mm = [2.82,2.76*ones(1,18),2.92]*1e5;
kk = [1.007,1.358,1.338,1.33,1.228,1.248,1.237,1.223,1.208,1.187,1.093,...
    1.054,1.036,0.908,0.897,0.882,0.798,0.734,0.598,0.449]*1.0e8;

if fun_par.stru_type == 1
    DoF = 3;                           % degree of freedom
    m = mm(1:DoF)*1e5*beta1;
    k = kk(1:DoF)*1e8*beta2;
    kexi = 0.05;                       % Rayleigh damping
    
elseif fun_par.stru_type == 2
    DoF = 10;                          % degree of freedom
    m = mm(1:DoF)*1e5*beta1; 
    k = kk(1:DoF)*1e8*beta2;
    kexi = 0.05;                       % Rayleigh damping  
    
elseif fun_par.stru_type == 3
    DoF = 20;
    m = mm(1:DoF)*1e5*beta1; 
    k = kk(1:DoF)*1e8*beta2;
    kexi = 0.05;                       % Rayleigh damping 

elseif fun_par.stru_type == 4
    DoF = 5;
    m = [3.57,3.36,3.14,2.99,2.65]*1e5*beta1;
    k =  [0.8,0.77,0.74,0.68,0.62]*1e8*beta2;
    kexi = 0.05;                       % Rayleigh damping
end
%% structure model
    % matrix of K M C
    K=zeros(DoF);
    for i=1:DoF-1
        K(i,i)=k(i)+k(i+1);
        K(i,i+1)=-k(i+1);
        K(i+1,i)=-k(i+1);
    end
    K(DoF,DoF) = k(DoF);
    M = diag(m);
    % modal analyze
    [phi,lamda] = eig(K,M);
    % rearrange the phi
    [lamda,lamda_index] = sort(diag(lamda),'ascend');
    phi= phi(:,lamda_index);
    % normalization
    temp = sqrt(sum(phi.^2))'*ones(1,DoF);
    phi = phi./temp;
    wn = sqrt(diag(lamda));
    wi = wn(1,1); wj = wn(2,2);
    ab = 2*wi*wj/(wj^2-wi^2)*[wj , -wi; -1/wj, 1/wi] * [kexi;kexi];
    C = ab(1)*M+ab(2)*K;   
    T_stru(iii,1:3) = 2*pi./[wn(1,1),wn(2,2),wn(3,3)];    
%%  Exciation
    % EIcenrto earthquake wave input
    Gacce = load('elcentro.txt');
    PGA = 0.7; %8度多遇，m/s^2
    Gacce = Gacce'/max(abs(Gacce))*PGA;
    dt_inter = 0.02;
%   t = 0:dt_inter:(length(Gacce)-1)*dt_inter;    
%% Slover
    % Newmark-beta 
    F = m(:)*Gacce; 
    beta=1/6;gama=1/2;
    Y0=zeros(size(M,1),2);
    x = Newmark_dsmj(K,M,C,F,dt_inter,Y0,beta,gama);

%% Perform function
    % 最大位移
    if fun_par.indexType == 1
        G(iii,1) = fun_par.threshold - max(abs(x(DoF,:)));
    elseif fun_par.indexType == 2 
        % 最大层间位移角
        rel_D(1,:) = x(1,:);
        rel_D(2:DoF,:) = x(2:DoF,:)-x(1:DoF-1,:);
        G(iii) = fun_par.threshold -max(max(abs(rel_D)));
    end    
    if mod(iii,100)==0
        disp(['Running:StruDynamicNF_2RV--',num2str(iii),'/',num2str(Ns)])
    end
end