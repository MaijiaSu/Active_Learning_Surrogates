function y = HighNonlinerar6dim(u,fun_par)
% HighNonlinerar6dim
% input m*n ，m为数据点数；n为x的维数
% fun_par =[] 无需求 
% output y为列向量

% 将标准正态变量u转化
    
    [Ns,Ndim] = size(u);
    
%     % 对数正态分布参数
    muY =[120 120 120 120 50 40];
    sigmaY = [12 12 12 12 15 12];
    
    % 相应正态分布参数
    deta = sigmaY./muY;
    sigmaX = sqrt(log(1+deta.^2));
    muX = log(muY)-sigmaX.^2/2;
    
    % 标转正态分布转换为变量X
    X = zeros(Ns,Ndim);
    for n = 1:Ndim
        X(:,n) = u(:,n)*sigmaX(n) + muX(n);
    end
    
    % 正态变量X转化为对数正态变量Y
    Y = exp(X);
    
    x1=Y(:,1);x2=Y(:,2);x3=Y(:,3);
    x4=Y(:,4);x5=Y(:,5);x6=Y(:,6);
  
    % 功能函数
  
    y = x1 + 2*x2 + 2*x3 + x4 - 5*x5 - 5*x6 + 0.001*sum(sin(100*Y'))';

return