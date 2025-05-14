function y=eg28Function(u,fun_par)
% dynamic response of a non-linear oscillator
% input m*n ，m为数据点数；n为x的维数
% fun_par =[] 无需求 
% output y为列向量

% 将标准正态变量u映射回原始空间
if strcmp(fun_par.StadNrom,'Yes')    
    [Ns,Ndim] = size(u);  
    muX = fun_par.muX;
    sigmaX = fun_par.sigmaX; 
    X = zeros(Ns,Ndim);
    for n = 1:Ndim
        X(:,n) = u(:,n)*sigmaX(n) + muX(n);
    end
    
else  % 不需转化
    
    X = u;

end
      
% -------------------------------------------------------------------------
  y=4.0*X(:,2)-3.9998*X(:,3)+4*X(:,4)-X(:,1);
% -------------------------------------------------------------------------

return