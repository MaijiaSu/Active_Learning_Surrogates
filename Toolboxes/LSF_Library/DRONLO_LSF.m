function y=DRONLO_LSF(u,fun_par)
% dynamic response of a non-linear oscillator
% input m*n ，m为数据点数；n为x的维数
% fun_par =[] 无需求 
% output y为列向量

% 将标准正态变量u转化

    
 [Ns,Ndim] = size(u);
  X = u;
  
  if Ns==1
      m=X(1);c1=X(2);c2=X(3);
      r=X(4);F1=X(5);t1=X(6);
  else
      m=X(:,1);c1=X(:,2);c2=X(:,3);
      r=X(:,4);F1=X(:,5);t1=X(:,6);
  end
  
  w0=sqrt((c1+c2)./m);
  y=3*r-abs(2*F1./m./w0.^2.*sin(w0.^1.*t1/2));


return