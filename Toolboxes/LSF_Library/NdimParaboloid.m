function y=NdimParaboloid(X,fun_par)
% N dimension Paraboloid
% input m*n ，m为数据点数；n为x的维数
% fun_par =[] 无需求 
% output y为列向量
  a = fun_par.a;
  b = fun_par.b;
  N=size(X,2);
  y = a*sum(X(:,2:end).^2,2)-b-X(:,1);
return