function y=NdimParaboloid(X,fun_par)
% N dimension Paraboloid
% input m*n ��mΪ���ݵ�����nΪx��ά��
% fun_par =[] ������ 
% output yΪ������
  a = fun_par.a;
  b = fun_par.b;
  N=size(X,2);
  y = a*sum(X(:,2:end).^2,2)-b-X(:,1);
return