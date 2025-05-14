function y=DRONLO_LSF(u,fun_par)
% dynamic response of a non-linear oscillator
% input m*n ��mΪ���ݵ�����nΪx��ά��
% fun_par =[] ������ 
% output yΪ������

% ����׼��̬����uת��

    
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