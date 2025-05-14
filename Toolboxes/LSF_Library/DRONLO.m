function y=DRONLO(u,fun_par)
% dynamic response of a non-linear oscillator
% input m*n ��mΪ���ݵ�����nΪx��ά��
% fun_par =[] ������ 
% output yΪ������

% ����׼��̬����uӳ���ԭʼ�ռ�
if strcmp(fun_par.StadNrom,'Yes')    
    [Ns,Ndim] = size(u);  
    muX = fun_par.muX;
    sigmaX = fun_par.sigmaX; 
    X = zeros(Ns,Ndim);
    for n = 1:Ndim
        X(:,n) = u(:,n)*sigmaX(n) + muX(n);
    end
    
else  % ����ת��
    
    X = u;

end
      
% -------------------------------------------------------------------------
  m=X(:,1);c1=X(:,2);c2=X(:,3);
  r=X(:,4);F1=X(:,5);t1=X(:,6); 
  w0=sqrt((c1+c2)./m);
  y=3*r-abs(2*F1./m./w0.^2.*sin(w0.^1.*t1/2));
% -------------------------------------------------------------------------

return