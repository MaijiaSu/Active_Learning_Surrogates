function y=eg28Function(u,fun_par)
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
  y=4.0*X(:,2)-3.9998*X(:,3)+4*X(:,4)-X(:,1);
% -------------------------------------------------------------------------

return