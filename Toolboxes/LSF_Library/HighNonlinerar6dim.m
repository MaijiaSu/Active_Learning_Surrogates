function y = HighNonlinerar6dim(u,fun_par)
% HighNonlinerar6dim
% input m*n ��mΪ���ݵ�����nΪx��ά��
% fun_par =[] ������ 
% output yΪ������

% ����׼��̬����uת��
    
    [Ns,Ndim] = size(u);
    
%     % ������̬�ֲ�����
    muY =[120 120 120 120 50 40];
    sigmaY = [12 12 12 12 15 12];
    
    % ��Ӧ��̬�ֲ�����
    deta = sigmaY./muY;
    sigmaX = sqrt(log(1+deta.^2));
    muX = log(muY)-sigmaX.^2/2;
    
    % ��ת��̬�ֲ�ת��Ϊ����X
    X = zeros(Ns,Ndim);
    for n = 1:Ndim
        X(:,n) = u(:,n)*sigmaX(n) + muX(n);
    end
    
    % ��̬����Xת��Ϊ������̬����Y
    Y = exp(X);
    
    x1=Y(:,1);x2=Y(:,2);x3=Y(:,3);
    x4=Y(:,4);x5=Y(:,5);x6=Y(:,6);
  
    % ���ܺ���
  
    y = x1 + 2*x2 + 2*x3 + x4 - 5*x5 - 5*x6 + 0.001*sum(sin(100*Y'))';

return