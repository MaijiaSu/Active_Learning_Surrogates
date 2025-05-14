function y = CLSS_2doscillator(u,fun_par)
% dynamic response of a non-linear oscillator
% input m*n ��mΪ���ݵ�����nΪx��ά��
% fun_par =[] ������ 
% output yΪ������

% ����׼��̬����uת��
    
%     mu_Fs = 15; 
% 
%     mu_Fs = fun_par.mu_Fs;    
%     
%     [Ns,Ndim] = size(u);
%     
%     % ������̬�ֲ�����
%     muY=[1.5 0.01 1 0.01 0.05 0.02 mu_Fs 100];
%     cov = [10 10 20 20 40 50 10 10]/100;
%     sigmaY = muY.*cov;
%     
%     % ��Ӧ��̬�ֲ�����
%     deta = sigmaY./muY;
%     sigmaX = sqrt(log(1+deta.^2));
%     muX = log(muY)-sigmaX.^2/2;
%     
%     % ��ת��̬�ֲ�ת��Ϊ����X
%     X = zeros(Ns,Ndim);
%     for n = 1:Ndim
%         X(:,n) = u(:,n)*sigmaX(n) + muX(n);
%     end
%     
%     % ��̬����Xת��Ϊ����Y
%     Y = exp(X);
      [Ns,Ndim] = size(u);
        if strcmp(fun_par.StadNrom,'Yes')
            % ������̬�ֲ�����
            muY = fun_par.muX;
            sigmaY = fun_par.sigmaX;
            % ��Ӧ��̬�ֲ�����
            deta = sigmaY./muY;
            sigmaX = sqrt(log(1+deta.^2));
            muX = log(muY)-sigmaX.^2/2;
            % ��ת��̬�ֲ�ת��Ϊ����X
            X = zeros(Ns,Ndim);
            for n = 1:Ndim
                X(:,n) = u(:,n)*sigmaX(n) + muX(n);
            end
            % ��̬����Xת��Ϊ����Y
            Y = exp(X);

        else  % ����ת��

            Y = u;

        end    
    
    mp=Y(:,1); ms=Y(:,2); kp=Y(:,3); ks=Y(:,4);
    kexip=Y(:,5); kexis=Y(:,6); Fs=Y(:,7); S0=Y(:,8);
  
    % ���ܺ���
    wp = sqrt(kp./mp); 
    ws = sqrt(ks./ms);
    wa = (wp+ws)*0.5;
    kexia = (kexip+kexis)*0.5;
    gama = ms./mp;
    thta = (wp-ws)./wa;
    
    temp1 = pi*S0./(4*kexis.*ws.^3);
    temp2 =  kexia.*kexis./(kexip.*kexis.*(4*kexia.^2+thta.^2)+gama.*kexia.^2);
    temp3 = (kexip.*wp.^3+kexis.*ws.^3).*wp./(4*kexia.*wa.^4);
    
    p = 3;
    y = Fs - p*ks .* sqrt(temp1.*temp2.*temp3);

return