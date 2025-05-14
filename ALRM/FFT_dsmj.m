function [P1 ph iFFT_error f]=FFT_dsmj(X,T,PTO)
%函数功能：采用FFT得到输入信号的频率与
%input  时域信号X，信号采用周期T,绘图参数
%output 单边幅值谱及单边相位谱，以及逆向傅里叶变化误差
%编于 2019/5/21 dsmj
%% 初始化
% % input
% Original signal X(t)
% T=0.02;                % Sampling period  
% 初始化
    L=length(X);          % Length of signal
    Fs = 1/T;             % Sampling frequency   
    t = (0:L-1)*T;        % Time vector
    f = Fs*(0:(L/2))/L;   % frequency domain
%% 算例
% clc,clear,close all
% PTO=[1 1]
% eg1
% Fs = 1000;            % Sampling frequency                    
% T = 1/Fs;             % Sampling period       
% L = 1000;             % Length of signal
% t = (0:L-1)*T;        % Time vector
% X=0.7*sin(2*pi*50*t) - sin(2*pi*150*t)+cos(2*pi*100*t);
% f = Fs*(0:(L/2))/L;   % frequency domain
% eg2
% wave=textread('elcentro.txt'); %地震波导入
% T=0.02;                         % Sampling period  
% X=wave(:,2);
% X=X/max(abs(X));
% L=length(X);          % Length of signal
% Fs = 1/T;             % Sampling frequency   
% t = (0:L-1)*T;        % Time vector
% f = Fs*(0:(L/2))/L;   % frequency domain
%% 进行傅里叶变化，并求单边幅值谱及单边相位谱
    X=X(:)'; %保证为行向量
    % 快速傅里叶变化
    Y = fft(X);
    % 单边幅值谱
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    % 单边相位谱
    ph=angle(Y(1:L/2+1));
%     ph(2:end-1)=ph(2:end-1);
%     ph=mod(ph+pi/2,pi)-pi/2;
%     ph=ph*180/pi;
   ph = rad2deg(ph);
    % ph(find(P1<1e-6))=0; %将幅值较小的相位谱除去
%% 还原信号
    X1=zeros(1,length(t));
    X2=zeros(1,length(t));
    n=0:L-1;
    for k=1:L/2+1  %注意k=1及k=L/2+1的频谱宽度为1/2L,其它频谱宽度为1/L
         if k==1||k==L/2+1
             VM=2;
         else
             VM=1;
         end
         X1=X1+(2*real(Y(k))*cos(2*pi*(k-1)/L*n)-2*imag(Y(k))*sin(2*pi*(k-1)/L*n))/(VM*L);
    %下面这种叠加方法无法测试成功
     %      XA(k)=sqrt((2*real(Y(k)))^2+(2*imag(Y(k)))^2)/(VM*L);
    %      Xphi(k)=atand(imag(Y(k))/real(Y(k)));
    %      X2=X2+XA(k)*cos(2*pi*(k-1)/L*n+Xphi(k)/180*pi);
    end
    iFFT_error=X-X1;  %逆向傅里叶误差
    
  
    
%% 绘图，采用参数PTO控制是否绘图
    if  PTO(1)==1
        % 原时域信号
        figure
        subplot(3,1,1)
        plot(t,X)
        title('Original signal X(t)')
        % 单边幅值谱
        subplot(3,1,2)
        plot(f,P1) 
        title('Single-Sided Amplitude Spectrum of X(t)')
        xlabel('f (Hz)')
        ylabel('|P1(f)|')
        % 单边相位谱
        subplot(3,1,3)
        plot(f,ph) 
        axis([0,Fs/2,-180,180])
        set(gca,'YTick',-180:90:180);
        title('Single-Sided Phase Spectrum of X(t)')
        xlabel('f(Hz)')
        ylabel('|Phase(°)|')
    end
    if PTO(2)==1
        figure
        plot(t,iFFT_error,'--')
        legend('信号逆向傅里叶变化误差')
    end
return