function [P1 ph iFFT_error f]=FFT_dsmj(X,T,PTO)
%�������ܣ�����FFT�õ������źŵ�Ƶ����
%input  ʱ���ź�X���źŲ�������T,��ͼ����
%output ���߷�ֵ�׼�������λ�ף��Լ�������Ҷ�仯���
%���� 2019/5/21 dsmj
%% ��ʼ��
% % input
% Original signal X(t)
% T=0.02;                % Sampling period  
% ��ʼ��
    L=length(X);          % Length of signal
    Fs = 1/T;             % Sampling frequency   
    t = (0:L-1)*T;        % Time vector
    f = Fs*(0:(L/2))/L;   % frequency domain
%% ����
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
% wave=textread('elcentro.txt'); %���𲨵���
% T=0.02;                         % Sampling period  
% X=wave(:,2);
% X=X/max(abs(X));
% L=length(X);          % Length of signal
% Fs = 1/T;             % Sampling frequency   
% t = (0:L-1)*T;        % Time vector
% f = Fs*(0:(L/2))/L;   % frequency domain
%% ���и���Ҷ�仯�����󵥱߷�ֵ�׼�������λ��
    X=X(:)'; %��֤Ϊ������
    % ���ٸ���Ҷ�仯
    Y = fft(X);
    % ���߷�ֵ��
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    % ������λ��
    ph=angle(Y(1:L/2+1));
%     ph(2:end-1)=ph(2:end-1);
%     ph=mod(ph+pi/2,pi)-pi/2;
%     ph=ph*180/pi;
   ph = rad2deg(ph);
    % ph(find(P1<1e-6))=0; %����ֵ��С����λ�׳�ȥ
%% ��ԭ�ź�
    X1=zeros(1,length(t));
    X2=zeros(1,length(t));
    n=0:L-1;
    for k=1:L/2+1  %ע��k=1��k=L/2+1��Ƶ�׿��Ϊ1/2L,����Ƶ�׿��Ϊ1/L
         if k==1||k==L/2+1
             VM=2;
         else
             VM=1;
         end
         X1=X1+(2*real(Y(k))*cos(2*pi*(k-1)/L*n)-2*imag(Y(k))*sin(2*pi*(k-1)/L*n))/(VM*L);
    %�������ֵ��ӷ����޷����Գɹ�
     %      XA(k)=sqrt((2*real(Y(k)))^2+(2*imag(Y(k)))^2)/(VM*L);
    %      Xphi(k)=atand(imag(Y(k))/real(Y(k)));
    %      X2=X2+XA(k)*cos(2*pi*(k-1)/L*n+Xphi(k)/180*pi);
    end
    iFFT_error=X-X1;  %������Ҷ���
    
  
    
%% ��ͼ�����ò���PTO�����Ƿ��ͼ
    if  PTO(1)==1
        % ԭʱ���ź�
        figure
        subplot(3,1,1)
        plot(t,X)
        title('Original signal X(t)')
        % ���߷�ֵ��
        subplot(3,1,2)
        plot(f,P1) 
        title('Single-Sided Amplitude Spectrum of X(t)')
        xlabel('f (Hz)')
        ylabel('|P1(f)|')
        % ������λ��
        subplot(3,1,3)
        plot(f,ph) 
        axis([0,Fs/2,-180,180])
        set(gca,'YTick',-180:90:180);
        title('Single-Sided Phase Spectrum of X(t)')
        xlabel('f(Hz)')
        ylabel('|Phase(��)|')
    end
    if PTO(2)==1
        figure
        plot(t,iFFT_error,'--')
        legend('�ź�������Ҷ�仯���')
    end
return