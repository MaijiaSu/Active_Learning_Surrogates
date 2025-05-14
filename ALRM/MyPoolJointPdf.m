function PoolPdf = MyPoolJointPdf(mu,sigma,TYPE,Pool)
%% ˵��
% 2020.04.28 (�����ƣ�Ŀǰ�൱��ֻ���˶�����̬�ֲ������)

% MyCdf�� �õ����ϸ����ܶȺ��� 

% input  ����������� mu,sigama 
%        ����������� TYPE
%        ������

% output ������������� Pool (��ģΪ NEP*Ndim)
%        ��������Ӧ�ĸ����ܶ� PoolPdf (��ģΪ NEP*1)

% ����������������������������������TYPE����������������������������������������
% TYPE = 1 : ��̫�ֲ�
% TYPE = 2 : ������̫�ֲ�
% TYPE = 3 : ����ֵI�ͷֲ�
% TYPE = 4 : ��СֵI�ͷֲ�
% TYPE = 5 : ���ȷֲ�
% TYPE =  'LHS_normal'    ���������������� - ��̬�ֲ�
% TYPE =  'LHS_uniform'   ���������������� - ���ȷֲ�
% TYPE =  'mvnrnd'        ����̬�ֲ������ϸ��ʷֲ�
% ������������������������������������������������������������������������������

% ������������������������������test example������������������������������������
% eg1
% niu=[0.36 0.18 20];     
% sigma=[0.036 0.018 5.0];
% TYPE=[1 2 4]; 
% NEP=1000;
% [Pool,PoolPdf] = RNgeneratorV2(mu,sigma,TYPE,NEP)
% eg2
% mu = [0.36 0.18];     
% sigma = [0.036 0.018];
% TYPE = ['LHS_uniform']; 
% NEP=10;
% [Pool,PoolPdf] = RNgeneratorV2(mu,sigma,TYPE,NEP)
% plot(Pool(:,1),Pool(:,2),'ko')
% ��������������������������������������������������������������������������������

% ��ʼ������
    Ndim=length(mu);   
    NEP = size(Pool,1);
    PoolPdf=zeros(NEP,Ndim);
%%  TYPE = num    

if ~ischar(TYPE)     %��TYPE�������ͷ��ַ�
    
    for i = 1:Ndim
        type = TYPE(i);
        switch type
            case 1   %��̬�������
                par = [mu(i),sigma(i)];
                PoolPdf(:,i) = normpdf(Pool(:,i),par(1),par(2));
                
            case 2   %������̬�������
                %����ת��
                deta = sigma(i)/mu(i);
                sLn = sqrt(log(1+deta^2));
                mLn = log(mu(i))-sLn^2/2;
                par = [mLn sLn];
                Pool(:,i) = lognrnd(par(1),par(2),NEP,1);
                PoolPdf(:,i) = lognpdf(Pool(:,i),par(1),par(2));
                
            case 3   %��ֵ1���������(����ֵ)
                %����ת��
                aEv = pi/sqrt(6)/sigma(i);
                uEv = psi(1)/aEv+mu(i); %-psi(1)Ϊŷ������
                par = [uEv aEv];
                Pool(:,i) = -1*evrnd(-par(1),1/par(2),NEP,1);
                PoolPdf(:,i) = evpdf(-Pool(:,i),-par(1),1/par(2));
                
            case 4    %��ֵ1���������(��Сֵ)
                %����ת��
                aEv = sqrt(6)*sigma(i)/pi;
                uEv = -psi(1)*aEv+mu(i); %-psi(1)Ϊŷ������
                par = [uEv aEv];
                Pool(:,i) = evrnd(par(1),par(2),NEP,1);
                PoolPdf(:,i) = evpdf(Pool(:,i),par(1),par(2));
                
            case 5  %���ȷֲ�
                a = mu(i)-sqrt(3)*sigma(i);
                b = mu(i)+sqrt(3)*sigma(i);
                Pool(:,i) = a+(b-a).*rand(NEP,1);
                PoolPdf(:,i) = 1/(b-a)*ones(NEP,1);
                
            otherwise
                error('type�������ͳ�����Χ')
        end
    end
    
    % �������ϸ����ܶ�
    PoolPdf = prod(PoolPdf')';
    
end