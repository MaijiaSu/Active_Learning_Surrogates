function [Xi,beta,DoE] = RsUsingJC(muX,sigmaX,TYPE,f,error_tol,X0,fname,fun_par)
%% 
% �������ܣ�������ζ���ʽ��Ӧ�棬�����޷�ֱ�Ӽ����ݶȵĹ��ܺ�������Ƶ�
% ����ʽ��ʽѡȡΪ���ζ���ʽ�����Һ��Խ���˻������ϵ������Ϊ(2n+1)��
% ���������ʽ������JC�����������
% input�� �������ͳ���� muX��sigmaX��TYPE (����XTRAN����ת���Ǳ�׼��̬�������) 
%         ���ܺ���fname
%         ��ʼ̩��չ���� X0
%         �㷨����: ̩��չ��ϵ�� f �� ��ֹ������� error_tol 
% ouput�� ��Ƶ�����X���ɿ���ָ�� beta
%% ��ʼ����Ϣ
Xi=X0;
Xi=Xi(:);             %ǿ��Ϊ������
error = 1;
Time = 0;
X_path = [];          % ����ÿһ���������չ����
beta_path = [];       % �ɿ���ָ��
time = 0;
DoE.X =[];
DoE.Y = [];
%% ����ѭ������
while error >= error_tol
    X_path  = [X_path  Xi];

    % 1.ȷ����Ӧ�浱ǰ�������Doe_RS
    Ndim = length(muX); 
    DOE_RS = []; DOE_RS = Xi;
    for i = 1:Ndim
        Xtemp1 = Xi; Xtemp2 = Xi;
        
        Xtemp1(i) = Xtemp1(i) + f*sigmaX(i);
        Xtemp2(i) = Xtemp2(i) - f*sigmaX(i);
        
        DOE_RS = [DOE_RS,Xtemp1,Xtemp2];
    end
    
    % 2.���ù��ܺ������㺯��ֵ
    g = fname(DOE_RS',fun_par);      
    DoE.X = [DoE.X;DOE_RS'];
    DoE.Y = [DoE.Y;g];
    
    % 3.�����Ӧ�溯������ϵ��
    %����ϵ������ 
    A = []; A = ones(2*Ndim+1,1);
    A = [A,DOE_RS',DOE_RS'.^2];
    %��ⷽ��
    RS_PAR = A\g;
    
    % 4.����JC��������Ƶ�
    a=RS_PAR(1); b=RS_PAR(2:Ndim+1); c=RS_PAR(Ndim+2:2*Ndim+1);
    g_JC = @(x) a+b'*x+x'*diag(c)*x;   % %����ֵ
    Dg_JC = @(x) b+2*c.*x;             %������
    [beta Xi_DP] = JC_dsmj(muX,sigmaX,Xi,TYPE,g_JC,Dg_JC);
    beta_path = [beta_path beta];
    
   % 5�������,���ж��Ƿ���ֹ����
    error = norm(Xi_DP-Xi)/norm(Xi);
   
   % 6.���µ�����
    Xi = Xi_DP;   %��������
    

end