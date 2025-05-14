function [beta,Xj,iter,DoE] = JC_dsmj(niu,sigma,X0,TYPE,g_JC,Dg_JC)
%JC_method_RS:��Ӧ�淨�е��õ�JC����
%�Ա�����Ϣ
%��ʼ��X0=Xi
niu = niu(:);
sigma = sigma(:);
X0=X0(:);  %��֤��������Ϊ������
xi=X0;
iter=0;
XX=[];Bbeta=[];AalphaX=[];
error=1;
%% ��������
    while error>1e-3
        iter=iter+1; %�������ͳ��
        normX=norm(xi);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%������ƫ������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        g = g_JC(xi);
        Dg = Dg_JC(xi);           
        
        DoE.X(iter,:) = xi;
        DoE.Y(iter,1) = g;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %�������̬������������̬��������ľ�ֵ�뷽��
        [niuX sigmaX]=Xtran(niu,sigma,xi,TYPE);
        %�������������λ��
        gs=Dg.*sigmaX;
        alphaX=-gs/norm(gs);
        beta=(g+Dg'*(niuX-xi))/norm(gs);   
        xj=niuX+beta*sigmaX.*alphaX;    
        %�м����洢
        XX=[XX;xj'];
        Bbeta=[Bbeta;beta];
        AalphaX=[AalphaX;alphaX'];
        %�������ٽ�������ѭ������
        error=norm(xj-xi)/norm(xi);
        xi=xj;
    end
    Xj=xj;
return