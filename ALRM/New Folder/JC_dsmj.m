function [beta Xj] = JC_dsmj(mu,sigma,X0,TYPE,g_JC,Dg_JC)
%JC_method_RS:��Ӧ�淨�е��õ�JC����
%�Ա�����Ϣ
%��ʼ��X0=Xi
X0=X0(:);  %��֤��������Ϊ������
xi=X0;
time=0;
XX=[];Bbeta=[];AalphaX=[];
error=1;
%% ��������
    while error>1e-6
        time=time+1; %�������ͳ��
        normX=norm(xi);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%������ƫ������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        g = g_JC(xi);
        Dg = Dg_JC(xi);           
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %�������̬������������̬��������ľ�ֵ�뷽��
        [niuX sigmaX]=Xtran(mu,sigma,xi,TYPE);
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
        
        % ǿ����ͣ
        if time > 100
            break
        end
    end
    Xj=xj;
return