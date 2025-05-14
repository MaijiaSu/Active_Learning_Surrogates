function [niuX sigmaX]=Xtran(niu,sigma,X0,TYPE)
%�������̫�ֲ������ĵ�����̬����ľ�ֵ�뷽��
%1.�Ƚ�niu��sigmaת��Ϊmatlab��Ӧ�ֲ�������Ĳ���
%2.�������̫�ֲ�������㴦���ۼƷֲ�������ܶ�
%3.���ݹ�ʽ���㵱����̫����ľ�ֵ�뷽��
    %��̫�ֲ�type=1
    %�����ֲ�type=2
    %��ֵ1�ͷֲ�type=3
    n=length(TYPE);
    for i=1:n
        type=TYPE(i);
        par=zeros(1,2);
        if type==1
            sigmaX=sigma;
            niuX=niu;
        elseif type==2
            %1����ת��
            deta=sigma(i)/niu(i);
            sLn=sqrt(log(1+deta^2));
            mLn=log(niu(i))-sLn^2/2;
            par=[mLn sLn];
            %2�������������ۼƷֲ�ֵ�������ܶȡ��溯��ֵ
            cdfX=logncdf(X0(i),par(1),par(2));
            pdfX=lognpdf(X0(i),par(1),par(2));
            invX=norminv(cdfX);
            %3���㵱������ı�׼�����ֵ
            sigmaX(i)=normpdf(invX)/pdfX;
            niuX(i)=X0(i)-invX*sigmaX(i);
        elseif type==3
            %1����ת��
            aEv=pi/sqrt(6)/sigma(i);
            uEv=psi(1)/aEv+niu(i); %-psi(1)Ϊŷ������
            par=[uEv aEv];
            % 5.8477   -46.6246
            %2�������������ۼƷֲ�ֵ�������ܶȡ��溯��ֵ
            cdfX=1-evcdf(-X0(i),-par(1),1/par(2));
            pdfX=evpdf(-X0(i),-par(1),1/par(2));
            invX=norminv(cdfX);
            %3���㵱������ı�׼�����ֵ
            sigmaX(i)=normpdf(invX)/pdfX;
            niuX(i)=X0(i)-invX*sigmaX(i);
        else
            disp('���������ֲ����ͣ�����')
        end
        %ǿ��ת��Ϊ������
        niuX = niuX(:);
        sigmaX = sigmaX(:);
    end