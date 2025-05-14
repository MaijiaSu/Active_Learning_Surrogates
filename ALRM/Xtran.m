function [niuX sigmaX]=Xtran(niu,sigma,X0,TYPE)
%计算非正太分布变量的当量正态化后的均值与方差
%1.先将niu与sigma转化为matlab相应分布需输入的参数
%2.计算非正太分布在验算点处的累计分布与概率密度
%3.根据公式计算当量正太化后的均值与方差
    %正太分布type=1
    %对数分布type=2
    %极值1型分布type=3
    n=length(TYPE);
    for i=1:n
        type=TYPE(i);
        par=zeros(1,2);
        if type==1
            sigmaX=sigma;
            niuX=niu;
        elseif type==2
            %1参数转化
            deta=sigma(i)/niu(i);
            sLn=sqrt(log(1+deta^2));
            mLn=log(niu(i))-sLn^2/2;
            par=[mLn sLn];
            %2计算在验算点的累计分布值、概率密度、逆函数值
            cdfX=logncdf(X0(i),par(1),par(2));
            pdfX=lognpdf(X0(i),par(1),par(2));
            invX=norminv(cdfX);
            %3计算当量化后的标准差与均值
            sigmaX(i)=normpdf(invX)/pdfX;
            niuX(i)=X0(i)-invX*sigmaX(i);
        elseif type==3
            %1参数转化
            aEv=pi/sqrt(6)/sigma(i);
            uEv=psi(1)/aEv+niu(i); %-psi(1)为欧拉常数
            par=[uEv aEv];
            % 5.8477   -46.6246
            %2计算在验算点的累计分布值、概率密度、逆函数值
            cdfX=1-evcdf(-X0(i),-par(1),1/par(2));
            pdfX=evpdf(-X0(i),-par(1),1/par(2));
            invX=norminv(cdfX);
            %3计算当量化后的标准差与均值
            sigmaX(i)=normpdf(invX)/pdfX;
            niuX(i)=X0(i)-invX*sigmaX(i);
        else
            disp('其它变量分布类型，暂无')
        end
        %强制转化为列向量
        niuX = niuX(:);
        sigmaX = sigmaX(:);
    end