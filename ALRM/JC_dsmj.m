function [beta Xj] = JC_dsmj(mu,sigma,X0,TYPE,g_JC,Dg_JC)
%JC_method_RS:响应面法中调用的JC函数
%自变量信息
%初始化X0=Xi
X0=X0(:);  %保证输入向量为列向量
xi=X0;
time=0;
XX=[];Bbeta=[];AalphaX=[];
error=1;
%% 迭代计算
    while error>1e-6
        time=time+1; %计算次数统计
        normX=norm(xi);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%导数及偏导计算%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        g = g_JC(xi);
        Dg = Dg_JC(xi);           
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %计算非正态变量在验算正态当量化后的均值与方差
        [niuX sigmaX]=Xtran(mu,sigma,xi,TYPE);
        %迭代计算验算点位置
        gs=Dg.*sigmaX;
        alphaX=-gs/norm(gs);
        beta=(g+Dg'*(niuX-xi))/norm(gs);   
        xj=niuX+beta*sigmaX.*alphaX;    
        %中间量存储
        XX=[XX;xj'];
        Bbeta=[Bbeta;beta];
        AalphaX=[AalphaX;alphaX'];
        %计算误差，再交换变量循环运行
        error=norm(xj-xi)/norm(xi);
        xi=xj;
        
        % 强制暂停
        if time > 100
            break
        end
    end
    Xj=xj;
return