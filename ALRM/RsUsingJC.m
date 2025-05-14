function [Xi,beta,DoE] = RsUsingJC(muX,sigmaX,TYPE,f,error_tol,X0,fname,fun_par)
%% 
% 函数功能：构造二次多项式响应面，计算无法直接计算梯度的功能函数的设计点
% 多项式形式选取为二次多项式，并且忽略交叉乘积项，待定系数个数为(2n+1)个
% 构造完多项式，采用JC法计算验算点
% input： 随机变量统计量 muX，sigmaX，TYPE (采用XTRAN函数转化非标准正态随机变量) 
%         功能函数fname
%         初始泰勒展开点 X0
%         算法参数: 泰勒展开系数 f ； 终止迭代误差 error_tol 
% ouput： 设计点坐标X，可靠度指标 beta
%% 初始化信息
Xi=X0;
Xi=Xi(:);             %强制为行向量
error = 1;
Time = 0;
X_path = [];          % 储存每一步计算出的展开点
beta_path = [];       % 可靠度指标
time = 0;
DoE.X =[];
DoE.Y = [];
%% 迭代循环计算
while error >= error_tol
    X_path  = [X_path  Xi];

    % 1.确定响应面当前的试验点Doe_RS
    Ndim = length(muX); 
    DOE_RS = []; DOE_RS = Xi;
    for i = 1:Ndim
        Xtemp1 = Xi; Xtemp2 = Xi;
        
        Xtemp1(i) = Xtemp1(i) + f*sigmaX(i);
        Xtemp2(i) = Xtemp2(i) - f*sigmaX(i);
        
        DOE_RS = [DOE_RS,Xtemp1,Xtemp2];
    end
    
    % 2.调用功能函数计算函数值
    g = fname(DOE_RS',fun_par);      
    DoE.X = [DoE.X;DOE_RS'];
    DoE.Y = [DoE.Y;g];
    
    % 3.求解响应面函数待定系数
    %方程系数矩阵 
    A = []; A = ones(2*Ndim+1,1);
    A = [A,DOE_RS',DOE_RS'.^2];
    %求解方程
    RS_PAR = A\g;
    
    % 4.采用JC法计算设计点
    a=RS_PAR(1); b=RS_PAR(2:Ndim+1); c=RS_PAR(Ndim+2:2*Ndim+1);
    g_JC = @(x) a+b'*x+x'*diag(c)*x;   % %标量值
    Dg_JC = @(x) b+2*c.*x;             %列向量
    [beta Xi_DP] = JC_dsmj(muX,sigmaX,Xi,TYPE,g_JC,Dg_JC);
    beta_path = [beta_path beta];
    
   % 5计算误差,并判断是否终止迭代
    error = norm(Xi_DP-Xi)/norm(Xi);
   
   % 6.更新迭代点
    Xi = Xi_DP;   %交换变量
    

end