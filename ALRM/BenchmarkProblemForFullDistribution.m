%%
% clc,clear
% ExampleOrder = [10	2	12	15	14	13	16	1	4	8	9	3	5	7	19	6];
% NN = 6;
% ExNo = ExampleOrder(NN);
% TestExample =  ['BatchProcess','Eg',num2str(ExNo)];

ProSys = [];
if strcmp(TestExample,'BatchProcessEg1')
    % RP-14
    % Defination of  Random varibales
    ProSys.muX = [75.0,39.0,1500.0,400.0,250000.0];
    ProSys.sigmaX = [2.887,0.1,350.0,0.1,35000.0];
    ProSys.Distri = [5,1,3,1,1];
    % Define Random model according to Parameters
    ProSys.RVPar.P1 = [70.0,39.0,1342.0,400.0,25e4];
    ProSys.RVPar.P2 = [80.0,0.1,272.9,0.1,35e3];
    ProSys.RVPar.isInputPar = [0,0,0,0,0];  % =1,'off'
    % Perfomance function
    ProSys.fname= @(X,P) X(:,1)-32/pi./X(:,2).^3.*...
        sqrt(X(:,3).^2.*X(:,4).^2/16+X(:,5).^2);
    ProSys.fun_par = [];
    ProSys.Ndim = 5;
    % Paramenter for UQ
    ProSys.ymin = 0.9154;
    ProSys.ymax = 44.2573;

elseif strcmp(TestExample,'BatchProcessEg2')
    % RP-53
    % Defination of  Random varibales
    ProSys.muX = [1.5,2.5];
    ProSys.sigmaX = [1.0,1.0];
    ProSys.Distri = [1,1];
    % Define Random model according to Parameters
    % ...
    % Perfomance function
    ProSys.fname= @(X,P) sin(5*X(:,1)/2)+2-(X(:,1).^2+4).*(X(:,2)-1)/20;
    ProSys.fun_par = [];
    ProSys.Ndim = 2;
    % Paramenter for UQ
    ProSys.ymin = -1.7079;
    ProSys.ymax = 3.4104;

elseif strcmp(TestExample,'BatchProcessEg3')
    % RP-38
    % Defination of  Random varibales
    ProSys.muX = [350,50.8,3.81,173,9.38,33.1,0.036];
    ProSys.sigmaX = [35,5.08,0.381,17.3,0.938,3.31,0.0036];
    ProSys.Distri = [1,1,1,1,1,1,1];
    % Define Random model according to Parameters
    % ...
    % Perfomance function
    ProSys.fname= @(X,P) 15.59*1e4-X(:,1).*X(:,2).^3/2./X(:,3).^3.*...
        ( X(:,4).^2-4*X(:,5).*X(:,6).*X(:,7).^2+X(:,4).*(X(:,6)+4*X(:,5)+2*X(:,6).*X(:,7)) )./...
        (X(:,4).*X(:,5).*(X(:,4)+X(:,6)+2*X(:,6).*X(:,7)));
    ProSys.fun_par = [];
    ProSys.Ndim = 7;
    % Paramenter for UQ
    ProSys.ymin = -6.1008e+04;
    ProSys.ymax = 1.4331e+05;

elseif strcmp(TestExample,'BatchProcessEg4')
    % Dynamic reliablity of nonlinear occilator 
    % Defination of  Random varibales
    ProSys.muX = [1 1 0.1 0.5 1 1];
    ProSys.sigmaX = [0.05 0.1 0.01 0.05 0.2 0.2];
    ProSys.Ndim = 6;
    ProSys.Distri = ones(1,ProSys.Ndim);
    % Define Random model according to Parameters
    % ...
    % Perfomance function
    ProSys.fname = @DRONLO;
    ProSys.fun_par.StadNrom = 'No';
    % Paramenter for UQ
    ProSys.ymin = -0.4273;
    ProSys.ymax = 1.4065;

elseif strcmp(TestExample,'BatchProcessEg5')
     % Defination of  Random varibales
     ProSys.Ndim = 8;
     ProSys.muX = [ 0.1	3698.252464	89335	1045	89.55	760	1400	10950];
     ProSys.sigmaX = [0.028867513	4890.907662	15164.10482	31.75426481	15.27091462	34.64101615	161.6580754	632.1985448];
     ProSys.Distri = [5,2,5,5,5,5,5,5];
     % Define Random model according to Parameters
     ProSys.RVPar.P1 = [0.05,7.71,63070,990,63.1,700,1120,9855];
     ProSys.RVPar.P2 = [0.15,1.0056,115600,1100,116,820,1680,12045];
     ProSys.RVPar.isInputPar = [1,1,1,1,1,1,1,1];  % =1,'on'
     % Perfomance function
     ProSys.fname = @BoreholeFun;
     ProSys.fun_par = [];
     % Paramenter for UQ
     ProSys.ymin =  11.9112;
     ProSys.ymax = 230.0507;

elseif strcmp(TestExample,'BatchProcessEg6')
    % Defination of  Random varibales
    ProSys.Ndim = 21;
    ProSys.muX = [1.33E+05	8.90E+04  7.12E+04	2.1738E+10	2.3796E+10	8.1344E-03...
        1.5090E-02	2.1375E-02	2.5961E-02	1.0812E-02	1.4105E-02	2.3279E-02  ...
        2.5961E-02	3.1256E-01	3.7210E-01	5.0606E-01	5.5815E-01	2.5302E-01 ...
        2.9117E-01	3.7303E-01	4.1860E-01];
    ProSys.sigmaX = [4.04E+04	3.56E+04	2.85E+04	1.9152E+09	1.9152E+09	1.0834E-03...
        1.2980E-03	2.5961E-03	3.0288E-03	2.5961E-03	3.4615E-03	5.6249E-03...
        6.4902E-03	5.5815E-03	7.4420E-03	9.3025E-02	1.1163E-01	9.3025E-02...
        1.0232E-01	1.2093E-01	1.9537E-01];
    ProSys.Distri = [ones(1,3)*2,ones(1,18)*6];
    % Perfomance function
    ProSys.fname = @FrameStru_3b5F;
    ProSys.fun_par = [];
    % Paramenter for UQ (10^6 MCS samples, Quantile of 1e-3 and 1-1e-3)
    ProSys.ymin =  6.8862e-03;
    ProSys.ymax =  4.1056e-02;

elseif strcmp(TestExample,'BatchProcessEg7')
    % Defination of  Random varibales
    ProSys.Ndim = 10;
    ProSys.muX = [2.0e-3, 2.1e11, 1.0e-3, 2.1e11, 5.0e4, 5.0e4, 5.0e4, 5.0e4, 5.0e4, 5.0e4];
    ProSys.sigmaX  = [2.0e-4, 2.1e10, 1.0e-4, 2.1e10, 7.5e3, 7.5e3, 7.5e3, 7.5e3, 7.5e3, 7.5e3];
    ProSys.Distri = [2,2,2,2,3,3,3,3,3,3];
    % Perfomance function
    ProSys.fun_par.DispThreshold = 12; %unit:cm
    ProSys.fname = @main_23BarsTruss;
    % Paramenter for UQ
    ProSys.ymin = -2.1908e-01;
    ProSys.ymax = 6.8021e+00;

elseif strcmp(TestExample,'BatchProcessEg8')
    % Defination of  Random varibales
    ProSys.Ndim = 6;
    ProSys.muX = [7.0e10,2.5e-3,0.524,0.90,8.0e4,7.0e4];
    ProSys.sigmaX  = [3.50e9,1.25e-4,0.010480,0.0225,6.4e3,5.60e3];
    ProSys.Distri = [1,1,1,1,1,1];
    % Define Random model according to Parameters
    %
    % Perfomance function
    ProSys.fun_par = [];
    ProSys.fname= @(X,P) 1-...
                sqrt(3*(1-0.3^2))./(pi*X(:,1).*X(:,2).^2.*cos(X(:,3).^2)).*...
                   (X(:,6)/0.66+X(:,5)./X(:,4)/0.41);    
    % Paramenter for UQ
    ProSys.ymin = 3.8951e-01;
    ProSys.ymax = 7.2532e-01;

elseif strcmp(TestExample,'BatchProcessEg9')
    % Defination of  Random varibales
    ProSys.Ndim = 6;
    ProSys.muX = [0.9377,2.2e5,2.1e4,0.29,24.0,8];
    ProSys.sigmaX  = [0.0459,5.0e3,1e3,5.8e-3,0.5,0.3];
    ProSys.Distri = [7,1,1,5,1,1];
    % Define Random model according to Parameters
    ProSys.RVPar.P1 = [0.958,2.2e5,2.1e4,0.28,24.0,8];
    ProSys.RVPar.P2 = [22.508,5.0e3,1e3,0.3,0.5,0.3];
    ProSys.RVPar.isInputPar = [1,0,0,1,0,0];  % =1,'off'
    % Perfomance function
    ProSys.fun_par = [];
    ProSys.fname= @(X,P) sqrt(X(:,1).*X(:,2)*3*385.82.*(X(:,5)-X(:,6)))./...
                          sqrt(X(:,4).*(X(:,3)*2*pi/60).^2.*(X(:,5).^3-X(:,6).^3))-0.37473;
    % Paramenter for UQ
    ProSys.ymin = -4.2136e-03;
    ProSys.ymax =  1.7049e-01;

elseif strcmp(TestExample,'BatchProcessEg10')
    % Defination of  Random varibales
    ProSys.Ndim = 2;
    ProSys.muX = [1,3];
    ProSys.sigmaX  = [1,3];
    ProSys.Distri = [1,1];
    % Define Random model according to Parameters
    % ...
    % Perfomance function
    ProSys.fun_par = [];
    ProSys.fname= @(X,P) 2.5-0.2357*(X(:,1)-X(:,2))+0.00463*(X(:,1)+X(:,2)-20).^4;
    % Paramenter for UQ
    ProSys.ymin = 11.8043;
    ProSys.ymax = 2.0610e+03;

elseif strcmp(TestExample,'BatchProcessEg11')  
    % Defination of  Random varibales
    ProSys.Ndim = 2;
    ProSys.muX = [78064.0,0.0104];
    ProSys.sigmaX  = [11710.0,0.00156];
    ProSys.Distri = [1,1];
    % Define Random model according to Parameters
    % ...
    % Perfomance function
    ProSys.fun_par = [];
    ProSys.fname= @(X,P) X(:,1).*X(:,2)-146.14;
    % Paramenter for UQ
    ProSys.ymin = 207.1328;
    ProSys.ymax = 1.2743e+03;

elseif  strcmp(TestExample,'BatchProcessEg12')     
    % Defination of  Random varibales
    ProSys.Ndim = 2;
    ProSys.muX = [0,0];
    ProSys.sigmaX = [1,1];
    ProSys.Distri = [1,1];
    % Define Random model according to Parameters
    % ...
    % Perfomance function
    ProSys.fname=@FBSS;
    ProSys.fun_par.type=1;ProSys.fun_par.a = 3;ProSys.fun_par.b = 7;
    ProSys.fun_par.f0 = 0;
    % Paramenter for UQ
    ProSys.ymin = -0.2303;
    ProSys.ymax =  3.2396;

elseif strcmp(TestExample,'BatchProcessEg13')   
    % Defination of  Random varibales
    ProSys.Ndim = 2;
%     ProSys.muX = [0,0];
%     ProSys.sigmaX = [1,1];
%     ProSys.Distri = [1,1];
    ProSys.muX = [0,0];
    po = 0.95;
    ProSys.sigmaX = [1,po;po,1];
    ProSys.Distri = 'mvnrnd';
    % Define Random model according to Parameters
    % ...
    % Perfomance function
    f0 = 0;
    ProSys.fname = @(x,fun_par) 10-(x(:,1).^2-5*cos(2*pi*x(:,1)))-...
        (x(:,2).^2-5*cos(2*pi*x(:,2)))-f0;
    ProSys.fun_par = [];
    % Paramenter for UQ
%     ProSys.ymin = -8.4074;
%     ProSys.ymax =  19.8020;
    ProSys.ymin = -12.5791;
    ProSys.ymax =  19.9331;

elseif strcmp(TestExample,'BatchProcessEg14')  
    % Defination of  Random varibales
    ProSys.Ndim = 2;
    ProSys.muX = [0,0];
    ProSys.sigmaX = [1,1];
    ProSys.Distri = [1,1];
    % Define Random model according to Parameters
    % ...
    % Perfomance function
    ProSys.fun_par.c = 3;
    ProSys.fname = @(x,fun_par) min([fun_par.c-1-x(:,2)+exp(-x(:,1).^2/10)+(x(:,1)/5).^4,...
        fun_par.c^2/2 - x(:,1).*x(:,2)]')';
    % Paramenter for UQ
    ProSys.ymin = -0.6955;
    ProSys.ymax =  5.6444;

elseif strcmp(TestExample,'BatchProcessEg15') 
    % Defination of  Random varibales
    ProSys.Ndim = 2;
    ProSys.muX = [0,0];
    ProSys.sigmaX = [1,1];
    ProSys.Distri = [1,1];
    % Define Random model according to Parameters
    % ...
    % Perfomance function
    ProSys.fun_par = [];
    ProSys.fname = @(X,fun_par) min([X(:,1)-X(:,2),X(:,1)+X(:,2)]')';
    % Paramenter for UQ
    ProSys.ymin =  -4.6649;
    ProSys.ymax = 2.6247;

elseif strcmp(TestExample,'BatchProcessEg16') 
    % Defination of  Random varibales
    ProSys.Ndim = 3;
    ProSys.a=[-pi,-pi,-pi];ProSys.b=[pi,pi,pi];
    ProSys.muX=[0,0,0];
    ProSys.sigmaX=[1.8138,1.8138,1.8138];
    ProSys.Distri=[5,5,5];
    % Define Random model according to Parameters
    ProSys.RVPar.P1 = [-pi,-pi,-pi];
    ProSys.RVPar.P2 = [pi,pi,pi];
    ProSys.RVPar.isInputPar = [1,1,1];  % =1,'on'
    % Perfomance function
    ProSys.fun_par.a = 7;
    ProSys.fun_par.b = 0.1;
    ProSys.fname = @(X,fun_par) sin(X(:,1))+fun_par.a*sin(X(:,2)).^2+...
        fun_par.b*X(:,3).^4.*sin(X(:,1));
    % Paramenter for UQ
    ProSys.ymin = -9.0611;
    ProSys.ymax = 16.1035;

elseif strcmp(TestExample,'BatchProcessEg17') 
    % Defination of  Random varibales
    ProSys.Ndim = 10;
    ProSys.muX = zeros(1,10);
    ProSys.sigmaX = ones(1,10);
    ProSys.Distri = ones(1,10);
    % Define Random model according to Parameters
    % ...
    % Perfomance function
    ProSys.fun_par = [];
    ProSys.fname = @(X,fun_par) min([X(:,1)-X(:,2),X(:,1)+X(:,2)]')';
    % Paramenter for UQ
    ProSys.ymin =  -4.6649;
    ProSys.ymax = 2.6247;

elseif strcmp(TestExample,'BatchProcessEg18') %
    % Defination of  Random varibales
    ProSys.Ndim = 2;
    ProSys.muX = exp(0.5)*ones(1,2);
    ProSys.sigmaX = exp(0.5)*sqrt(exp(1)-1)*ones(1,2);
    ProSys.Distri = 2*ones(1,2);
    % Define Random model according to Parameters
    ProSys.RVPar.P1 = [0,0];
    ProSys.RVPar.P2 = [0.5,0.5];
    ProSys.RVPar.isInputPar = [1,1];  % =1,'on'
    % Perfomance function
    ProSys.fun_par = [];
    %     ProSys.fname = @(X,fun_par) X(:,1).^2+X(:,1).*X(:,2)+X(:,2);
    % ProSys.fname = @(X,fun_par) min([X(:,1)-X(:,2),X(:,1)+X(:,2)]')';
    ProSys.fname = @(X,fun_par) X(:,1)+X(:,2);
    % Paramenter for UQ
    ProSys.ymin =  0.2540;
    ProSys.ymax =  28.6746;
 
elseif strcmp(TestExample,'BatchProcessEg19')  % 10F-frame-structural 
    % Defination of  Random varibales
    ProSys.fun_par.stru_type = 2; %=1,5F; =2,10F,=3,2F;
    mm = [2.82,2.76*ones(1,18),2.92]*1e5;
    kk = [1.007,1.358,1.338,1.33,1.228,1.248,1.237,1.223,1.208,1.187,1.093,...
        1.054,1.036,0.908,0.897,0.882,0.798,0.734,0.598,0.449]*1.0e8;
    ProSys.Ndim = 20;
    ProSys.muX = [mm(1:ProSys.Ndim/2),kk(1:ProSys.Ndim/2)];
    ProSys.sigmaX = [mm(1:10)*0.02,kk(1:10)*0.05];
    ProSys.Distri = ones(1,ProSys.Ndim)*2;

    % Define Random model according to Parameters
    % ...
    % Perfomance function
    ProSys.fname = @StruDynamicNF_nRV;
    ProSys.fun_par.indexType = 2; % story dirif
    ProSys.fun_par.stru_type = 2;
    ProSys.fun_par.StadNrom = 'No';
    ProSys.fun_par.threshold = 0;
    % Paramenter for UQ
    ProSys.ymin =   -0.0174;
    ProSys.ymax = -0.0121;

elseif strcmp(TestExample,'BatchProcessEg20')  % test 
    % Adaption of BatchProcessEg2 (from 2-dim to 4-dim)
    % Defination of  Random varibales
    ProSys.muX = [0.75,0.75,1.25,1.25];
    ProSys.sigmaX = [sqrt(2)/2,sqrt(2)/2,sqrt(2)/2,sqrt(2)/2];
    ProSys.Distri = [1,1,1,1];
    % Define Random model according to Parameters
    % ...
    % Perfomance function
    ProSys.fname= @(X,P) sin(5*(X(:,1)+X(:,2))/2)+2-((X(:,1)+X(:,2)).^2+4).*((X(:,3)+X(:,4))-1)/20;
    ProSys.fun_par = [];
    ProSys.Ndim = 4;
    % Paramenter for UQ
    ProSys.ymin = -1.7079;
    ProSys.ymax = 3.4104;

elseif strcmp(TestExample,'BatchProcessEg21')  % test 
    % Modified EASOM FUNCTION
    % Defination of  Random varibales
    ProSys.muX = [0,0];
    ProSys.sigmaX = [1,1];
    ProSys.Distri = [1,1];
    % Define Random model according to Parameters
    % ...
    % Perfomance function
    ProSys.fname= @(X,P) cos(10*X(:,1)).*cos(10*X(:,2)).*exp(10*(-(X(:,1)).^2-(X(:,2)).^2));
    ProSys.fun_par = [];
    ProSys.Ndim = 2;
    % Paramenter for UQ
    ProSys.ymin = -0.4213;
    ProSys.ymax = 0.8846;

elseif strcmp(TestExample,'BatchProcessEg22')  % test 
    ProSys.muX = [-5,-5];
    ProSys.sigmaX = [1,1];
    ProSys.Distri = [1,1];
    % Define Random model according to Parameters
    % ...
    % Perfomance function
    ProSys.fname= @(X,P) sin(X(:,2)).*exp((1-cos(X(:,1))).^2)+...
                            cos(X(:,1)).*exp((1-sin(X(:,2))).^2)+(X(:,1)-X(:,2)).^2;
    ProSys.fun_par = [];
    ProSys.Ndim = 2;
    % Paramenter for UQ
    ProSys.ymin = -103.2580;
    ProSys.ymax = 201.9534;

elseif strcmp(TestExample,'BatchProcessEg23')  % test
    % RP-53
    % Defination of  Random varibales
    ProSys.muX = [1.5,2.5];
    ProSys.sigmaX = [1.0,1.0];
    ProSys.Distri = [1,1];
    % Define Random model according to Parameters
    % ...
    % Perfomance function
    ProSys.fname= @(X,P) sin(5*X(:,1)/2)+2-(X(:,1).^2+4).*(X(:,2)-1)/20+...
                         5*exp(-pdist2(X,[1,2],'euclidean'))+...
                         6*exp(-pdist2(X,[3,1],'euclidean'))+...
                         4*exp(-pdist2(X,[-1,4],'euclidean'))+...
                         -4.5*exp(-pdist2(X,[2,5],'euclidean'))+...
                          -5*exp(-pdist2(X,[-1,1],'euclidean'));
%         ProSys.fname= @(X,P) 0+...
%                          5*exp(-pdist2(X,[0,0],'euclidean'));
    ProSys.fun_par = [];
    ProSys.Ndim = 2;
    % Paramenter for UQ
    ProSys.ymin = -1.7079;
    ProSys.ymax = 3.4104;

elseif strcmp(TestExample,'BatchProcessEg24')  % test
    % RP-53
    % Defination of  Random varibales
    ProSys.Ndim = 5;
    ProSys.muX =  zeros(1,ProSys.Ndim);
    ProSys.sigmaX = ones(1,ProSys.Ndim);
    ProSys.Distri = ones(1,ProSys.Ndim);
    % Define Random model according to Parameters
    % ...
    % Perfomance function
    %     ProSys.fname= @testfunction;
    rng(1);
    CenterP =  unifrnd(-3,3,10,5);
    ScaleFactor = unifrnd(0.5,1,10,1).*transpose(sum(transpose(CenterP)));
    rng('shuffle')
    ProSys.fname= @(X,P) transpose(sum(transpose(X))) +...
        ScaleFactor(1)*exp(-pdist2(X,CenterP(1,:),'euclidean'))+...
        ScaleFactor(2)*exp(-pdist2(X,CenterP(2,:),'euclidean'))+...
        ScaleFactor(3)*exp(-pdist2(X,CenterP(3,:),'euclidean'))+...
        ScaleFactor(4)*exp(-pdist2(X,CenterP(4,:),'euclidean'))+...
        ScaleFactor(5)*exp(-pdist2(X,CenterP(5,:),'euclidean'))+...
        ScaleFactor(6)*exp(-pdist2(X,CenterP(6,:),'euclidean'))+...
        ScaleFactor(7)*exp(-pdist2(X,CenterP(7,:),'euclidean'))+...
        ScaleFactor(8)*exp(-pdist2(X,CenterP(8,:),'euclidean'))+...
        ScaleFactor(9)*exp(-pdist2(X,CenterP(9,:),'euclidean'))+...
        ScaleFactor(10)*exp(-pdist2(X,CenterP(10,:),'euclidean'));


    ProSys.fun_par = [];

    % Paramenter for UQ
    ProSys.ymin = 0.6997;
    ProSys.ymax = 20.4919;

elseif strcmp(TestExample,'BatchProcessEg25')  % test
    % RP-53
    % Defination of  Random varibales
    ProSys.Ndim = 4;
     ProSys.muX = [1.5,2.5,zeros(1,2)];
    ProSys.sigmaX = [1.0,1.0,ones(1,2)];
    ProSys.Distri = ones(1,4);
    % Define Random model according to Parameters
    % ...
    % Perfomance function
    ProSys.fname= @(X,P) sin(5*X(:,1)/2)+2-(X(:,1).^2+4).*(X(:,2)-1)/20;
    ProSys.fun_par = [];
    % Paramenter for UQ
    ProSys.ymin = -1.7079;
    ProSys.ymax = 3.4104;

elseif strcmp(TestExample,'BatchProcessEg26')  % test
    % Defination of  Random varibales
    ProSys.Ndim = 2;
    ProSys.muX = [0,0];
    ProSys.sigmaX = [1,1];
    ProSys.Distri = ones(1,2);
    % Define Random model according to Parameters
    % ...
    % Perfomance function
    CenterP =  [0,0;
        -1,-1;
        -1,1;
        1,-1;
        1,1];
%     ScaleFactor = [1,2,2,1,1]*5;
ScaleFactor = [1,0,0,0,0]*5;
    ProSys.fname= @(X,P) X(:,1)-X(:,2)+...
        ScaleFactor(1)*exp(-pdist2(X,CenterP(1,:)./1,'euclidean'))+...
        ScaleFactor(2)*exp(-pdist2(X,CenterP(2,:)./1,'euclidean'))+...
        ScaleFactor(3)*exp(-pdist2(X,CenterP(3,:)./1,'euclidean'))+...
        ScaleFactor(4)*exp(-pdist2(X,CenterP(4,:)./1,'euclidean'))+...
        ScaleFactor(5)*exp(-pdist2(X,CenterP(5,:)./1,'euclidean'));
    ProSys.fun_par = [];

        % Paramenter for UQ
    ProSys.ymin = -4.1800;
    ProSys.ymax = 4.8518;

end
 

%%


%%
% check the mean and sigama
% mean_est = mean(SBM.SamplePool);
% std_est = std(SBM.SamplePool);
% (mean_est-ProSys.muX)./ProSys.muX
% (std_est-ProSys.sigmaX)./ProSys.sigmaX

%%
if 0
close all
NofSamples = 1e5;
 [SBM.SamplePool,~] = RNgeneratorV4(ProSys,NofSamples);
G = ProSys.fname(SBM.SamplePool,ProSys.fun_par);
Pf = numel(find(G<=0))/NofSamples

G_sort = sort(G,'ascend');
ymin =  G_sort(floor(NofSamples.*1e-3))
ymax =  G_sort(floor(NofSamples.*(1-1e-3)))
histogram(G,100,'Normalization','cdf');
[Nc,Xc] = histcounts(G,100);
XC = (Xc(1:end-1)+Xc(2:end))/2;
plot(XC,Nc)
y99 =  G_sort(floor(NofSamples.*(0.99)))
y999 =  G_sort(floor(NofSamples.*(0.999)))
y9999 =  G_sort(floor(NofSamples.*(0.9999)))
end
%%
% clc,clear,close all

% TestExample = 'eg12';  
% LimtStateFunction_select
% TestExample =  ['BatchProcess','Eg',num2str(3)];
% BenchmarkProblemForFullDistribution

    % Defination of  Random varibales
%     ProSys.Ndim = 2;
%      ProSys.muX = [0,0];
%     ProSys.sigmaX = [1,1];
%     ProSys.Distri = ones(1,2);
%     % Define Random model according to Parameters
%     % ...
%     % Perfomance function
%     CenterP =  [0,0;
%                 -1,-1;
%                 -1,1;
%                 1,-1;
%                 1,1];
%     ScaleFactor = [1,2,2,1,1]*5;
%     ProSys.fname= @(X,P) X(:,1)-X(:,2)+...
%                  ScaleFactor(1)*exp(-pdist2(X,CenterP(1,:)./1,'euclidean'))+...
%                  ScaleFactor(2)*exp(-pdist2(X,CenterP(2,:)./1,'euclidean'))+...
%                  ScaleFactor(3)*exp(-pdist2(X,CenterP(3,:)./1,'euclidean'))+...
%                  ScaleFactor(4)*exp(-pdist2(X,CenterP(4,:)./1,'euclidean'))+...
%                  ScaleFactor(5)*exp(-pdist2(X,CenterP(5,:)./1,'euclidean'));
%     ProSys.fun_par = [];
%     % Paramenter for UQ
%     ProSys.ymin = -1.7079;
%     ProSys.ymax = 3.4104;
% ----------------------------------------------------------------------
% figure
% vm = 5;
% for ii = 1:ProSys.Ndim  
%     if ProSys.Distri(ii) == 5
%         bound(ii,:) = [ProSys.RVPar.P1(ii),ProSys.RVPar.P2(ii)];
%     elseif ProSys.Distri(ii) == 6 || ProSys.Distri(ii) == 2
%         bound(ii,:) = [0,ProSys.muX(ii)+vm*ProSys.sigmaX(ii)];
%     else
%         bound(ii,:) = [ProSys.muX(ii)-vm*ProSys.sigmaX(ii),ProSys.muX(ii)+vm*ProSys.sigmaX(ii)];
%     end
% end
% % bound = [ProSys.muX-vm*ProSys.sigmaX;ProSys.muX+vm*ProSys.sigmaX]'
% density = 100;
% Myfun = @(X) ProSys.fname(X,ProSys.fun_par)
% % data = PlotHighDimSurface(Myfun,bound,density)
% % colormap(gcf,'hot')
% Surface3D(bound,density,Myfun);
% box on
% ----------------------------------------------------------------------

% NofSamples = 1e6;
%  [SBM.SamplePool,~] = RNgeneratorV4(ProSys,NofSamples);
%  hold on
% plot(SBM.SamplePool(:,1),SBM.SamplePool(:,2),'.','Color',[0.9,0.9,0.9])
   