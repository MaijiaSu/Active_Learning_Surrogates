function G = IsoSysOf4F(u,fun_par)
%% 转化数据
% 将标准正态变量u映射回原始空间
    [Ns,Ndim] = size(u);
if strcmp(fun_par.StadNrom,'Yes')     
    % 对数正态分布参数
     muY = fun_par.muX;
     sigmaY = fun_par.sigmaX; 
    % 相应正态分布参数
    deta = sigmaY./muY;
    sigmaX = sqrt(log(1+deta.^2));
    muX = log(muY)-sigmaX.^2/2;    
    % 标转正态分布转换为变量X
    X = zeros(Ns,Ndim);
    for n = 1:Ndim
        X(:,n) = u(:,n)*sigmaX(n) + muX(n);
    end
    % 正态变量X转化为变量Y
    u = exp(X);
    
else  % 不需转化
    
    u = u;
    
end

u = [u(:,1:5),zeros(size(u,1),1),u(:,6:7)];

for iii = 1:size(u,1)
%% System Par
% random variable, unit: m, N, s
    X = u(iii,:);
% Floor Parameter
    mf = X(:,1);
    kf = X(:,2);
    cf = X(:,3);
% Isolation layer parmeter
    dy = X(:,4);     % yield disaplacement
    fy = X(:,5);     % yield force
% Equipment Floor parmeter
    m1 = 500;
    m2 = 100;
    k1 = 2500;
    k2 = 105;
    c1 = 350;
    c2 = 200;
% Excitation Parmeter
    T = X(:,7);  % Pulse excitation period 
    A = X(:,8);  % Pulse ampltitude
%% Structure system with isolation layer
% Structure system 
    NofFloor = 4;
    m = ones(1,NofFloor)*mf;
    k = [fy/dy,ones(1,NofFloor-1)*kf];
    c = [cf,ones(1,NofFloor)*cf];
    M_stru = diag(m);
    for i=1:NofFloor-1
        K_stru(i,i)=k(i)+k(i+1);
        K_stru(i,i+1)=-k(i+1);
        K_stru(i+1,i)=-k(i+1);
    end
    K_stru(NofFloor,NofFloor) = k(NofFloor);
    for i=1:NofFloor-1
        C_stru(i,i)=c(i)+c(i+1);
        C_stru(i,i+1)=-c(i+1);
        C_stru(i+1,i)=-c(i+1);
    end
    C_stru(NofFloor,NofFloor) = c(NofFloor);
% Equipment isolation system
    M_iso = [m1 , 0
             0 , m2];
    K_iso = [k1+k2 , 0;
               0   , k2];
    C_iso = [c1+c2 , 0;
               0   , c2];       
% Assemble the two system
    NofIsoFloor = 2; % number of equipment isolation Floor 
% Martix M
    M = [M_stru, zeros(NofFloor,2);
         zeros(2,NofFloor), M_iso];
% Martix K 
    interfaceK = zeros(NofFloor,2);
    interfaceK(NofIsoFloor,1) = -k1;
    K = [K_stru, interfaceK;
         interfaceK', K_iso];
% Martix C
    interfaceC = zeros(NofFloor,2);
    interfaceC(NofIsoFloor,1) = -c1;
    C = [C_stru, interfaceC;
         interfaceC', C_iso];
% Time history of exciation
%     dt_inter = 0.01;
%     t = 0:dt_inter:10;
%     Gacce = zeros(size(t));
%     tempT = 0:dt_inter:T;
%     Gacce(1:size(tempT,2)) = A*sin(2*pi/T*tempT);
    
    % EIcenrto earthquake wave input
    Gacce = load('elcentro.txt');
    Gacce = Gacce'/max(abs(Gacce))*A*2;
    dt_inter = 0.02;
    t = 0:dt_inter:(length(Gacce)-1)*dt_inter;
    
%% Slover
    % Newmark-beta 
    mm = [ones(1,NofFloor)*mf,m1,m2];
    F=zeros(size(M,1),length(Gacce)); %单位N
    for i = 1:6
        F(i,:) = -mm(i)*Gacce;
    end    
    beta=1/6;gama=1/2;
    Y0=zeros(size(M,1),2);
    [x,v,a]=Newmark_dsmj(K,M,C,F,dt_inter,Y0,beta,gama);

%% Perform function
    index1 = 0.04 - max(max(abs(x(2:4,:)-x(1:3,:))'));
    index2 = 0.5 - max(abs(Gacce+a(6,:)));
    index3 = 0.25 - max(abs(x(5,:)-x(2,:)));
    G(iii,1) = 12.5*index1 + index2 + 2*index3;
    
    if mod(iii,100)==0
        disp(['Running:IsoSysOf4F--',num2str(iii),'/',num2str(Ns)])
    end
end

% disaplacement
% figure
% subplot(3,1,1)
% plot(t,Gacce)
% subplot(3,1,2)
% plot(t,x)
% hl = legend('isoF','1F','2F','3F','Eiso','Equipment');
% hl.Box = 'off';
% subplot(3,1,3)
% plot(t,a(1:6,:)+Gacce)
% hl = legend('isoF','1F','2F','3F','Eiso','Equipment')
% hl.Box = 'off';

% 
% figure
% subplot(3,1,1)
% plot(t,[x(1,:);x(2:4,:)-x(1:3,:)])
% hl = legend('Drift isoF','Drift 1F','Drift 2F','Drift 3F');
% hl.Box = 'off';
% subplot(3,1,2)
% plot(t,a(1:6,:)+Gacce)
% hl = legend('abs acce of isoF','1F','2F','3F','Eiso','Equipment')
% hl.Box = 'off';

% % modal analyze
% [phi,lamda] = eig(K,M);
% % normalization
% temp = sqrt(sum(phi.^2))'*ones(1,size(M,1));
% phi = phi./temp';  
% wn = sqrt(diag(lamda));

% 计算加速度放大系数
% w = 2*pi/T;
% r = (w/wn(1))^2;
% kexi = cf/ (2*sum(diag(M))*wn(1));
% Ra = sqrt( (1+4*kexi.^2*r)/((1-r)^2+4*kexi.^2*r) );

%% sapce state
% DOF = size(M,1);
% As = [zeros(DOF,DOF),eye(DOF,DOF);
%       -M\K,     -M\C  ];
% Bs = [zeros(DOF,1);
%     -ones(DOF,1)];
% Cs = [ -M\K, -M\C]; 
% Ds = zeros(DOF ,1);
% 
% %
% Tsim = 10;
% dtint = 2*pi/wn(DOF)/40;
% 
% %
% sim('Sim_IsoSysOf4F')
% 
% subplot(2,1,1)
% plot(u_simulink.time,u_simulink.data)
% subplot(2,1,2)
% plot(t,a(1:6,:)+Gacce)
% hl = legend('abs acce of isoF','1F','2F','3F','Eiso','Equipment');
% hl.Box = 'off';
