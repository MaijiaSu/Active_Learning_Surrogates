function G = StructModel_BLWN(X,fun_par)

Ns = size(X,1);
G = zeros(Ns,1);

% system parmeter
DoF=3;                             % degree of freedom
k=[200*1e3,200*1e3,150*1e3];       % story stiffness,unit: N/m
m=[1200,1200,1200];                % story mass,unit: kg
model_kexi=[0.02 0.02 0.04];       % moadl dapming

%%
% matrix of K M C
    K=zeros(DoF);
    for i=1:DoF-1
        K(i,i)=k(i)+k(i+1);
        K(i,i+1)=-k(i+1);
        K(i+1,i)=-k(i+1);
    end
    K(DoF,DoF) = k(DoF);
    
    M = diag(m);

% modal analyze
    [phi,lamda] = eig(K,M);
% rearrange the phi
    [lamda,lamda_index] = sort(diag(lamda),'ascend');
    phi= phi(:,lamda_index);
% normalization
    temp = sqrt(sum(phi.^2))'*ones(1,DoF);
    phi = phi./temp;  
    wn = sqrt(diag(lamda));
    Me = phi'*M*phi;
    Ce = 2*diag(model_kexi).*wn.*Me;
    C = phi*Ce*phi'; 

% sapce state
    As = [zeros(3,3),eye(3,3);
          -M\K,     -M\C  ];
    Bs = [zeros(3,1);
        -ones(3,1)];
    Cs = [eye(DoF)+diag(-ones(DoF-1,1),-1),zeros(DoF)];
    % namely
    % Cs = [1 0 0 0 0 0;
    %       -1 1 0 0 0 0;
    %       0 -1 1 0 0 0];
    Ds = zeros(DoF,1);

% excitation
    Tsim = 10;
    dtint = 2*pi/wn(DoF,DoF)/40;
    BLWNpower = 1;
    dtBLWN = floor(0.1/dtint)*dtint;
    BLWN = [];
    Nstep = ceil(Tsim/dtBLWN);
    t_BLWN = 0:dtBLWN:dtBLWN*Nstep;
    t = 0:dtint:Tsim;   
    
for ii = 1:Ns
    
    BLWN = X(ii,:);               % random vector
    
%% slover
% simulink
%     sim('Sim_excitedBLWN.slx')
    
% Newmark-beta

    Gace = interp1(t_BLWN,[0,BLWN],t);  
    for i = 1:DoF
        F(i,:) = -m(i)*Gace;
    end
    beta=1/6;gama=1/2;
    Y0=zeros(DoF,2);
    [x,v,a]=Newmark_dsmj(K,M,C,F,dtint,Y0,beta,gama);
    rel_D(1,:) = x(1,:);
    rel_D(2:3,:) = x(2:3,:)-x(1:2,:);
%% Perform function
    G(ii) = max(max(abs(rel_D)));
    if mod(ii,100) == 0
        disp([num2str(ii),'/',num2str(Ns)]);
    end  
end

end

%% Figrue
% subplot(311)
% plot(simout.time,simout.data(:,1))
% ylabel('r_d1/(m)')
% subplot(312)
% plot(simout.time,simout.data(:,2))
% ylabel('r_d2/(m)')
% subplot(313)
% plot(simout.time,simout.data(:,1))
% ylabel('r_d3/(m)')
% xlabel('time/s')
% %  BLWN
% figure
% plot(BLWN)
% ylabel('ground acceleration/(m/s^2)')
% xlabel('time/s')

%% Monte Carlo
% Tsim = 10
% Ncall = 1000
% for ii=1:Ncall
%     tic
%     disp([num2str(ii) '/' num2str(Ncall)])
%     rand_seed=ceil(100000*(0.1+0.9*rand(1)));
%     % 
%     sim('Sim_excitedBLWN.slx')
%     RU(:,:,ii) = simout.data;
%     BLWNhis(:,ii) = BLWN.data;
%     toc
% end
% 
% for ii = 1:DoF
%     data = squeeze(RU(:,ii,:));
%     Max_RD = max(abs(data));
%     figure
%     hist(Max_RD,30)
% end   