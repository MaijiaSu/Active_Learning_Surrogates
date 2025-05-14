function [G,dy] = SFstru_sinInput(CandidatePool,fun_par)

% X = normrnd(0,1,1,4);
fun_par = [];
Ns = size(CandidatePool,1);

for jj = 1:Ns
    %% uncertainty variables
    X = CandidatePool(jj,:);
    
    %% Structure Parameter
    % parameter of Shear-frame structure
    k = [3.0,2.8,1.5]*1e8; % N/m
    m = [1,1,1]*1e6;    % kg
    kexi = 0.05;        % Rayleigh damping
    DoF = 3;
    % Bouc-wen model
    alpha1 = 0.1; alpha2 = 0.1;alpha3 = 0.1;
    uy1 = 0.04; uy2 = 0.04; uy3 = 0.04;
    n1 = 5; n2 = 5; n3 = 5;
    A1 = 1; A2 = 1; A3 = 1;
    yita1 = 1/(2*uy1^n1); yita2 = 1/(2*uy2^n2); yita3 = 1/(2*uy3^n3);
    gamma1 = yita1; gamma2 = yita2; gamma3 = yita3;
    
    %% matrix of K M C
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
    wi = wn(1,1); wj = wn(2,2);
    ab = 2*wi*wj/(wj^2-wi^2)*[wj , -wi; -1/wj, 1/wi] * [kexi;kexi];
    C = ab(1)*M+ab(2)*K;
    %% STATE-SPACE MATRICES
    kalpha1 = k.*[alpha1,alpha2,alpha3];       % kalha2 = alpha*k
    kalpha2 = k.*[1-alpha1,1-alpha2,1-alpha3]; % kalha2 = (1-alpha)*k
    Kalpha1=zeros(DoF); Kalpha2=zeros(DoF);
    for i=1:DoF-1
        Kalpha1(i,i)=kalpha1(i)+kalpha1(i+1);
        Kalpha1(i,i+1)=-kalpha1(i+1);
        Kalpha1(i+1,i)=-kalpha1(i+1);
        
        Kalpha2(i,i)=kalpha2(i)+kalpha2(i+1);
        Kalpha2(i,i+1)=-kalpha2(i+1);
        Kalpha2(i+1,i)=-kalpha2(i+1);
    end
    Kalpha1(DoF,DoF) = kalpha1(DoF);
    Kalpha2(DoF,DoF) = kalpha2(DoF);
    
    Abw = [zeros(DoF) ,   eye(DoF);
        -inv(M)*Kalpha1, -inv(M)*C];
    Bbw = [zeros(DoF,2*DoF);
        inv(M),-inv(M)*Kalpha2];
    % output of relative displacment and velocity
    temp_matrix = [1 0 0;
        -1 1 0;
        0 -1 1];
    Cbw = [temp_matrix,zeros(DoF);
        zeros(DoF),temp_matrix];
    Dbw = zeros(2*DoF,2*DoF);
    %% Input of sysetem
    % Input 1
    % Tsim = 20;
    % dt = 0.002;
    % Fs = 1/dt;
    % t = 0:dt:Tsim;
    % S0=0.05; %medium excitation
    % g = 9.81; % Gravity [m/sec2]
    % NoisePower = 2*pi*S0/g^2;
    % rand_seed=ceil(100000*(0.1+0.9*rand(1)));
    
    % Input 2
    Tsim = 10;
    dt = 0.01;
    t = 0:dt:Tsim;
    sinWave =1/6*(X(1)*sin(2*pi*t)+X(2)*sin(4*pi*t)+X(3)*sin(8*pi*t)+X(4)*sin(16*pi*t));
    f = m'*sinWave;
    sinInput.time=t;
    sinInput.signals.values=f';
    sinInput.signals.dimensions=DoF;
    %% Model construction and solver
    options = simset('SrcWorkspace','current');
    t = sim('ShearFrameStru_3F_sinInput',[],options);
    
    %% results
    d1 = y.data(:,1);
    d2 = y.data(:,2);
    d3 = y.data(:,3);
    dy = [d1';d2';d3'];
    G(jj,1) = max(max(abs(dy)));
    
    if mod(jj,100)==0
        disp(['Running:SFstru_sinInput--',num2str(jj),'/',num2str(Ns)])
    end
end



% z = [z1.data,z2.data,z3.data]';
% % Calculate story forces
% R = Kalpha1*dy+Kalpha2*z;
% % Calculate restoring force due to the stiffness of the ith story
% F = diag(kalpha1)*dy+diag(kalpha2)*z;
end
%%
% figure
% for ii = 1:3
%     subplot(3,1,ii)
%     plot(dy(i,:),F(i,:))
%     hold on
% end
% figure
% for ii = 1:3
% %     subplot(3,1,ii)
%     plot(t,dy(ii,:))
%     hold on
% end
% legend('1F','2F','3F')


