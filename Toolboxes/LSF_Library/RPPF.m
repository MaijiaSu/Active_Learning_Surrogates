function y=RPPF(u,fun_par)
% dynamic response of a non-linear oscillator
% input m*n ，m为数据点数；n为x的维数
% fun_par =[] 无需求 
% output y为列向量


    [Ns,Ndim] = size(u);
    
%     kk = 1/1.5;     % kkk = 1/1.5, pf = 0.0033;
    kk = 1;
    muX=[1 1 1 1 0.7/kk 1.0/kk];
    sigmaX=[0.15 0.15 0.15 0.15 0.119/kk 0.5/kk];
    
    X = zeros(Ns,Ndim);
    for n = 1:Ndim
        X(:,n) = u(:,n)*sigmaX(n) + muX(n);
    end

  M1=X(:,1);M2=X(:,2);M3=X(:,3);
  M4=X(:,4);H=X(:,5);V=X(:,6);
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   G1 = M1+2*M3+2*M4-H-V;
   G2 = M2+2*M3+M4-V;
   G3 = M1+M2+M4-H;
   G4 = M1+2*M2+2*M3-H+V;
   y=min([G1,G2,G3,G4]')';
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return