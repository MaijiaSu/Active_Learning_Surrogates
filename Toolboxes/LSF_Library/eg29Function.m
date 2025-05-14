function y=eg29Function(u,fun_par)
% dynamic response of a non-linear oscillator
% input m*n ，m为数据点数；n为x的维数
% fun_par =[] 无需求 
% output y为列向量
      
M1 = u(:,1);
M2 = u(:,2);
M3 = u(:,3);
S = u(:,4);

% -------------------------------------------------------------------------
g1 = 2*M1+2*M3-4.5*S;
g2 = 2*M1+M2+M3-4.5*S;
g3 = M1+M2+2*M3-4.5*S;
g4 = M1+2*M2+M3-4.5*S;
y = min([g1,g2,g3,g4]');
y = y(:);
% -------------------------------------------------------------------------

return