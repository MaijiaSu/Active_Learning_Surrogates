function y = JumpDiscotinuityNdim2n(X,fun_par)

x1 = X(:,1);
x2 = X(:,2);

 y = -(1.8*x1-0.1*x2.^3-5);
 
%  y = -(1*x1.^2-2*x2.^2-9);
 
 index = find(x1<0);
 
 y(index) = -(1.2*x1(index).^2-2*x2(index).^2-12);
 
end