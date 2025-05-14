function bPDF = MistureBananaPDF(X,Par)
if size(X,2) == 1
    X = X';
end

x1 = X(:,1);
x2 = X(:,2);

% Normalization coefficient
Z = 1;
a = Par(1);
b1 = Par(2);
b2 = Par(3);
po = Par(4);

temp1 = x1.^2/a^2+a^2*(x2-b1*x1.^2/a^2-b1*a^2).^2-2*po*x1/a.*(x2-b1*x1.^2/a^2-b1*a^2);
temp2 = x1.^2/a^2+a^2*(x2-b2*x1.^2/a^2-b2*a^2).^2-2*po*x1/a.*(x2-b2*x1.^2/a^2-b2*a^2);

temp3 = -1/2/(1-po^2)*(temp1);
temp4 = -1/2/(1-po^2)*(temp2);


bPDF = 1/Z*(exp(temp3)+exp(temp4));


% % % % % x1 = X(:,1);
% % % % % x2 = X(:,2);
% % % % % 
% % % % % % Normalization coefficient
% % % % % Z = 2.7320;
% % % % % % Par = [1.15,0.5,0.9];
% % % % % a = Par(1);
% % % % % b = Par(2);
% % % % % po = Par(3);
% % % % % 
% % % % % temp1 = x1.^2/a^2+a^2*(x2-b*x1.^2/a^2-b*a^2).^2-2*po*x1/a.*(x2-b*x1.^2/a^2-b*a^2);
% % % % % temp2 = -1/2/(1-po^2)*(temp1);
% % % % % bPDF = 1/Z*exp(temp2);

end
