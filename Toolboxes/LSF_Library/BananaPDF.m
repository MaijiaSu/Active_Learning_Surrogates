function bPDF = BananaPDF(X,Par)
if size(X,2) == 1
    X = X';
end

x = X(:,1);
y = X(:,2);

% Normalization coefficient
Z = 2.7320;
% Par = [1.15,0.5,0.9];
a = Par(1);
b = Par(2);
po = Par(3);

temp1 = x.^2/a^2+a^2*(y-b*x.^2/a^2-b*a^2).^2-2*po*x/a.*(y-b*x.^2/a^2-b*a^2);
temp2 = -1/2/(1-po^2)*(temp1);
bPDF = 1/Z*exp(temp2);
end
