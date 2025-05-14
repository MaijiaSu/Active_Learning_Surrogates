function G = BoreholeFun(X,P)
% borehole-function
rw = X(:,1);
r = X(:,2);
Tu = X(:,3);
Hu = X(:,4);
Tl = X(:,5);
Hl = X(:,6);
L = X(:,7);
kw = X(:,8);

G = 2*pi*Tu.*(Hu-Hl)./log(r./rw)./...
            (1+2*L.*Tu./log(r./rw)./rw.^2./kw+Tu./Tl);   

end

