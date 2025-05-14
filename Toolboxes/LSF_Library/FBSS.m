function y=FBSS(X,fun_par)
% Four branch series system
% input m*n ，m为数据点数；n为x的维数
    FBSS_type=fun_par.type;
    a=fun_par.a;
    b=fun_par.b;
    f0 = fun_par.f0;
    x1=X(:,1);x2=X(:,2);
    f1=a+0.1*(x1-x2).^2-(x1+x2)./sqrt(2);
    f2=a+0.1*(x1-x2).^2+(x1+x2)./sqrt(2);
    f3=(x1-x2)+b./sqrt(2);
    f4=(x2-x1)+b./sqrt(2);
    switch FBSS_type
        case 1   
            y=min([f1';f2';f3';f4']);
        case 2
             y=f1';
        case 3
             y=f2';
        case 4
             y=f3';
        case 5
             y=f4';
        case 6
            y1=min([f1';f2';f3';f4']);
            y2=max([f1';f2';f3';f4']);
            n=length(y1);
            for i=1:n
                if y1(i)>0
                    y(i)=y2(i);
                elseif y1(i)<=0
                    y(i)=y1(i);
                end
            end
    end
    y=y(:)-f0;
return