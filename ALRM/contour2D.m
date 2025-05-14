function [X1, X2, YX] = contour2D(bound,gap,MyFun)
    xx = linspace(bound(1),bound(3),gap);
    yy = linspace(bound(2),bound(4),gap);
    [X1,X2] = meshgrid(xx,yy);
    X = [reshape(X1,gap*gap,1),reshape(X2,gap*gap,1)];
    YX = MyFun(X);
    YX = reshape(YX,gap,gap);
    contour(X1, X2, YX,50)
    % mix-max
%     caxis([(min(YX(:))),max(YX(:))])
% max-min
%     caxis([-1,1]*max(abs(caxis)))
    polarmap([cmap('red',1,0,0);flipud(cmap('red',20,30,5))])
    polarmap(300,0.9)
    alpha(0.6)
end
