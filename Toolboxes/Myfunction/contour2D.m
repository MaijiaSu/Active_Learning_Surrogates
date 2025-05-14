%     [YX_predict,YX_Gmse] = SurrModelPar.ModelPredictor(SurrModelPar,X,);  
function [X1, X2, YX] = contour2D(bound,gap,MyFun)
    set(0,'DefaultFigureColor',[1 1 1])
    set(groot, 'defaultAxesTickLabelInterpreter','latex');
    set(groot, 'defaultLegendInterpreter','latex');
    set(groot, 'defaultTextInterpreter','latex');
    set(0,'defaulttextinterpreter','latex','DefaultAxesFontSize',15)
% bound = [ProSys.muX-vm*diag(ProSys.sigmaX)',ProSys.muX+vm*diag(ProSys.sigmaX)'];
% gap = 100;
xx = linspace(bound(1),bound(3),gap);
yy = linspace(bound(2),bound(4),gap);
[X1,X2] = meshgrid(xx,yy);
X = [reshape(X1,gap*gap,1),reshape(X2,gap*gap,1)];
YX = MyFun(X);
YX = reshape(YX,gap,gap);
% figure(1)
contour(X1, X2, YX,20)
polarmap([cmap('red',1,0,0);flipud(cmap('red',20,30,5))])
polarmap(300,0.9)
alpha(1)

end
