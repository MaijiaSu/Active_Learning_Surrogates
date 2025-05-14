%     [YX_predict,YX_Gmse] = SurrModelPar.ModelPredictor(SurrModelPar,X,);  
function [X1, X2, YX] = Surface3D(bound,gap,MyFun)
%     set(0,'DefaultFigureColor',[1 1 1])
%     set(groot, 'defaultAxesTickLabelInterpreter','latex');
%     set(groot, 'defaultLegendInterpreter','latex');
%     set(groot, 'defaultTextInterpreter','latex');
%     set(0,'defaulttextinterpreter','latex','DefaultAxesFontSize',15)
% bound = [ProSys.muX-vm*diag(ProSys.sigmaX)',ProSys.muX+vm*diag(ProSys.sigmaX)'];
% gap = 100;



 xx = linspace(bound(1,1),bound(1,2),gap);
    yy = linspace(bound(2,1),bound(2,2),gap);
[X1,X2] = meshgrid(xx,yy);
X = [reshape(X1,gap*gap,1),reshape(X2,gap*gap,1)];
YX = MyFun(X);
YX = reshape(YX,gap,gap);
% figure(1)
surfc(X1,X2,YX);
polarmap([cmap('red',1,0,0);flipud(cmap('red',20,30,5))])
polarmap(300,0.9)
shading interp
alpha 0.5

% mix-max
    caxis([(min(YX(:))),max(YX(:))])
% max-min
%     caxis([-1,1]*max(abs(caxis)))


% polarmap([cmap('red',1,0,0);flipud(cmap('red',20,30,5))])
% polarmap(300,0.9)
% alpha(1)
end
%     saveas(gcf,[figpath,'\',figname,'.png'])
%     saveas(gcf,[figpath,'\',figname,'.fig'])
%     saveas(gcf,[figpath,'\',figname,'.eps'])
