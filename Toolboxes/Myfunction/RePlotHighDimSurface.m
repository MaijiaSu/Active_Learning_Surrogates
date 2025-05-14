function RePlotHighDimSurface(data,bound)

Ndim = size(bound,1);
FigID = [];
FigNov = [];
for ii = 1:Ndim
    for jj = 1:ii-1
        FigID = [FigID;ii,jj];
        FigNov = [FigNov;(ii-1)*Ndim+jj];
    end
end
NofFig = numel(FigNov)+Ndim;
%%
set(0,'DefaultFigureColor',[1 1 1])
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
set(0,'defaulttextinterpreter','latex','DefaultAxesFontSize',6)
fig_width = 300*Ndim; fig_width = min(fig_width,1200);
fig_height = 300*Ndim; fig_height = min(fig_height,1200);
f1 = figure('OuterPosition',[-fig_width/2 -fig_height/2 fig_width fig_height]);
for ii = 1:NofFig-Ndim
    subplot(Ndim,Ndim,FigNov(ii))
    fig_id = FigID(ii,1:2);
    
    surfc(data(fig_id(1),fig_id(2)).X1,data(fig_id(1),fig_id(2)).X2,data(fig_id(1),fig_id(2)).YY);
    polarmap([cmap('red',1,0,0);flipud(cmap('red',20,30,5))])
    polarmap(300,0.9)
    shading interp
    alpha 0.7
%     view(0,90)
    
    xlabel(['$x_',num2str(fig_id(1)),'$']);
    ylabel(['$x_',num2str(fig_id(2)),'$']);
end
tightfig;
FigPos = get(f1,'position');
set (f1,'position',[FigPos(3:4)/2,FigPos(3:4)])