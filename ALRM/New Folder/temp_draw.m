clc,clear,close all
TestExample = 'eg1';                           % select the example from LSF_library, LimtStateFunction_select.m;
SBM.method = 'MCS_UQ';
LimtStateFunction_select
figpath = 'C:\software\matlab2016\bin\dsmj\¿É¿¿¶È\ALRM\Result_Projection of the function';
figname = TestExample;
save_figures = 1;
%%
Ndim = ProSys.Ndim;
density = 100;
VM = 5;
% bound = ones(3,1)*[-pi,pi];
muX = ProSys.muX';
sigmaX = diag(ProSys.sigmaX);

bound = [(muX-VM*sigmaX),(muX+VM*sigmaX)];

for ii = 1:Ndim
    xx(ii,:) = linspace(bound(ii,1),bound(ii,2),density);
end  
%     xx = linspace(PDC.bound(1),PDC.bound(3),PDC.gap);
%     yy = linspace(PDC.bound(2),PDC.bound(4),PDC.gap);
%     [X1,X2] = meshgrid(xx,yy);
%     X = [reshape(X1,PDC.gap*PDC.gap,1),reshape(X2,PDC.gap*PDC.gap,1)];

%% distribution of Fig 
FigID = [];
FigNov = [];
for ii = 1:Ndim
    for jj = 1:ii-1
        FigID = [FigID;ii,jj];
        FigNov = [FigNov;(ii-1)*Ndim+jj];
    end
end
NofFig = numel(FigNov)+Ndim;

%% Figure setting
set(0,'DefaultFigureColor',[1 1 1])
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
set(0,'defaulttextinterpreter','latex','DefaultAxesFontSize',12)
fig_width = 300*Ndim; fig_width = min(fig_width,1200);
fig_height = 300*Ndim; fig_height = min(fig_height,1200);
f1 = figure('OuterPosition',[-fig_width/2 -fig_height/2 fig_width fig_height]);

%% drawing graph at the lower triangle
for ii = 1:NofFig-Ndim
    subplot(Ndim,Ndim,FigNov(ii))
    fig_id = FigID(ii,1:2);
    [X1,X2] = meshgrid(xx(fig_id(1),:),xx(fig_id(2),:));
    X = [reshape(X1,density*density,1),reshape(X2,density*density,1)];
    X_input = ones(density*density,1)*ProSys.muX;
    X_input(:,fig_id(1)) =  X(:,1);
    X_input(:,fig_id(2)) =  X(:,2);  
    Y = ProSys.fname(X_input,ProSys.fun_par);
    YY = reshape(Y,density,density);
    % plot the contour plan
    surfc(X1,X2,YY);
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
set (f1,'position',[-FigPos(3:4)/2,FigPos(3:4)]) 

%%  drawing graph at the diagonal
%     ax = gca;
% outerpos = ax.OuterPosition;
% ti = ax.TightInset; 
% left = outerpos(1) + 1*ti(1);
% bottom = outerpos(2) + 1*ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];
if save_figures == 1
    saveas(f1,[figpath,'\',figname,'.png'])
end
