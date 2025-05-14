function data = PlotHighDimSurface(Myfun,bound,density)
% This function is to illustrate the response surface in multiple-dimesion
% problem 
% input Variables
% Myfun: Function handle, !!! Myfun should support vectorization input
% bound: bound of drawing range of each dimension; 
% bound = [x1_min,x1_max;
%          x2_min,x2_max;
%          ...];
% density: The density of discretization, i.e Grid density: 

%% Evaluation Points
Ndim = size(bound,1);
CenterPints = mean(bound,2);
for ii = 1:Ndim
    xx(ii,:) = linspace(bound(ii,1),bound(ii,2),density);
end  

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
set(0,'defaulttextinterpreter','latex','DefaultAxesFontSize',6)
fig_width = 300*Ndim; fig_width = min(fig_width,1200);
fig_height = 300*Ndim; fig_height = min(fig_height,1200);
% f1 = figure('OuterPosition',[-fig_width/2 -fig_height/2 fig_width fig_height]);

%% drawing graph at the lower triangle
for ii = 1:NofFig-Ndim
    subplot(Ndim,Ndim,FigNov(ii))
    fig_id = FigID(ii,1:2);
    disp(['fig_id:',num2str(fig_id(1)),'-',num2str(fig_id(2))]);
    [X1,X2] = meshgrid(xx(fig_id(1),:),xx(fig_id(2),:));
    X = [reshape(X1,density*density,1),reshape(X2,density*density,1)];
    X_input = ones(density*density,1)*CenterPints';
    X_input(:,fig_id(1)) =  X(:,1);
    X_input(:,fig_id(2)) =  X(:,2);  
    try
        [Y,Add_data] = Myfun(X_input);
    catch
        Y = Myfun(X_input);
        Add_data = [];
    end
    YY = reshape(Y,density,density);
    % plot the contour plan
    surfc(X1,X2,YY);
%     polarmap([cmap('red',1,0,0);flipud(cmap('red',20,30,5))])
%     polarmap(300,0.9)
    shading interp
    alpha 0.7
    view(0,90)
    xlabel(['$x_',num2str(fig_id(1)),'$']);
    ylabel(['$x_',num2str(fig_id(2)),'$']);
    box on
% Save data    
    data(fig_id(1),fig_id(2)).X1 = X1;
    data(fig_id(1),fig_id(2)).X2 = X2;
    data(fig_id(1),fig_id(2)).X_input = X_input;
    data(fig_id(1),fig_id(2)).YY = YY;
    data(fig_id(1),fig_id(2)).Add_data = Add_data;
end
tightfig;
% FigPos = get(f1,'position');
% set(f1,'position',[-FigPos(3:4)/2,FigPos(3:4)]) 