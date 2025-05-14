function  ALRM_PostProcessingV1(ALRMResult,PTO)
% PTO = [1     0     1     1];
% PostProcessing of Pf estiamation 
% PTO = [1,0,1,1,1]
%PTO1: 2D-plot of distribution of DoEs and LSF
%PTO2: 3D-plot of perfomance function (the ture and the predicted one)
%PTO3: CDF and CCDF
%PTO4: Plot the points (G_predict,G_ture)
%PTO5: CDF and CCDF in different iteration
%% Figure setting
set(0,'DefaultFigureColor',[1 1 1])
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultTextInterpreter','latex');
set(0,'defaulttextinterpreter','latex','DefaultAxesFontSize',15);

%% Input Parameter
SurrModelPar = ALRMResult.SurrModelPar;
ProSys = ALRMResult.ProSys;
SBM = ALRMResult.SBM;
DoE = SurrModelPar.DoE;
% text_ALRM = ALRMResult.text_ALRM;
ALSMTimeHis = ALRMResult.ALSMTimeHis;
ALSMPar = ALRMResult.ALSMPar;

%%  Figrue Par
% figure bound = [x1min x2min,x1max x2max]
vm = 5;
if size(ProSys.sigmaX,1)>1
    PDC.bound = [ProSys.muX-vm*diag(ProSys.sigmaX)',ProSys.muX+vm*diag(ProSys.sigmaX)'];
else
    PDC.bound = [ProSys.muX-vm*ProSys.sigmaX,ProSys.muX+vm*ProSys.sigmaX];
end
if ProSys.Distri == 5
    PDC.bound = [ProSys.a,ProSys.b];
end
PDC.gap = 100;
% Picture Drawing Control
PDC.PTO = PTO;

%% plot the contour plane only when dimension=2
P_Nov = 1;
% close all
if ProSys.Ndim == 2 
    xx = linspace(PDC.bound(1),PDC.bound(3),PDC.gap);
    yy = linspace(PDC.bound(2),PDC.bound(4),PDC.gap);
    [X1,X2] = meshgrid(xx,yy);
    X = [reshape(X1,PDC.gap*PDC.gap,1),reshape(X2,PDC.gap*PDC.gap,1)];   
    % Predict
    [YX_predict,YX_Gmse] = SurrModelPar.ModelPredictor(SurrModelPar,X);  
    YX_predict = reshape(YX_predict,PDC.gap,PDC.gap);
    YX_Gmse = reshape(YX_Gmse,PDC.gap,PDC.gap);
    % Exact
    if ALSMTimeHis.false_ture_statistics == 1
        YX_true = ProSys.fname(X,ProSys.fun_par);
        YX_true = reshape(YX_true,PDC.gap,PDC.gap);
    else
        YX_true  = YX_predict;
    end
end  

%% -----------------------------------------------------------------------------------------------
% ============================== Picture #1 ===============================
% --------------------------------------------------------------------------------------------------------------
if ProSys.Ndim==2&&PDC.PTO(1)==1  
        fg(P_Nov)=figure(P_Nov);P_Nov = P_Nov+1;
        set(gcf,'OuterPosition',[50 50 500 500]); 
       
        % 1. Plot the distribution of Candidate Pool
        hold on
        pl{1} = plot(SBM.SamplePool(:,1),SBM.SamplePool(:,2),'o','MarkerSize',1);
        pl{1}.Color = [0.83 0.82 0.78];
        xlabel('$x_1$')
        ylabel('$x_2$')
               
        % 2. Plot the Contour Plan
        % Exact
        pl{2} = contour(X1, X2, YX_true,[0,0],'color','g','linewidth',2,'LineStyle','-');
        % Predicted
        pl{3} = contour(X1, X2, YX_predict,[0,0],'color','r','linewidth',2,'LineStyle','--');
        
        % 3. Plot the DoE
        % Initial DoE
        hold on
        pl{4} = plot(DoE.X(1:ALSMPar.IniDoE.N0,1),DoE.X(1:ALSMPar.IniDoE.N0,2),'o','color',[178,34,34]/255);
        pl{4}.MarkerSize = 5; 
        pl{4}.MarkerEdgeColor = [178,34,34]/255; pl{4}.MarkerFaceColor = [178,34,34]/255; 
         % Added best DoE
        hold on
        pl{5} = plot(DoE.X(ALSMPar.IniDoE.N0+1:end,1),DoE.X(ALSMPar.IniDoE.N0+1:end,2),'b+');
        pl{5}.MarkerSize = 5;
      
         
        % legend       
        temp_str{1} = 'Contour';
        temp_str{2} = 'Sampling Points';
        temp_str{3} = 'Initial DoE';
        temp_str{4} = 'Added best DoE';    
        hl = legend(temp_str,'Location','best');
        hl.Box = 'off';

        % grid
        set(gca,'Box','on')
        grid off
end

%% -----------------------------------------------------------------------------------------------
% ============================ Picture #2 =================================
% --------------------------------------------------------------------------------------------------------------    
if ProSys.Ndim==2&&PDC.PTO(2)==1  
        fg(P_Nov)=figure(P_Nov);P_Nov = P_Nov+1;
        subplot(1,2,1)
        % 3D plot of exact model
         surfc(X1,X2,YX_true);
          polarmap([cmap('red',1,0,0);flipud(cmap('red',20,30,5))])
         polarmap(300,0.9)
         shading interp
         alpha 0.85
        title('Original model')
        xlabel('$x_1$');
        ylabel('$x_2$');
        zlabel('$\mathcal{M}(x_1,x_2)$');
        temprange = get(gca,'Zlim');
        box on
        set(gcf,'Units','centimeters','Position',[3,3,20,8]); 
        subplot(1,2,2)
        % 3D plot of metamode
        
        YX_predict(find(YX_predict>temprange(2)))=nan;
         YX_predict(find(YX_predict<=temprange(1)))=nan;
        surfc(X1,X2,YX_predict);
         polarmap([cmap('red',1,0,0);flipud(cmap('red',20,30,5))])
         polarmap(300,0.9)
         shading interp
         alpha 0.85
         hold on;
         plot3(DoE.X(:,1),DoE.X(:,2),DoE.Y, 'o', 'MarkerEdgeColor',[178,34,34]/255, 'MarkerFaceColor',[178,34,34]/255, 'MarkerSize',3);
         hold on
         ax = gca;
         temp = ax.ZLim(1);
         plot3(DoE.X(:,1),DoE.X(:,2),DoE.Y.*0+temp, 'o', 'MarkerEdgeColor',[178,34,34]/255, 'MarkerFaceColor',[178,34,34]/255, 'MarkerSize',3);
         box on
         xlabel('$x_1$');
         ylabel('$x_2$');
         zlabel('$\hat \mathcal{M}(x_1,x_2)$');
         title('Metamodel')      
          
        % 3D plot of variance surface
        figure(P_Nov);P_Nov = P_Nov+1;
        surfc(X1,X2,sqrt(YX_Gmse));
        polarmap([cmap('red',1,0,0);flipud(cmap('red',20,30,5))])
        polarmap(300,0.9)
        shading interp
        alpha 0.85
        title('Variance')
end
%% -----------------------------------------------------------------------------------------------
% ============================ Picture #3 =================================
% ------------------------------------------------------------------------------------------------------------
if PDC.PTO(3)==1   
    fg(P_Nov)=figure(P_Nov);P_Nov = P_Nov+1;
    set(gcf,'Units','centimeters','Position',[0,0,15,16]);
    subplot(2,1,1)
    nn = numel(ALSMTimeHis.Pf_predict);
    plot(1:nn,ALSMTimeHis.Pf_predict,'--ob');
    grid on
    xlabel('Iter')
    ylabel('$\hat{P}_f$')
    title('Metamodel')
    
    subplot(2,1,2)
    nn = numel(ALSMTimeHis.Pf_ture);
    plot(1:nn,ALSMTimeHis.Pf_ture,'ob--')
    grid on
    xlabel('Iter')
    ylabel('${P}_f$')
    title('Crude Numerical Method')
end
%% -----------------------------------------------------------------------------------------------
% ============================ Picture #4 =================================
% -------------------------------------------------------------------------------------------------------------- 
%% Plot the points (G_predict,G_ture)
if PDC.PTO(4)==1&&ALSMTimeHis.false_ture_statistics == 1
    fg(P_Nov)=figure(P_Nov);P_Nov = P_Nov+1;
    % Preparation of data
    if ALSMTimeHis.false_ture_statistics == 0
        tic
         Y_true = ProSys.fname(SBM.SamplePool,ProSys.fun_par);
        load data_Ytrue
        toc
    else
        Y_true = SBM.G_ture;
    end    
    Y_predict = SBM.G_predict;
    
    % Plot the pairs
    plot(Y_predict,Y_true,'.','MarkerEdgeColor',[0,113,188]/255, 'MarkerFaceColor',[0,113,188]/255)
    alpha(0.1)
    y_max = max(Y_true);
    y_min = min(Y_true);
    hold on
    plot([y_min,y_max],[y_min,y_max],'-','Color',[178,34,34]/255,'LineWidth',2)
    hold off
    axis([y_min,y_max,y_min,y_max])
    grid on  
    xlabel('$\hat{\mathcal{M}}(\mathbf{x})$','FontSize',15);
    ylabel('$\mathcal{M}(\mathbf{x})$','FontSize',15);
end
return
%% 
% % save picture
% figpath = 'C:\software\matlab2016\bin\dsmj\¿É¿¿¶È\ALRM\Disscusion_Correlation\AL';
% figpath = cd;
% TestExample = 'eg4_100intervals'; 
% figname1 = [TestExample,'-','SampleDistri and Contour'];
% figname2 = [TestExample,'-','Original Model and metamodel'];
% figname3 = [TestExample,'-','CDF and CCDF'];
% figname4 = [TestExample,'-','Gpredicct and Gture'];
% figname5 = [TestExample,'-','iteration of CDFandCCDF'];
% 
% saveas(fg(1),[figpath,'\',figname1,'.png']); 
% saveas(fg(2),[figpath,'\',figname2,'.png']);
% saveas(fg(3),[figpath,'\',figname3,'.png']);saveas(fg(3),[figpath,'\',figname3,'.fig']);
% saveas(fg(4),[figpath,'\',figname4,'.png']);saveas(fg(4),[figpath,'\',figname4,'.fig']);
% saveas(fg(5),[figpath,'\',figname5,'.png']);

%%
[XL,yl,XS,YS,beta,PCTVAR,MSE,stats] = plsregress(SamplePool,Y,10);