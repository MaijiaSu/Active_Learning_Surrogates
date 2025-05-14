function  ALRM_PostProcessingV2(ALRMResult,PDC)
% PostProcessing of CDF/CCDF estiamation 
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

if ~isfield(PDC,'bound')
    vm = 5;
    if size(ProSys.sigmaX,1)>1
        PDC.bound = [ProSys.muX-vm*diag(ProSys.sigmaX)',ProSys.muX+vm*diag(ProSys.sigmaX)'];
    else
        PDC.bound = [ProSys.muX-vm*ProSys.sigmaX,ProSys.muX+vm*ProSys.sigmaX];
    end
    if ProSys.Distri == 5
        PDC.bound = [ProSys.a,ProSys.b];
    end
end
if ~isfield(PDC,'gap')
    PDC.gap = 100;
end
% Picture Drawing Control
if ~isfield(PDC,'PTO')
    PDC.PTO = [1,1,1,1,1,1,1];
end 
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
        set(gcf,'OuterPosition',[50 50 1000 500]); 
        % 2. Plot the Contour Plane
        % The exact Contour Plan
        subplot(1,2,1)
        contour(X1, X2, YX_true,20)
        polarmap([cmap('red',1,0,0);flipud(cmap('red',20,30,5))])
        polarmap(300,0.9)
        alpha(0.7)
%       Plot the distribution of Candidate Pool
        hold on
        pl(1) = plot(SBM.SamplePool(:,1),SBM.SamplePool(:,2),'o','MarkerSize',1);
        pl(1).Color = [0.83 0.82 0.78];
        xlabel('$x_1$')
        ylabel('$x_2$')
        title('Original Model')
          
        subplot(1,2,2)
%         YX_predict(find(YX_predict>=20))=nan;
%         YX_predict(find(YX_predict<=-50))=-nan;
        contour(X1, X2, YX_predict,20)
        % Plot the distribution of Candidate Pool
        hold on
        pl(2) = plot(SBM.SamplePool(:,1),SBM.SamplePool(:,2),'o','MarkerSize',1);
        pl(2).Color = [0.83 0.82 0.78];
        hold on
        pl(3) = plot(DoE.X(1:ALSMPar.IniDoE.N0,1),DoE.X(1:ALSMPar.IniDoE.N0,2),'o','color',[178,34,34]/255);
        pl(3).MarkerSize = 3; 
        pl(3).MarkerEdgeColor = [178,34,34]/255; pl(3).MarkerFaceColor = [178,34,34]/255; 
        hold on
         
        if size(DoE.X,1)> ALSMPar.IniDoE.N0
            pl(4) = plot(DoE.X(ALSMPar.IniDoE.N0+1:end,1),DoE.X(ALSMPar.IniDoE.N0+1:end,2),'b+');
            pl(4).MarkerSize = 5;
        end

        xlabel('$x_1$')
        ylabel('$x_2$')
        title('Metamodel')
         
        % legend       
        temp_str{1} = 'Contour';
        temp_str{2} = 'Sampling Points';
        temp_str{3} = 'Initial DoE';
        temp_str{4} = 'Added best DoE';    
        hl = legend(temp_str,'Location','best');
        hl.Box = 'on';

        % grid
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
        
        % 3D plot of metamodel
        subplot(1,2,2)
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
%         figure(P_Nov);P_Nov = P_Nov+1;
%         surfc(X1,X2,sqrt(YX_Gmse));
% %         polarmap([cmap('red',1,0,0);flipud(cmap('red',20,30,5))])
% %         polarmap(300,0.9)
%         shading interp
%         alpha 0.85
%         title('Variance')
elseif (ProSys.Ndim>2&&PDC.PTO(2)==1)&&0
    % For the case the target function doesn't support Vectorized Operation
    % data = PlotHighDimSurface(@(X) FunTran_SupportVectorized(X,Myfun),bound,density)
    % For the case the target function does support Vectorized Operation
    % % Bound, density = 100;
    density = 30;
    VM = 3;
    muX = ProSys.muX';
    sigmaX = ProSys.sigmaX;
    bound = [(muX(:)-VM*sigmaX(:)),(muX(:)+VM*sigmaX(:))];
    
    Myfun_predict = @(X) SurrModelPar.ModelPredictor(SurrModelPar,X);
    Myfun_ture = @(X) ProSys.fname(X,ProSys.fun_par);
    figure
    data1 = PlotHighDimSurface(Myfun_predict,bound,density);
    figure
    data2 = PlotHighDimSurface(Myfun_ture,bound,density);
end
%% -----------------------------------------------------------------------------------------------
% ============================ Picture #3 =================================
% ------------------------------------------------------------------------------------------------------------
if PDC.PTO(3)==1   
    ymin = ProSys.ymin;
    ymax = ProSys.ymax;
    BW = (ymax-ymin)/SBM.NofInterval;
    yy = ymin:BW:ymax;
    mu_FY = ALSMTimeHis.CDF(end,:);
    % CDF
     fg(P_Nov)=figure(P_Nov);P_Nov = P_Nov+1;
     set(gcf,'Units','centimeters','Position',[3,3,22,10]); 
     subplot(1,2,1)
    %  3 sigma confidence interval 
    maxY = ALSMTimeHis.CDF3(end,:);
    minY = ALSMTimeHis.CDF1(end,:);  
    maxY(find(maxY<=0))=1/SBM.NofSamples;
    minY(find(minY<=0))=1/SBM.NofSamples;
    yFill = [maxY, fliplr(minY)];
    xFill = [yy, fliplr(yy)];  
    ll1 = fill(xFill, yFill, [0.83 0.82 0.78]); 
    ll1.LineStyle = 'None';
    % Exact (MCS) and predicted CDF
    hold on
    plot(yy,SBM.CDF,'-.','LineWidth',2)
    hold on
    plot(yy,SBM.CDF_ture,'-','LineWidth',2)
    ax1 = gca;
    ax1.YScale = 'log';
    ax1.XLim = [ymin,ymax];
    grid on
    box on
    hl1 = legend('$Metamodel^{+, -}$','$Metamodel$','$Exact(MCS)$');
    hl1.Box = 'off';
    hl1.Position = [0.1979 0.1557 0.2512 0.1937];
    set(ax1,'XColor','b','YColor','b','Box','on');
    set(gca,'fontsize',15,'fontname','Times');
    
    % title
    title('CDF')

    % CCDF
    subplot(1,2,2)
    %  3 sigma confidence interval 
    maxY = 1-ALSMTimeHis.CDF1(end,:);
    minY = 1-ALSMTimeHis.CDF3(end,:);
    maxY(find(maxY<=0))=1/SBM.NofSamples;
    minY(find(minY<=0))=1/SBM.NofSamples;
    yFill = [maxY, fliplr(minY)];
    xFill = [yy, fliplr(yy)];  % fliplr函数：左右翻转
    ll2 = fill(xFill, yFill, [0.83 0.82 0.78]); 
    ll2.LineStyle = 'None';
    % Exact (MCS) and predicted CDF
    hold on
    plot(yy,SBM.CCDF,'-.','LineWidth',2)
    hold on
    plot(yy,SBM.CCDF_ture,'-','LineWidth',2)
    ax2 = gca;
    ax2.YScale = 'log';
    ax2.XLim = [ymin,ymax];
    grid on
    box on
    hl2 = legend('$Metamodel^{+, -}$','$Metamodel$','$Exact(MCS)$');
    hl2.Box = 'off';
    hl2.Position = [0.5960 0.1610 0.2512 0.1937];
    set(ax2,'XColor','b','YColor','b','Box','off');
    set(gca,'fontsize',15,'fontname','Times');
    % title
    title('CCDF')
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
    hold on
    plot(DoE.Y,DoE.Y,'oc')
    hold off
    axis([y_min,y_max,y_min,y_max])
    grid on  
    xlabel('$\hat{\mathcal{M}}(\mathbf{x})$','FontSize',15);
    ylabel('$\mathcal{M}(\mathbf{x})$','FontSize',15);
end
%% 
%% -----------------------------------------------------------------------------------------------
% ============================ Picture #5 =================================
% -------------------------------------------------------------------------------------------------------------- 
if PDC.PTO(5)==1   
    ymin = ProSys.ymin;
    ymax = ProSys.ymax;
    BW = (ymax-ymin)/SBM.NofInterval;
    yy = ymin:BW:ymax;
    
    fg(P_Nov)=figure(P_Nov);P_Nov = P_Nov+1;
    set(gcf,'Units','centimeters','Position',[-11,-10,22,20]); 
    
    nn = 1:ceil((size(ALSMTimeHis.CDF,1)-1)/3):size(ALSMTimeHis.CDF,1);
    nn = [nn,size(ALSMTimeHis.CDF,1)]; nn = unique(nn);
%     nn = 1:5:size(ALSMTimeHis.std_FY,1)
%     nn = [1:10:70,70]
% nn = size(ALSMTimeHis.std_FY,1)-5:size(ALSMTimeHis.std_FY,1)
    for ii = 1:length(nn)
        CDF_i = ALSMTimeHis.CDF(nn(ii),:);
        CDF1_i = ALSMTimeHis.CDF1(nn(ii),:);
        CDF3_i = ALSMTimeHis.CDF3(nn(ii),:);
        % CDF
        subplot(length(nn),2,2*ii-1)
        %  3 sigma confidence interval
        maxY = CDF3_i;
        minY =  CDF1_i;
        maxY(find(maxY<=0))=1/SBM.NofSamples;
        minY(find(minY<=0))=1/SBM.NofSamples;
        yFill = [maxY, fliplr(minY)];
        xFill = [yy, fliplr(yy)];  
        ll1 = fill(xFill, yFill, [0.83 0.82 0.78]);
        ll1.LineStyle = 'None';
        % Exact (MCS) and predicted CDF
        hold on
        plot(yy,CDF_i,'-.')
        hold on
        plot(yy,SBM.CDF_ture)
        ax1 = gca;
        ax1.YScale = 'log';
        ax1.XLim = [ymin,ymax];
        grid on
        ylabel(['iter=',num2str(nn(ii))])
        
        % CCDF
        subplot(length(nn),2,2*ii)
        %  3 sigma confidence interval
        maxY = 1-CDF1_i;
        minY = 1-CDF3_i;
        maxY(find(maxY<=0))=1/SBM.NofSamples;
        minY(find(minY<=0))=1/SBM.NofSamples;
        yFill = [maxY, fliplr(minY)];
        xFill = [yy, fliplr(yy)];  % fliplr函数：左右翻转
        ll1 = fill(xFill, yFill, [0.83 0.82 0.78]);
        ll1.LineStyle = 'None';
        % Exact (MCS) and predicted CDF
        hold on
        semilogy(yy,(1-CDF_i),'-.','LineWidth',2)
        hold on
        semilogy(yy,SBM.CCDF_ture,'LineWidth',2)
        ax1 = gca;
        ax1.YScale = 'log';
%        ax1.YScale = 'Linear';
        ax1.XLim = [ymin,ymax];
        grid on        
    end
%% -----------------------------------------------------------------------------------------------
% ============================ Picture #6 =================================
% --------------------------------------------------------------------------------------------------------------    
if PDC.PTO(6)==1
    fg(P_Nov)=figure(P_Nov);P_Nov = P_Nov+1;
    subplot(2,1,1)
    set(gcf,'Units','centimeters','Position',[0,0,12,14]);
    nn = numel(ALSMTimeHis.errorW_y);
    plot(1:nn,ALSMTimeHis.errorW_y,'-.')
    if nn==1  nn = 2;  end
    %% 
    xlim([1,nn])
    xlabel('iter')
    ylabel('${{\varepsilon}_W}$')
    grid on
    
    subplot(2,1,2)
    set(gcf,'Units','centimeters','Position',[0,0,15,16]);
    nn = numel(ALSMTimeHis.W_y  );
    plot(1:nn,ALSMTimeHis.W_y ,'-.')
    xlim([1,nn])
    xlabel('iter')
    ylabel('${\hat{\varepsilon}_W}$')
    grid on       
end
%% -----------------------------------------------------------------------------------------------
% ============================ Picture #7 =================================
% --------------------------------------------------------------------------------------------------------------   
if PDC.PTO(7)==1
    figure
    error_wy = abs(SBM.CDF_ture-SBM.CDF)./min([SBM.CDF_ture;1-SBM.CDF_ture]);
    error_wy(isnan(error_wy))=0;error_wy(isinf(error_wy))=0;
    plot(yy,error_wy)
    xlabel('$y$')
    ylabel('${{\varepsilon}_{W_y}}$')
    grid on
end

end
% % save picture
% figpath = 'C:\software\matlab2016\bin\dsmj\可靠度\ALRM\Disscusion_Correlation\AL';
% figpath = cd;
% TestExample = 'eg4_100intervals'; 
% figname1 = [TestExample,'-','SampleDistri and Contour'];
% figname2 = [TestExample,'-','Original Model and metamodel'];
% figname3 = [TestExample,'-','CDF and CCDF'];
% figname4 = [TestExample,'-','Gpredicct and Gture'];
% figname5 = [TestExample,'-','iteration of CDFandCCDF'];
% saveas(fg(1),[figpath,'\',figname1,'.png']); 
% saveas(fg(2),[figpath,'\',figname2,'.png']);
% saveas(fg(3),[figpath,'\',figname3,'.png']);saveas(fg(3),[figpath,'\',figname3,'.fig']);
% saveas(fg(4),[figpath,'\',figname4,'.png']);saveas(fg(4),[figpath,'\',figname4,'.fig']);
% saveas(fg(5),[figpath,'\',figname5,'.png']);
return