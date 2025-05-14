% 

%% 2维时，绘制登高线数据
P_Nov = 1;
close all
if ProSys.Ndim == 2 %当且仅当自变量维度为2
    xx = linspace(PDC.bound(1),PDC.bound(3),PDC.gap);
    yy = linspace(PDC.bound(2),PDC.bound(4),PDC.gap);
    [X1,X2] = meshgrid(xx,yy);
    
    X = [reshape(X1,PDC.gap*PDC.gap,1),reshape(X2,PDC.gap*PDC.gap,1)];
    
    % 预测
    YX_predict = SurrModelPar.ModelPredictor(SurrModelPar,X);  
    YX_predict = reshape(YX_predict,PDC.gap,PDC.gap);
    
    % 实际
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
        figure(P_Nov);P_Nov = P_Nov+1;
        % 1.绘制SamplePool
        pl(1) = plot(SBM.SamplePool(:,1),SBM.SamplePool(:,2),'o','MarkerSize',1);
        pl(1).Color = [0.83 0.82 0.78];

        hold on
        
        % 2. 绘制失效面
        % 真实失效面登高线绘制
        contour(X1, X2, YX_true, [0 0],'color','g','linewidth',1,'ShowText','on')
        % 绘图预测的等高线
        hold on
        contour(X1, X2, YX_predict, [0 0],'color','r','linewidth',1,'ShowText','on')
        
        % 3. 把状态误判的点标注出来
        hold on
        % 绘制弃真的点
        try
            pl(2) = plot(SBM.SamplePool(ALSMTimeHis.S2,1),SBM.SamplePool(ALSMTimeHis.S2,2),'vm', 'MarkerSize',8);
        catch
            pl(2) = plot(ProSys.muX(1),ProSys.muX(1),'vm', 'MarkerSize',0.1);
        end
        % 绘制取伪的点
        hold on
        try
            pl(3) = plot(SBM.SamplePool(ALSMTimeHis.S3,1),SBM.SamplePool(ALSMTimeHis.S3,2),'xm', 'MarkerSize',8);
        catch
            pl(3) = plot(ProSys.muX(1),ProSys.muX(1),'xm', 'MarkerSize',0.1);
        end
        
        % 4. 绘制初始实验点
         hold on
         pl(4) = plot(DoE.X(1:ALSMPar.IniDoE.N0,1),DoE.X(1:ALSMPar.IniDoE.N0,2),'ms');
         pl(4).MarkerSize = 5;
         pl(4).MarkerEdgeColor = 'y'; pl(4).MarkerFaceColor = 'y';
         
        % 5. 绘制最佳实验点
         hold on
         try
             pl(5) = plot(DoE.X(ALSMPar.IniDoE.N0+1:end,1),DoE.X(ALSMPar.IniDoE.N0+1:end,2),'bv');
         catch 
             pl(5) = plot(ProSys.muX(1),ProSys.muX(1),'bv','MarkerSize',0.1);
         end         
         pl(5).MarkerSize = 5;
         pl(5).MarkerEdgeColor = 'b'; pl(5).MarkerFaceColor = 'b';
       
         % ------------------- SDMCS Boundary------------------------------
         if strcmp(SBM.method,'SDMCS')
             rs2 = chi2inv(SBM.P_Divde_CDF(1:SBM.NoR),2); % r^2
             ff = @(x1,x2) x1.^2+x2.^2;
             Z1 = ff(X1,X2);
             ZZ = reshape(Z1,300,300);
             hold on
             [~,p6] = contour(X1,X2,ZZ,[rs2]);
             p6.LineWidth = 1;
             p6.LineStyle = '--';
             p6.ShowText = 'off';
         end
        %------------------------------------------------------------------
        
        title([text_ALRM{1},'(','NofDoE= ',num2str(size(DoE.X,1)),')'])
        % 6. 画布设置
        % 设置画布大小以及字体大小
        %图片大小 单位：cm
        width = 15;
        height = 12;
        FontSize = 12;       
        set(gcf,'Units','centimeters','Position',[6 6 width height]);  %图片大小设置
        
        % 设置坐标轴字体
        set(gca,'fontsize',FontSize,'fontname','Times');
        xlabel('\it\fontname{Times New Roman}x\rm\fontname{Times New Roman}_1')
        ylabel('\it\fontname{Times New Roman}x\rm\fontname{Times New Roman}_2')
         
        % 7. 设置图签
        % 打开网格
        grid on

        % 设置图签       
        temp_str{1} = 'SamplesPool';
        temp_str{2} = 'Ture LSF';
        temp_str{3} = 'Predicted LSF';
        temp_str{4} = 'False positive classified points';
        temp_str{5} = 'False negative classified points';
        temp_str{6} = 'Initial DoE';
        temp_str{7} = 'Added best DoE';    
        hl = legend(temp_str);
        hl.Box = 'off';
        hl.FontName = 'Times';
        hl.FontSize = 12;
        
        % 补充画PMC_DoE
        if ALSMPar.IniDoE.GenType == 3
            for i = 1:size(DoE.PMC_DoE,2)
                hold on
                iPMC = DoE.PMC_DoE{i};
                pp(1) = plot(DoE.X(iPMC,1),DoE.X(iPMC,2),'-.c');  %先画线，后标注
                % pp(2) = plot(DoE.X(iPMC(1:end),1),DoE.X(iPMC(1:end),2),'^k');
                % pp(3) = plot(DoE.X(iPMC(1),1),DoE.X(iPMC(1),2),'^m');
                text(DoE.X(iPMC(1),1),DoE.X(iPMC(1),2),['PMC-' num2str(i)])
            end
        end
end

%% -----------------------------------------------------------------------------------------------
% ============================ Picture #2 =================================
% --------------------------------------------------------------------------------------------------------------    
if ProSys.Ndim==2&&PDC.PTO(2)==1  
        % 绘制预测值图
        figure(P_Nov);P_Nov = P_Nov+1;
        mesh(X1, X2, YX_predict)
        title('预测响应面')
        hold on
        
        % 绘制真实响应值图
        figure(P_Nov);P_Nov = P_Nov+1;
        mesh(X1, X2, YX_true)
        title('实际响应面')
        hold on
        
        % 绘制预测值方差图
%         figure(P_Nov);P_Nov = P_Nov+1;
%         mesh(X1, X2, YXmse)
%         title('预测值方差')
%         hold on
end

%% -----------------------------------------------------------------------------------------------
% ============================ Picture #3 =================================
% -------------------------------------------------------------------------------------------------------------- 
if PDC.PTO(3)==1   
       figure(P_Nov);P_Nov = P_Nov+1;
       % 1.失效概率迭代进程
       subplot(3,1,1)
       N = length(ALSMTimeHis.Pf_predict);
       plot(1:N,ALSMTimeHis.Pf_predict,'o-');
       temp_str =[];
       temp_str{1} = ['Pf(','ALRM',')'];
       text = [];
       text = [temp_str{1},'=',num2str(ALSMTimeHis.Pf_predict(end))];
       
       if ALSMTimeHis.false_ture_statistics == 1
           hold on
           plot(1:N,ALSMTimeHis.Pf_ture,'g--')
           temp_str{2} = ['Pf(',SBM.method,')'];
           text = [text, ' ', temp_str{2},num2str(ALSMTimeHis.Pf_ture(end))];
       end
%        set(gca,'YScale','Log')
       hl = legend(temp_str);
       hl.Box = 'off';
       grid on
       set(gca,'fontsize',12,'fontname','Times');
       title(text)
       
        % 2.失效点相对判断数量    
        subplot(3,1,2)  
        if ALSMTimeHis.false_ture_statistics == 1
            N = length( ALSMTimeHis.Pf_Error);
            plot(1:N,1+ ALSMTimeHis.Pf_Error,'o-');
            hold on
            plot(1:N,1+ALSMTimeHis.MCP_state(:,2),'--g')  %弃真集合
            hold on
            plot(1:N,1+ALSMTimeHis.MCP_state(:,3),'-.c')  %取伪集合
            hold on
            plot(1:N,1*ones(1,N),'--r')
%             set(gca,'YScale','Log')
            
            temp_str =[];
            temp_str{1} = 'Relative Error of number of failur samples';
            temp_str{2} = 'False positive classified points';
            temp_str{3} = 'False negative classified points';
            temp_str{4} = 'Exact Prediction';
            
            hl = legend(temp_str);
            hl.Box = 'off';
            grid on
            set(gca,'fontsize',8,'fontname','Times');
        else
%             N = length(ALSMTimeHis.nf);
%             plot(1:N,ALSMTimeHis.nf,'o-');
%             grid on
%             hl = legend('Number of failure points');hl.Box = 'off';
%             set(gca,'fontsize',8,'fontname','Times');
                plot(0,0)
        end
        
        % 3.学习函数终止条件
        subplot(3,1,3)     
        if strcmp(ALSMPar.Stopcon.type,'minU')       
            N = length(ALSMTimeHis.LFindex);
            plot(1:N,ALSMTimeHis.LFindex,'--go')
            hold on 
            plot(1:N,ALSMPar.Stopcon.minU*ones(1,N),'r')
            hl = legend('U','Threshold');
            hl.Box = 'off';
            
        elseif strcmp(ALSMPar.Stopcon.type,'maxEFF')  
            N = length(ALSMTimeHis.LFindex);
            plot(1:N,ALSMTimeHis.LFindex,'--go')
            hold on
            plot(1:N,ALSMPar.Stopcon.maxEFF*ones(1,N),'-.r')
            set(gca,'YScale','Log','YGrid','on','XGrid','on')
            hl = legend('EFF','Threshold');
            hl.Box = 'off';
            
        elseif strcmp(ALSMPar.Stopcon.type,'SC_UQLib')||strcmp(ALSMPar.Stopcon.type,'Modifie_SC_UQLib')  
            N = length(ALSMTimeHis.SC1);
%             ALSMTimeHis.SC1(find(ALSMTimeHis.SC1==0)) = 1e-6;
%             ALSMTimeHis.SC2(find(ALSMTimeHis.SC2==0)) = 1e-6;
            plot(1:N,ALSMTimeHis.SC1,':go')
            hold on
            plot(1:N,ALSMTimeHis.SC2,':bX')
            hold on
            plot(1:N,ALSMPar.Stopcon.SC1*ones(1,N),'g')
            hold on
            plot(1:N,ALSMPar.Stopcon.SC2*ones(1,N),'b')
            
            set(gca,'YScale','Log','YGrid','off','XGrid','on')
            hl = legend('SC1','SC2','Threshold1','Threshold2');
            hl.Box = 'off';
            
        elseif strcmp(ALSMPar.Stopcon.type,'CompositeSC')  
            
            N = length(ALSMTimeHis.LFindex);
            % 坐标轴1
            plot(1:N,ALSMTimeHis.LFindex,'-b.')           
            ax1 = gca;
            hold on
            plot(1:N,ALSMPar.Stopcon.minU*ones(1,N),'--b')
            hl = legend('U','Threshold of U');
            hl.Box = 'off';
            set(ax1,'XColor','b','YColor','b','Box','off');
            set(gca,'fontsize',12,'fontname','Times');
            
            % 坐标轴2
            ax2 = axes('Position',get(ax1,'Position'),...
                'XAxisLocation','top',...
                'YAxisLocation','right',...
                'Color','none',...
                'XColor','g','YColor','g');
            N = length(ALSMTimeHis.SC1);
            hold on
            plot(1:N,ALSMTimeHis.SC1,'-.go')
            hold on
            plot(1:N,ALSMTimeHis.SC2,'-.gX')
            hold on
            plot(1:N,ALSMPar.Stopcon.SC1*ones(1,N),'--g')
            set(gca,'YScale','Log')
            hl = legend('SC1','SC2','Threshold of SCofUQLib');
            hl.Box = 'off';
            set(gca,'fontsize',12,'fontname','Times');
        end        
        
        % 3. 标记终止时刻
        for StopN = ALSMTimeHis.StopLearningTimes
            plot([StopN,StopN],[1e-6,1],'--r')
        end
        
        % 4. 画布设置          
        width = 18;
        height = 12;
        set(gcf,'Units','centimeters','Position',[6 6 width height]);  %图片大小设置       
end
%% figure Ture-predict 四象限对比
if PDC.PTO(4)==1&&ALSMTimeHis.false_ture_statistics == 1
    figure(P_Nov);P_Nov = P_Nov+1;
    
    if ALSMTimeHis.false_ture_statistics == 0
        tic
%         Y_true = ProSys.fname(SBM.SamplePool,ProSys.fun_par);
        load data_Ytrue
        toc
    else
        Y_true = SBM.G_ture;
    end
    
    Y_predict = SurrModelPar.ModelPredictor(SurrModelPar,SBM.SamplePool(1:size(Y_true,1),:));
    
    % 统计
    T1 = sign(Y_predict);
    T2 = sign(Y_true);
    temp1 = T1-T2;
    temp2 = T1+T2;
    
    
    index{1} = find(temp2==2);
    index{2} = find(temp1==-2);
    index{3} = find(temp2==-2);
    index{4} = find(temp1==2);
    
    state(1) = length(find(temp2==2));
    state(2) = length(find(temp1==-2));
    state(3) = length(find(temp2==-2));
    state(4) = length(find(temp1==2));
    
    x_max = max(abs(Y_predict));
    y_max = max(abs(Y_true));
    tempcolor = {'b.','rx','b.','rx'};
    for n = 1:4
        hold on
        plot(Y_predict(index{n}),Y_true(index{n}),tempcolor{n})
    end
    
%     plot(Y_predict,Y_true,'+')
%     hold on
    hold on
    plot([0,0],[-y_max,y_max],'m--')
    hold on
    plot([-x_max,x_max],[0,0],'m--')
    hold on
    plot([-x_max,x_max],[-x_max,x_max],'-.g')
    axis([-x_max,x_max,-y_max,y_max])
    hold on 
    
    N0 = ALSMPar.IniDoE.N0;
    plot(DoE.Y(1:N0),DoE.Y(1:N0),'co')
    hold on
    plot(DoE.Y(N0+1:end),DoE.Y(N0+1:end),'gs')
    
    xlabel('Y\_predict')
    ylabel('Y\_ture')
    title([text_ALRM{1},'(','NofDoE= ',num2str(size(DoE.X,1)),')',...
           '  ','NofState=',num2str(state)])
end

