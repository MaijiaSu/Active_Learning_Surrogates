         % Calculation of three kinds of models
         G1 = SBM.G_predict-2*SBM.Gmse;
         G2 = SBM.G_predict;
         G3 = SBM.G_predict+2*SBM.Gmse;  
         % Calculation of CDF of three models
         ymin = ProSys.ymin;
         ymax = ProSys.ymax;
         BW = (ymax-ymin)/100;
         Edges = [-inf,ymin:BW:ymax,+inf];   
         CDF1 = histcounts(G1,Edges,'Normalization','CDF');CDF1(end)=[];
         CDF2 = histcounts(G2,Edges,'Normalization','CDF');CDF2(end)=[];
         CDF3 = histcounts(G3,Edges,'Normalization','CDF');CDF3(end)=[];
         
%          figure
%          xx = Edges(2:end-1);
%          maxY = CDF1;minY = CDF3;
%          yFill = [maxY, fliplr(minY)];
%          xFill = [xx, fliplr(xx)];
%          fill(xFill, yFill, [0.83 0.82 0.78]);
%          set(gca,'YScale','Log')
         %% another calculation method for error measure
          yy = ymin:BW:ymax;
          for ii = 1:length(yy)
              temp11(:,ii) = normcdf((yy(ii)-SBM.G_predict)./SBM.Gmse);
          end
          mu_FY = mean(temp11);
          std_FY = sqrt(mean((1-temp11).*temp11));
         
          figure
          xx = ymin:BW:ymax;
          maxY2 = mu_FY-1*std_FY;minY2 =mu_FY+1*std_FY;
          yFill = [maxY2, fliplr(minY2)];
          xFill = [xx, fliplr(xx)];
          fill(xFill, yFill, [0.83 0.82 0.78]);
          hold on
          plot(ymin:BW:ymax,mu_FY,'r')
          set(gca,'YScale','Log')
          
          hold on
          plot(Edges(2:end-1),CDF1,'r--','LineWidth',2)
          hold on
          plot(Edges(2:end-1),CDF3,'r--','LineWidth',2)
          
         %%   comparsion of error measure
         
%          figure
%          plot(Edges(2:end-1),(CDF1-CDF3),'--g')
% %          set(gca,'YScale','Log')
%          figure
%          plot(ymin:BW:ymax,std_FY,'--g')
         
        % method #1
          CCDF2=1-CDF2;CCDF2(find(CCDF2<0))=0;
          w_y = abs((CDF3-CDF1))./min([CDF2;CCDF2]);
          w_y(isnan(w_y))=0;w_y(isinf(w_y))=0;

         
        % method #2  
          figure
          w_y2 = std_FY./min([mu_FY;1-mu_FY]);
          hold on 
          
%           figure
%           plot(Edges(2:end-1),w_y)
%           figure
%           plot(ymin:BW:ymax,w_y2)

          figure
          plot(Edges(2:end-1),abs((CDF3-CDF1)/max(abs((CDF3-CDF1)))),'r')
          hold on
          plot(ymin:BW:ymax,std_FY/max(std_FY),'g')
          %% 
          xc = [0:0.02:1];
          xcc = 0.01:0.02:1;
          for ii = 1:length(yy)
              pdfN(ii,:) = histcounts(temp11(:,ii),[0:0.02:1],'Normalization','count');
          end
          
          figure
          plot(yy,pdfN(:,1),'r')
          hold on
         plot(yy,pdfN(:,end),'--g') 
%           for ii = 1:length(yy)
% %               plot([0:0.02:1],pdfN(ii,5))
%               plot3(ones(1,101)*,[0:0.02:1],pdfN(ii,:))
%               hold on
%           end
% hist( temp(:,85))
%% MCS to calculate the distribution of FY
close all
 Edges = [-inf,ymin:BW:ymax,+inf]; 
for ii = 1:1000
    GG1 = normrnd(SBM.G_predict,10*SBM.Gmse);
    CDFii(ii,:) = histcounts(GG1,Edges,'Normalization','CDF');
end
CCDFii = 1-CDFii;

% CDF
figure
for ii = 1:1000
    plot(Edges(2:end), CDFii(ii,:))
    hold on
end  
       set(gca,'YScale','Log')
hold on
plot(Edges(2:end-1),CDF1,'r--')
hold on
plot(Edges(2:end-1),CDF3,'r--')

% figure
% hist(CDFii(:,10),20)
CDF_ave = mean(CDFii);
CDF_std = std(CDFii);

figure
plot(Edges(2:end), CDF_ave)
hold on 
plot(ymin:BW:ymax,mu_FY,'r')

figure
plot(Edges(2:end),CDF_std/max(CDF_std),'b')
hold on 
plot(ymin:BW:ymax,std_FY/max(std_FY),'r')
figure
plot(Edges(2:end-1),CDF1-CDF3,'r--')



% CCDF
% figure
% for ii = 1:1000
%     plot(Edges(2:end), CCDFii(ii,:))
%     hold on
% end  
%        set(gca,'YScale','Log')  
% hold on
% plot(Edges(2:end-1),1-CDF1,'r--','LineWidth',2)
% hold on
% plot(Edges(2:end-1),1-CDF3,'r--','LineWidth',2)

% G1 = SBM.G_predict-2*SBM.Gmse;
%          G2 = SBM.G_predict;

