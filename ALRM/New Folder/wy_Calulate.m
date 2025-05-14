function [w_y,CDF1,CDF2,CDF3] = wy_Calulate(ProSys,SBM,G_predict,Gstd)

ymin = ProSys.ymin;
ymax = ProSys.ymax;
BW = (ymax-ymin)/SBM.NofInterval;
yy = ymin:BW:ymax;
Edges = [-inf,ymin:BW:ymax,+inf];
G1 = G_predict-2*Gstd;
G2 = G_predict;
G3 = G_predict+2*Gstd;

CDF1 = histcounts(G1,[-inf,yy],'Normalization','cumcount')/SBM.NofSamples;
CDF2 = histcounts(G2,[-inf,yy],'Normalization','cumcount')/SBM.NofSamples;
CDF3 = histcounts(G3,[-inf,yy],'Normalization','cumcount')/SBM.NofSamples;

CCDF2=1-CDF2;
w_y = abs((CDF3-CDF1))./min([CDF2;CCDF2]);
w_y(isnan(w_y))=0;w_y(isinf(w_y))=0;

end