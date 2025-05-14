function CSModel = MyCrossval(SurrModelPar)   

MoedlTrain = SurrModelPar.MoedlTrain;
DoE = SurrModelPar.DoE;
ks = SurrModelPar.ks;
% 1. Divide the DoEs into ks CV-subsets 
n = size(DoE.Y,1);
temp=ks-mod(n,ks);
temp=[randperm(n),zeros(1,temp)];
CVSet=reshape(temp,[ks,length(temp)/ks]);

% 2. Train ks cross-validation models based on CV-subsets
for i=1:ks
    temp=CVSet(i,find(CVSet(i,:)>0));
    iDoE=DoE.X; iDoE(temp,:)=[];
    iG=DoE.Y;   iG(temp,:)=[];
    CSModel{i} = MoedlTrain(iDoE,iG);
end

return