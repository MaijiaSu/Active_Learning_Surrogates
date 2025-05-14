function  MoedlTrain = UQLabMetaModelTrain(DoE,G,MetaOpts)   
    MetaOpts.ExpDesign.X = DoE;
    MetaOpts.ExpDesign.Y = G;
    MoedlTrain = uq_createModel(MetaOpts);    
end

