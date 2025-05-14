% example
load data1
theta = [10 10]; %优化参数初值
lob = [1e-1 1e-1]; upb = [20 20]; %优化参数取值范围
% 模型建立
[dmodel, perf] = ...
dacefit(S, G, @regpoly0, @corrgauss, theta, lob, upb)
% 模型预测
X = gridsamp([0 0;100 100], 40);
[YX MSE] = predictor(X, dmodel);