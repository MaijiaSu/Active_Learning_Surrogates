% example
load data1
theta = [10 10]; %�Ż�������ֵ
lob = [1e-1 1e-1]; upb = [20 20]; %�Ż�����ȡֵ��Χ
% ģ�ͽ���
[dmodel, perf] = ...
dacefit(S, G, @regpoly0, @corrgauss, theta, lob, upb)
% ģ��Ԥ��
X = gridsamp([0 0;100 100], 40);
[YX MSE] = predictor(X, dmodel);