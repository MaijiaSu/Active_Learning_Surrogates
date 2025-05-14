addpath('ALRM')
addpath(genpath('Toolboxes'))
clc,clear,close all

%%
rng(1)
% the benchmark example in the file "BenchmarkProblemForFullDistribution"
BenchmarkID = [10 2	12	15	14	13	16	1	4	8	9	3	5	7	19	6];


for NN = BenchmarkID
    TestExample =  ['BatchProcessEg',num2str(NN)]
    BenchmarkProblemForFullDistribution

    SBM.method = 'MCS_UQ';
    SBM.iniNoS = 1e5;             % the size of the population of initial sampling-based-method
    SBM.NofInterval = 100;

    % Number of Solutions: 18
    NofSolution = 18;
    % Number of repeated runs
    NofRun = 10;

    %% Solutions
    % Metamodel = {'UQLab_Kriging','PCE','PCK'}
    % Initial_ED = {'iniRandom','iniMDS'}
    % Learning_fun = {'MoV','TwoStepLF','MoV-distance'}
    % No.  Metamodel Initial_ED	 Learning_fun
    % 1	       1	     1	          1
    % 2  	   1	     1	          2
    % 3	       1	     1	          3
    % 4	       1	     2	          1
    % 5	       1	     2	          2
    % 6	       1	     2	          3
    % 7	       2	     1	          1
    % 8	       2	     1	          2
    % 9  	   2	     1	          3
    % 10	   2	     2	          1
    % 11	   2	     2	          2
    % 12	   2	     2	          3
    % 13	   3	     1	          1
    % 14	   3	     1	          2
    % 15	   3	     1	          3
    % 16	   3	     2	          1
    % 17	   3	     2	          2
    % 18	   3	     2	          3
    SolutionMethods = [1    1	1	1	1	1	2	2	2	2	2	2	3	3	3	3	3	3
        1	1	1	2	2	2	1	1	1	2	2	2	1	1	1	2	2	2
        1	2	3	1	2	3	1	2	3	1	2	3	1	2	3	1	2	3]';

    %% Batch Process
    tic
    for ii = 1:NofSolution
        MetamodelStype = SolutionMethods(ii,1);
        AS1 = SolutionMethods(ii,2);
        AS2 = SolutionMethods(ii,3);
        SurrModelPar = []; ALSMPar = [];
        SolutionStrategies % Return Setting of 'SurrModelPar' and 'ALSMPar'

        if NN == 12||NN == 14||NN == 13
            ALSMTimeHis.MaxNoDoE = 300
        end
        if NN == 19||NN == 6
            NofRun = 10;
        end
        %     ALSMTimeHis.File_iniState = 'eg19_iniState';
        %     ALSMTimeHis.false_ture_statistics = 1;

        for jj = 1:NofRun
            disp(['--------------------','Solution-',num2str(ii),'--------------------------'])
            disp(['--------------------','RUN-',num2str(jj),'--------------------------'])
            tstart = tic;
            ALRMResult = mainALRM...
                (TestExample,SBM,SurrModelPar,ALSMPar,ALSMTimeHis,ProSys);
            Runtime = toc(tstart);

            % Save the data of interest
            Result{ii}.DoE{jj,1} =      ALRMResult.SurrModelPar.DoE;
            Result{ii}.NofDoE(jj,1) =   ALRMResult.NofDoE;
            Result{ii}.errorW_y{jj} =   ALRMResult.ALSMTimeHis.errorW_y;
            Result{ii}.w_y{jj,1} =      ALRMResult.ALSMTimeHis.w_y;
            Result{ii}.W_y{jj,1} =      ALRMResult.ALSMTimeHis.W_y;
            Result{ii}.CDF{jj,1} =      ALRMResult.ALSMTimeHis.CDF;
            Result{ii}.CDF1{jj,1} =     ALRMResult.ALSMTimeHis.CDF1;
            Result{ii}.CDF3{jj,1} =     ALRMResult.ALSMTimeHis.CDF3;
            Result{ii}.CDF_ture{jj,1} = ALRMResult.ALSMTimeHis.CDF_ture;
            Result{ii}.Runtime{jj,1} = Runtime;
        end

    end
    toc
    %%
    SaveFileName = [TestExample,'-','Result'];
    save(SaveFileName,'Result','SBM','SurrModelPar','ALSMPar','ALSMTimeHis','ProSys')
end