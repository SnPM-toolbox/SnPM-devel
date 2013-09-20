% Compare results (t-map and positive effects filtered map T>3) of
% interactive and batch multi-subject regression
classdef test_multisubsimpleregression < matlab.unittest.TestCase & generic_test_snpm
    properties
    end
    
    methods (TestMethodSetup)
        function create_basis_matlabbatch(testCase)
            %matlabbatch{1}.cfg_snpm.Design.Corr.dir = {batchResDir};
            testCase.matlabbatch{1}.cfg_snpm.Design.Corr.P = {
                                                         fullfile(testCase.testDataDir, 'su_control01', 'cn_sess1', 'con_0001.img,1')
                                                         fullfile(testCase.testDataDir, 'su_control02', 'cn_sess1', 'con_0001.img,1')
                                                         fullfile(testCase.testDataDir, 'su_control03', 'cn_sess1', 'con_0001.img,1')
                                                         fullfile(testCase.testDataDir, 'su_control04', 'cn_sess1', 'con_0001.img,1')
                                                         fullfile(testCase.testDataDir, 'su_control05', 'cn_sess1', 'con_0001.img,1')
                                                         };
            testCase.matlabbatch{1}.cfg_snpm.Design.Corr.CovInt = [1 3 5 0 2];
        end
    end
    
    
    
    methods (Test)
        % No covariate, no variance smoothing
        function test_multisubsimpleregression_1(testCase)
            testCase.testName = 'multisubsimpleregression_1';
        end

        % With variance smoothing
        function test_multisubsimpleregression_var(testCase)
            testCase.testName = 'multisubsimpleregression_var';
            
            testCase.matlabbatch{1}.cfg_snpm.Design.Corr.vFWHM = [6 6 6];
        end

        % With approximate test
        function test_multisubsimpleregression_approx(testCase)
            testCase.testName = 'multisubsimpleregression_approx';
            
            rand('seed',200);
            
            testCase.matlabbatch{1}.cfg_snpm.Design.Corr.P(end+1:end+8) = {
                 fullfile(testCase.testDataDir, 'su_control06', 'cn_sess1', 'con_0001.img,1')
                 fullfile(testCase.testDataDir, 'su_control07', 'cn_sess1', 'con_0001.img,1')
                 fullfile(testCase.testDataDir, 'su_control08', 'cn_sess1', 'con_0001.img,1')
                 fullfile(testCase.testDataDir, 'su_control09', 'cn_sess1', 'con_0001.img,1')
                 fullfile(testCase.testDataDir, 'su_control10', 'cn_sess1', 'con_0001.img,1')
                 fullfile(testCase.testDataDir, 'su_control11', 'cn_sess1', 'con_0001.img,1')
                 fullfile(testCase.testDataDir, 'su_control12', 'cn_sess1', 'con_0001.img,1')
                 fullfile(testCase.testDataDir, 'su_control13', 'cn_sess1', 'con_0001.img,1')
                 };
            testCase.matlabbatch{1}.cfg_snpm.Design.Corr.CovInt = [1 3 5 0 2 6 7 2 1 -1 2 3 1];
            testCase.matlabbatch{1}.cfg_snpm.Design.Corr.nPerm = 100;
        end
    end
    
    methods (TestMethodTeardown)
        function complete_batch(testCase)
            % Find the result directory for the batch execution and the
            % corresponding result directory computed manually using the
            % original spm2-like interface
            testCase.batchResDir = fullfile(testCase.parentDataDir, 'results', 'batch', testCase.testName);
            testCase.interResDir = fullfile(spm_str_manip(testCase.batchResDir,'hh'), 'interactive', testCase.testName);
            testCase.matlabbatch{1}.cfg_snpm.Design.Corr.dir = {testCase.batchResDir};
        end
    end
end
