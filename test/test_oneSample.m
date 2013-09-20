% Compare results (t-map and positive effects filtered map T>3) of
% interactive and batch one-sample t-tests
classdef test_oneSample < matlab.unittest.TestCase & generic_test_snpm
    properties
    end
    
    methods (TestMethodSetup)
        function create_basis_matlabbatch(testCase)
            testCase.matlabbatch{1}.cfg_snpm.Design.OneSampT.P = {
                     fullfile(testCase.testDataDir, 'su_control01', 'cn_sess1', 'con_0001.img,1')
                     fullfile(testCase.testDataDir, 'su_control02', 'cn_sess1', 'con_0001.img,1')
                     fullfile(testCase.testDataDir, 'su_control03', 'cn_sess1', 'con_0001.img,1')
                     fullfile(testCase.testDataDir, 'su_control04', 'cn_sess1', 'con_0001.img,1')
                     fullfile(testCase.testDataDir, 'su_control05', 'cn_sess1', 'con_0001.img,1')
                     };
        end
        
    end

    methods (Test)
        % No covariate, no variance smoothing
        function test_onesample_1(testCase)
            testCase.testName = 'onesample_1';
        end

        % With 1 covariate
        function test_onesample_cov(testCase)
            testCase.testName = 'onesample_cov';
            
            testCase.matlabbatch{1}.cfg_snpm.Design.OneSampT.mcov.c = [1 5 2 21 0];
            testCase.matlabbatch{1}.cfg_snpm.Design.OneSampT.mcov.cname = 'age';
        end

        % With 3 covariates
        function test_onesample_cov3(testCase)
            testCase.testName = 'onesample_cov3';
            
            testCase.matlabbatch{1}.cfg_snpm.Design.OneSampT.mcov(1).c = [1 1 2 3 1];
            testCase.matlabbatch{1}.cfg_snpm.Design.OneSampT.mcov(1).cname = 'age';
            testCase.matlabbatch{1}.cfg_snpm.Design.OneSampT.mcov(2).c = [0 21 15 18 3];
            testCase.matlabbatch{1}.cfg_snpm.Design.OneSampT.mcov(2).cname = 'height';
            testCase.matlabbatch{1}.cfg_snpm.Design.OneSampT.mcov(3).c = [-1 -0.5 -1 1 0];
            testCase.matlabbatch{1}.cfg_snpm.Design.OneSampT.mcov(3).cname = 'width';
        end

        % With variance smoothing
        function test_onesample_var(testCase)
            testCase.testName = 'onesample_var';
            
            testCase.matlabbatch{1}.cfg_snpm.Design.OneSampT.vFWHM = [6 6 6];
        end

        % With approximate test
        function test_onesample_approx(testCase)
            testCase.testName = 'onesample_approx';
            
            rand('seed',200);
            testCase.matlabbatch{1}.cfg_snpm.Design.OneSampT.P(end+1:end+8) = {
                 fullfile(testCase.testDataDir, 'su_control06', 'cn_sess1', 'con_0001.img,1')
                 fullfile(testCase.testDataDir, 'su_control07', 'cn_sess1', 'con_0001.img,1')
                 fullfile(testCase.testDataDir, 'su_control08', 'cn_sess1', 'con_0001.img,1')
                 fullfile(testCase.testDataDir, 'su_control09', 'cn_sess1', 'con_0001.img,1')
                 fullfile(testCase.testDataDir, 'su_control10', 'cn_sess1', 'con_0001.img,1')
                 fullfile(testCase.testDataDir, 'su_control11', 'cn_sess1', 'con_0001.img,1')
                 fullfile(testCase.testDataDir, 'su_control12', 'cn_sess1', 'con_0001.img,1')
                 fullfile(testCase.testDataDir, 'su_control13', 'cn_sess1', 'con_0001.img,1')
                 };
            testCase.matlabbatch{1}.cfg_snpm.Design.OneSampT.nPerm = 100;
        end

    end
    
    methods
        function complete_batch(testCase)
            % Find the result directory for the batch execution and the
            % corresponding result directory computed manually using the
            % original spm2-like interface
            testCase.batchResDir = fullfile(testCase.parentDataDir, 'results', 'batch', testCase.testName);
            testCase.interResDir = fullfile(spm_str_manip(testCase.batchResDir,'hh'), 'interactive', testCase.testName);
            testCase.matlabbatch{1}.cfg_snpm.Design.OneSampT.dir = {testCase.batchResDir};
        end
    end
end