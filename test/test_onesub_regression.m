% Perform non-regression tests on one subject regression in SnPM. 
% Check that results obtained using the batch version are identical to the 
% results computed manually (using the interactive GUI).
%_______________________________________________________________________
% Copyright (C) 2013 The University of Warwick
% Id: test_onesub_regression.m  SnPM13 2013/10/12
% Camille Maumet
classdef test_onesub_regression < generic_test_snpm
    properties
    end
    
    methods (TestMethodSetup)
        function create_basis_matlabbatch(testCase)
            testCase.compaWithSpm = false;
            
            % Fill the design part in the batch
            nScans = 10;
            for i = 1:nScans
                testCase.matlabbatch{1}.spm.tools.snpm.des.Corr1S.P{i,1} = ...
                 fullfile(testCase.testDataDir, ['test_data_', num2str(i, '%02.0f'),'.nii,1']);
            end
            testCase.matlabbatch{1}.spm.tools.snpm.des.Corr1S.CovInt = [1 2 0 -1 5 4 3 3 1 -2];
            testCase.matlabbatch{1}.spm.tools.snpm.des.Corr1S.xblock = 2;
        end
    end
    
    methods (Test)
        % No covariate, no variance smoothing
        function test_onesub_regression_1(testCase)
            % Nominal test
            testCase.testName = 'onesub_regression_1';
        end

        % With 1 covariate
        function test_onesub_regression_cov(testCase)
            % This was introduced in SnPM13 and can therefore not be tested
            % for non-regression with SnPM8            
            testCase.checks = false;
            
            testCase.testName = 'onesub_regression_cov';
            
            testCase.matlabbatch{1}.spm.tools.snpm.des.Corr1S.cov.c = [1 5 2 21 0 4 5 7 8 3];
            testCase.matlabbatch{1}.spm.tools.snpm.des.Corr1S.cov.cname = 'age';
        end

        % With 3 covariates
        function test_onesub_regression_cov3(testCase)
            % This was introduced in SnPM13 and can therefore not be tested
            % for non-regression with SnPM8            
            testCase.checks = false;
            
            testCase.testName = 'onesub_regression_cov3';
            
            testCase.matlabbatch{1}.spm.tools.snpm.des.Corr1S.cov(1).c = [1 1 2 3 1 2 4 3 1 2];
            testCase.matlabbatch{1}.spm.tools.snpm.des.Corr1S.cov(1).cname = 'age';
            testCase.matlabbatch{1}.spm.tools.snpm.des.Corr1S.cov(2).c = [0 21 15 18 3 0 14 5 5 2];
            testCase.matlabbatch{1}.spm.tools.snpm.des.Corr1S.cov(2).cname = 'height';
            testCase.matlabbatch{1}.spm.tools.snpm.des.Corr1S.cov(3).c = [-1 -0.5 -1 1 0 1 1.5 0.5 -0.5 1];
            testCase.matlabbatch{1}.spm.tools.snpm.des.Corr1S.cov(3).cname = 'width';
        end
        
        % With variance smoothing
        function test_onesub_regression_var(testCase)
            testCase.testName = 'onesub_regression_var';
            
            testCase.matlabbatch{1}.spm.tools.snpm.des.Corr1S.vFWHM = [8.5 8.5 8.5];
        end

        % With approximate test
        function test_onesub_regression_approx(testCase)
            testCase.testName = 'onesub_regression_approx';
            
            try
                % Syntax for newest Matlab versions
                rng(200);
            catch
                % Old syntax
                rand('seed',200);
            end
            
            testCase.matlabbatch{1}.spm.tools.snpm.des.Corr1S.nPerm = 25;
        end
    end
    
    methods (TestMethodTeardown)
        function complete_batch(testCase)
            % Find the result directory for the batch execution and the
            % corresponding result directory computed manually using the
            % original spm2-like interface
            testCase.batchResDir = fullfile(testCase.parentDataDir, 'results', 'batch', testCase.testName);
            testCase.interResDir = fullfile(spm_str_manip(testCase.batchResDir,'hh'), 'GT', testCase.testName);
            testCase.matlabbatch{1}.spm.tools.snpm.des.Corr1S.dir = {testCase.batchResDir};
        end
    end
end