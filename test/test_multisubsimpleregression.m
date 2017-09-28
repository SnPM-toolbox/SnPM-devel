% Perform non-regression tests on multi-subject regressions in SnPM. 
% Check that results obtained using the batch version are identical to the 
% results computed manually (using the interactive GUI).
%_______________________________________________________________________
% Copyright (C) 2013 The University of Warwick
% Id: test_multisubsimpleregression.m  SnPM13 2013/10/12
% Camille Maumet
classdef test_multisubsimpleregression < generic_test_snpm
    properties
    end
    
    methods (TestMethodSetup)
        function create_basis_matlabbatch(testCase)
            testCase.compaWithSpm = false;
            
            nSubjects = 4;
            for i = 1:nSubjects
                testCase.matlabbatch{1}.spm.tools.snpm.des.Corr.P{i,1} = ...
                     fullfile(testCase.testDataDir, ['test_data_', num2str(i, '%02.0f'), '.nii,1']);
            end
            testCase.matlabbatch{1}.spm.tools.snpm.des.Corr.CovInt = [1 3 5 0];
        end
    end
    
    
    
    methods (Test)
        % No covariate, no variance smoothing
        function test_multisubsimpleregression_1(testCase)
            testCase.testName = 'multisubsimpleregression_1';
        end

            % With 1 covariate
        function test_multisubsimpleregression_cov(testCase)
            % This was introduced in SnPM13 and can therefore not be tested
            % for non-regression with SnPM8            
            testCase.checks = false;
            
            testCase.testName = 'multisubsimpleregression_cov';
            
            testCase.matlabbatch{1}.spm.tools.snpm.des.Corr.cov.c = [1 5 2 21];
            testCase.matlabbatch{1}.spm.tools.snpm.des.Corr.cov.cname = 'age';
        end

        % With 3 covariates
        function test_multisubsimpleregression_cov3(testCase)
            % This was introduced in SnPM13 and can therefore not be tested
            % for non-regression with SnPM8        
            % TODO: we should create this test data now based on SnPM13?
            testCase.checks = false;
            
            testCase.testName = 'multisubsimpleregression_cov3';
            
            testCase.matlabbatch{1}.spm.tools.snpm.des.Corr.cov(1).c = [1 1 2 3];
            testCase.matlabbatch{1}.spm.tools.snpm.des.Corr.cov(1).cname = 'age';
            testCase.matlabbatch{1}.spm.tools.snpm.des.Corr.cov(2).c = [0 21 15 18];
            testCase.matlabbatch{1}.spm.tools.snpm.des.Corr.cov(2).cname = 'height';
            testCase.matlabbatch{1}.spm.tools.snpm.des.Corr.cov(3).c = [-1 -0.5 -1 1];
            testCase.matlabbatch{1}.spm.tools.snpm.des.Corr.cov(3).cname = 'width';
        end
        
        % With variance smoothing
        function test_multisubsimpleregression_var(testCase)
            testCase.testName = 'multisubsimpleregression_var';
            
            testCase.matlabbatch{1}.spm.tools.snpm.des.Corr.vFWHM = [8 8 8];
        end

        % With approximate test
        function test_multisubsimpleregression_approx(testCase)
            testCase.testName = 'multisubsimpleregression_approx';
            
            try
                % Syntax for newest Matlab versions
                rng(200);
            catch
                % Old syntax
                rand('seed',200);
            end
            testCase.matlabbatch{1}.spm.tools.snpm.des.Corr.nPerm = 14;
        end
    end
    
    methods (TestMethodTeardown)
        function complete_batch(testCase)
            % Find the result directory for the batch execution and the
            % corresponding result directory computed manually using the
            % original spm2-like interface
            testCase.batchResDir = fullfile(testCase.parentDataDir, 'results', 'batch', testCase.testName);
            testCase.interResDir = fullfile(spm_str_manip(testCase.batchResDir,'hh'), 'GT', testCase.testName);
            testCase.matlabbatch{1}.spm.tools.snpm.des.Corr.dir = {testCase.batchResDir};
        end
    end
end
