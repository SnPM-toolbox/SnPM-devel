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
            testCase.matlabbatch{1}.spm.tools.snpm.des.Corr1S.P = {
                 fullfile(testCase.testDataDir, 'PET_motor', 's8np01160em01R.img,1')
                 fullfile(testCase.testDataDir, 'PET_motor', 's8np01160em02R.img,1')
                 fullfile(testCase.testDataDir, 'PET_motor', 's8np01160em03R.img,1')
                 fullfile(testCase.testDataDir, 'PET_motor', 's8np01160em04R.img,1')
                 fullfile(testCase.testDataDir, 'PET_motor', 's8np01160em05R.img,1')
                 fullfile(testCase.testDataDir, 'PET_motor', 's8np01160em06R.img,1')
                 fullfile(testCase.testDataDir, 'PET_motor', 's8np01160em07R.img,1')
                 fullfile(testCase.testDataDir, 'PET_motor', 's8np01160em08R.img,1')
                 fullfile(testCase.testDataDir, 'PET_motor', 's8np01160em09R.img,1')
                 fullfile(testCase.testDataDir, 'PET_motor', 's8np01160em10R.img,1')
                 fullfile(testCase.testDataDir, 'PET_motor', 's8np01160em11R.img,1')
                 fullfile(testCase.testDataDir, 'PET_motor', 's8np01160em12R.img,1')
                 };
            testCase.matlabbatch{1}.spm.tools.snpm.des.Corr1S.CovInt = [1 2 0 -1 5 4 3 3 1 -2 0 1.5];
            testCase.matlabbatch{1}.spm.tools.snpm.des.Corr1S.xblock = 2;
        end
    end
    
    methods (Test)
        % No covariate, no variance smoothing
        function test_onesub_regression_1(testCase)
            % Nominal test
            testCase.testName = 'onesub_regression_1';
        end

        % With variance smoothing
        function test_onesub_regression_var(testCase)
            testCase.testName = 'onesub_regression_var';
            
            testCase.matlabbatch{1}.spm.tools.snpm.des.Corr1S.vFWHM = [8.5 8.5 8.5];
        end

        % With approximate test
        function test_onesub_regression_approx(testCase)
            testCase.testName = 'onesub_regression_approx';
            
            rand('seed',200);
            
            testCase.matlabbatch{1}.spm.tools.snpm.des.Corr1S.nPerm = 50;
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