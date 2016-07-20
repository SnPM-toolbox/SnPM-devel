% Perform non-regression tests on 1 subject, 2 conditions tests in SnPM. 
% Check that results obtained using the batch version are identical to the 
% results computed manually (using the interactive GUI).
%_______________________________________________________________________
% Copyright (C) 2013 The University of Warwick
% Id: test_onesub_twocondrepl.m  SnPM13 2013/10/12
% Camille Maumet
classdef test_onesub_twocondrepl < generic_test_snpm
    properties
    end
    
    methods (TestMethodSetup)
        function create_basis_matlabbatch(testCase)
            testCase.compaWithSpm = false;
            
            % Fill the design part in the batch
            for i = 1:12
                testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampTss.P{i} = ...
                     fullfile(testCase.testDataDir, ['test_data_', num2str(i, '%02.0f') '.nii,1']);
            end
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampTss.condidx = [1 2 1 2 1 2 1 2 1 2 1 2];
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampTss.Tss_repc = 6;
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampTss.TwosampTss_Block = 4;
        end
    end
    
    methods (Test)
        % No covariate, no variance smoothing
        function test_onesub_twocondrepl_1(testCase)
            % Nominal test
            testCase.testName = 'onesub_twocondrepl_1';
        end
        
        % No covariate, no variance smoothing
        function test_onesub_twocondrepl_1_other_design(testCase)
            % Nominal test with alternative design
            testCase.testName = 'onesub_twocondrepl_1_other_design';
            % Discard last two volumes
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampTss.P(end) = [];
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampTss.P(end) = [];
            % Change conditions
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampTss.condidx = [1 2 1 2 1 2 1 2 1 2];
            % Change number of replications
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampTss.Tss_repc = 5;
            % Change exchangeability block size
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampTss.TwosampTss_Block = 10;
        end

        % With variance smoothing
        function test_onesub_twocondrepl_var(testCase)
            testCase.testName = 'onesub_twocondrepl_var';
            
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampTss.vFWHM = [12 12 12];
        end
        
        % With approx (does not check equal as do not exist in SnPM8)
        function test_onesub_twocondrepl_approx(testCase)
            testCase.checks = false;
            testCase.testName = 'onesub_twocondrepl_approx';
            
            try
                % Syntax for newest Matlab versions
                rng(200);
            catch
                % Old syntax
                rand('seed',200);
            end
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampTss.nPerm = 15;
        end
    end
    
    methods (TestMethodTeardown)
        function complete_batch(testCase)
            % Find the result directory for the batch execution and the
            % corresponding result directory computed manually using the
            % original spm2-like interface
            testCase.batchResDir = fullfile(testCase.parentDataDir, 'results', 'batch', testCase.testName);
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampTss.dir = {testCase.batchResDir};
        end
    end
end