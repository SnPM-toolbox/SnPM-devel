% Perform non-regression tests on between-group ANOVA in SnPM. 
% Check that results obtained using the batch version are identical to the 
% results computed manually (using the interactive GUI).
%_______________________________________________________________________
% Copyright (C) 2013 The University of Warwick
% Id: test_ANOVAbetween.m  SnPM13 2013/10/12
% Camille Maumet
classdef test_ANOVAbetween < generic_test_snpm
    properties
    end
    
    methods (TestMethodSetup)
        function create_basis_matlabbatch(testCase)
            testCase.compaWithSpm = false;
            testCase.stattype = 'f';
            
            % Fill the design part in the batch
            testCase.matlabbatch{1}.spm.tools.snpm.des.ANOVAbtwnS.group(1).scans = {};
            testCase.matlabbatch{1}.spm.tools.snpm.des.ANOVAbtwnS.group(2).scans = {};
            testCase.matlabbatch{1}.spm.tools.snpm.des.ANOVAbtwnS.group(3).scans = {};
            for gr1 = 1:2
                testCase.matlabbatch{1}.spm.tools.snpm.des.ANOVAbtwnS.group(1).scans{end+1,1} = ...
                    fullfile(testCase.testDataDir, ['test_data_' num2str(gr1, '%02d') '.nii,1']);
            end
            for gr2 = 3:4
                testCase.matlabbatch{1}.spm.tools.snpm.des.ANOVAbtwnS.group(2).scans{end+1,1} = ...
                    fullfile(testCase.testDataDir, ['test_data_' num2str(gr2, '%02d') '.nii,1']);
            end
            for gr3 = 5:6
                testCase.matlabbatch{1}.spm.tools.snpm.des.ANOVAbtwnS.group(3).scans{end+1,1} = ...
                    fullfile(testCase.testDataDir, ['test_data_' num2str(gr3, '%02d') '.nii,1']);
            end
            
            testCase.matlabbatch{1}.spm.tools.snpm.des.ANOVAbtwnS.nullHypAllZero = false;
        end
    end
    
    methods (Test)
        % Null hypothesis: all group means are equal
        function test_ANOVAbetween_1(testCase)
            % Nominal test
            testCase.testName = 'ANOVAbetween_1';
        end

        % Null hypothesis: all group means are zero
        function test_ANOVAbetween_allzero(testCase)
            testCase.testName = 'ANOVAbetween_allzero';
            
            testCase.matlabbatch{1}.spm.tools.snpm.des.ANOVAbtwnS.nullHypAllZero = true;
        end
        
        % With variance smoothing
        function test_ANOVAbetween_var(testCase)
            testCase.testName = 'ANOVAbetween_var';
            
            testCase.matlabbatch{1}.spm.tools.snpm.des.ANOVAbtwnS.nullHypAllZero = true;
            testCase.matlabbatch{1}.spm.tools.snpm.des.ANOVAbtwnS.vFWHM = [8 8 8];
        end

        % With approximate test
        function test_ANOVAbetween_approx(testCase)
            testCase.testName = 'ANOVAbetween_approx';           
            try
                % Syntax for newest Matlab versions
                rng(200);
            catch
                % Old syntax
                rand('seed',200);
            end
            testCase.matlabbatch{1}.spm.tools.snpm.des.ANOVAbtwnS.nPerm = 15;
        end
    end
    
    methods (TestMethodTeardown)
        function complete_batch(testCase)
            % Find the result directory for the batch execution and the
            % corresponding result directory computed manually using the
            % original spm2-like interface
            testCase.batchResDir = fullfile(testCase.parentDataDir, 'results', 'batch', testCase.testName);
            testCase.matlabbatch{1}.spm.tools.snpm.des.ANOVAbtwnS.dir = {testCase.batchResDir};
        end
    end
end
