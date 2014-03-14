% Perform non-regression tests on within-subject ANOVA in SnPM. 
% Check that results obtained using the batch version are identical to the 
% results computed manually (using the interactive GUI).
%_______________________________________________________________________
% Copyright (C) 2013 The University of Warwick
% Id: test_multisub_withinsubanova.m  SnPM13 2013/10/12
% Camille Maumet
classdef test_multisub_withinsubanova < generic_test_snpm
    properties
    end
    
    methods (TestMethodSetup)
        function set_stat(testCase)
            testCase.stattype = 'f';
        end
        
        function create_basis_matlabbatch(testCase)
            testCase.compaWithSpm = false;
            
            nSubjects = 5;
            nScansPerSub = 2;
            
            % Fill the design part in the batch
            for i = 1:nSubjects
                for j = 1:nScansPerSub
                    testCase.matlabbatch{1}.spm.tools.snpm.des.ANOVAwithinS.fsubject(i).scans{j} = ...
                     fullfile(testCase.testDataDir, ['test_data_', num2str((i-1)*nScansPerSub+j, '%02.0f'), '.nii,1']);
                end
            end
        end
    end
    
    methods (Test)
        % No covariate, no variance smoothing
        function test_multisub_withinsubanova_1(testCase)
            % Nominal test
            testCase.testName = 'multisub_withinsubanova_1';
        end

        % With variance smoothing
        function test_multisub_withinsubanova_var(testCase)
            testCase.testName = 'multisub_withinsubanova_var';
            
            testCase.matlabbatch{1}.spm.tools.snpm.des.ANOVAwithinS.vFWHM = [8 8 8];
        end

        % With approximate test
        function test_multisub_withinsubanova_approx(testCase)
            testCase.testName = 'multisub_withinsubanova_approx';
            
            rand('seed',200);            
            testCase.matlabbatch{1}.spm.tools.snpm.des.ANOVAwithinS.nPerm = 13;
        end
    end
    
    methods (TestMethodTeardown)
        function complete_batch(testCase)
            % Find the result directory for the batch execution and the
            % corresponding result directory computed manually using the
            % original spm2-like interface
            testCase.batchResDir = fullfile(testCase.parentDataDir, 'results', 'batch', testCase.testName);
            testCase.matlabbatch{1}.spm.tools.snpm.des.ANOVAwithinS.dir = {testCase.batchResDir};
        end
    end
end