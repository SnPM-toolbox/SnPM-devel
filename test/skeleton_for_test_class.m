% Compare results (t-map and positive effects filtered map T>3) of
% interactive and batch XXX
classdef skeleton_for_test_class < generic_test_snpm
    properties
    end
    
    methods (TestMethodSetup)
        function create_basis_matlabbatch(testCase)
            % Fill the design part in the batch
        end
    end
    
    methods (Test)
        % No covariate, no variance smoothing
        function test_testclassname_1(testCase)
            % Nominal test
            testCase.testName = 'testclassname_1';
        end

        % With variance smoothing
        function test_testclassname_var(testCase)
            testCase.testName = 'testclassname_var';
            
        end

        % With approximate test
        function test_testclassname_approx(testCase)
            testCase.testName = 'testclassname_approx';
            
            rand('seed',200);
        end
    end
    
    methods (TestMethodTeardown)
        function complete_batch(testCase)
            % Find the result directory for the batch execution and the
            % corresponding result directory computed manually using the
            % original spm2-like interface
            testCase.batchResDir = fullfile(testCase.parentDataDir, 'results', 'batch', testCase.testName);
            testCase.interResDir = fullfile(spm_str_manip(testCase.batchResDir,'hh'), 'interactive', testCase.testName);
            testCase.matlabbatch{1}.spm.tools.snpm.des.DESIGN_NAME.dir = {testCase.batchResDir};
        end
    end
end
