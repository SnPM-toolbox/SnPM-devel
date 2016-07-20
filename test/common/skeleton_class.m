% Skeleton for non-regression tests classes. 
%_______________________________________________________________________
% Copyright (C) 2013 The University of Warwick
% Id: skeleton_class.m  SnPM13 2013/10/12
% Camille Maumet
classdef skeleton_class %< generic_test_snpm
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
            try
                % Syntax for newest Matlab versions
                rng(200);
            catch
                % Old syntax
                rand('seed',200);
            end
        end
    end
    
    methods (TestMethodTeardown)
        function complete_batch(testCase)
            % Find the result directory for the batch execution and the
            % corresponding result directory computed manually using the
            % original spm2-like interface
            testCase.batchResDir = fullfile(testCase.parentDataDir, 'results', 'batch', testCase.testName);
            testCase.interResDir = fullfile(spm_str_manip(testCase.batchResDir,'hh'), 'GT', testCase.testName);
            testCase.matlabbatch{1}.spm.tools.snpm.des.DESIGN_NAME.dir = {testCase.batchResDir};
        end
    end
end
