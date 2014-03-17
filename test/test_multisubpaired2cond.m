% Perform non-regression tests on multi-sub paired 2 condition tests SnPM. 
% Check that results obtained using the batch version are identical to the 
% results computed manually (using the interactive GUI).
%_______________________________________________________________________
% Copyright (C) 2013 The University of Warwick
% Id: test_multisubpaired2cond.m  SnPM13 2013/10/12
% Camille Maumet
classdef test_multisubpaired2cond < generic_test_snpm
    properties
    end
    
    methods (TestMethodSetup)
        function create_basis_matlabbatch(testCase)
            testCase.compaWithSpm = false;
            
            nSubjects = 5;
            for i = 1:nSubjects
                for j = 1:2
                    % In SPM12, scand must be one by line (not one by
                    % colmun) otherwise error.
                    testCase.matlabbatch{1}.spm.tools.snpm.des.PairT.fsubject(i).scans{j,1} = ...
                      fullfile(testCase.testDataDir, ['test_data_', num2str((i-1)*2+j, '%02.0f'), '.nii,1']);
                    testCase.matlabbatch{1}.spm.tools.snpm.des.PairT.fsubject(i).scindex = [1 2];
                end
            end
        end
    end

    methods (Test)
        % No variance smoothing, no approximate test, all sessions in same order
        function test_multisubpaired2cond_1(testCase)
            testCase.testName = 'multisubpaired2cond_1';
        end

        % Session order inverted for two subjects
        function test_multisubpaired2cond_chgorder(testCase)
            testCase.testName = 'multisubpaired2cond_chgorder';
        
            testCase.matlabbatch{1}.spm.tools.snpm.des.PairT.fsubject(1).scindex = testCase.matlabbatch{1}.spm.tools.snpm.des.PairT.fsubject(2).scindex(2:-1:1);
            testCase.matlabbatch{1}.spm.tools.snpm.des.PairT.fsubject(2).scindex = testCase.matlabbatch{1}.spm.tools.snpm.des.PairT.fsubject(4).scindex(2:-1:1);
        end

        % With variance smoothing
        function test_multisubpaired2cond_var(testCase)
            testCase.testName = 'multisubpaired2cond_var';
            
            testCase.matlabbatch{1}.spm.tools.snpm.des.PairT.vFWHM = [8 8 8];
        end

        % With approximate test
        function test_multisubpaired2cond_approx(testCase)
            testCase.testName = 'multisubpaired2cond_approx';
            rand('seed',200);
            
            testCase.matlabbatch{1}.spm.tools.snpm.des.PairT.nPerm = 14;
        end

    end
    
    methods    
        function complete_batch(testCase)
            % Find the result directory for the batch execution and the
            % corresponding result directory computed manually using the
            % original spm2-like interface
            testCase.batchResDir = fullfile(testCase.parentDataDir, 'results', 'batch', testCase.testName);
            testCase.interResDir = fullfile(spm_str_manip(testCase.batchResDir,'hh'), 'GT', testCase.testName);
            testCase.matlabbatch{1}.spm.tools.snpm.des.PairT.dir = {testCase.batchResDir};
        end
    end
end