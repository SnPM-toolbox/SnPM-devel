% Perform non-regression tests on two-sample two conditions tests in SnPM. 
% Check that results obtained using the batch version are identical to the 
% results computed manually (using the interactive GUI).
%_______________________________________________________________________
% Copyright (C) 2013 The University of Warwick
% Id: test_twosample_twocond.m  SnPM13 2013/10/12
% Camille Maumet
classdef test_twosample_twocond < generic_test_snpm
    properties
    end
    
    methods (TestMethodSetup)
        function create_basis_matlabbatch(testCase)
            testCase.compaWithSpm = false;
            
            nSubjectsPerGroup = 3;
            nScansPerSub = 2;
            jj = 1;
            for g = 1:2
                for s = 1:nSubjectsPerGroup
                    for i = 1:nScansPerSub
                        % Warning: error in spm12b if files are defined by rows
                        % instead of colmuns.
                        testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampPairT.(['scans' num2str(g)]).fsubject(s).scans{i,1} = ...
                            fullfile(testCase.testDataDir, ['test_data_' num2str(jj, '%02d') '.nii,1']);
                        testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampPairT.(['scans' num2str(g)]).fsubject(s).scindex = [1 2];
                        jj = jj + 1;
                    end
                end
            end
        end
    end
    
    methods (Test)
        % No covariate, no variance smoothing
        function test_twosample_twocond_1(testCase)
            % Nominal test
            testCase.testName = 'twosample_twocond_1';
        end
        
        % Change the acquisition order of conditions for part of the
        % subjects
        function test_twosample_twocond_chgorder(testCase)
            testCase.testName = 'twosample_twocond_chgorder';        
            
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampPairT.scans1.fsubject(1).scindex = [2 1];
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampPairT.scans1.fsubject(2).scindex = [2 1];
        end    
        
         % With 1 covariate
        function test_twosample_twocond_cov(testCase)
            testCase.testName = 'twosample_twocond_cov';
            
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampPairT.cov.c = [1 5 2 21 0 5 4 8 7 1 1 5];
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampPairT.cov.cname = 'Age';
        end

        % With 3 covariates
        function test_twosample_twocond_cov3(testCase)
            testCase.testName = 'twosample_twocond_cov3';
            
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampPairT.cov(1).c = [1 1 2 3 1 5 4 6 3 1 1 2];
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampPairT.cov(1).cname = 'Age';
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampPairT.cov(2).c = [0 21 15 18 3 4 22 1 5 4 21 15];
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampPairT.cov(2).cname = 'Height';
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampPairT.cov(3).c = [-1 -0.5 -1 1 0 2 0.1 1 -1 2 -0.5 -1];
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampPairT.cov(3).cname = 'Width';
        end

        % With variance smoothing
        function test_twosample_twocond_var(testCase)
            testCase.testName = 'twosample_twocond_var';
            
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampPairT.vFWHM = [8.5 8.5 8.5];
        end

        % With approximate test
        function test_twosample_twocond_approx(testCase)
            testCase.testName = 'twosample_twocond_approx';
            rand('seed',200);
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampPairT.nPerm = 15;
        end
    end
    
    methods (TestMethodTeardown)
        function complete_batch(testCase)
            % Find the result directory for the batch execution and the
            % corresponding result directory computed manually using the
            % original spm2-like interface
            testCase.batchResDir = fullfile(testCase.parentDataDir, 'results', 'batch', testCase.testName);
            testCase.interResDir = fullfile(spm_str_manip(testCase.batchResDir,'hh'), 'GT', testCase.testName);
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampPairT.dir = {testCase.batchResDir};
        end
    end
end
