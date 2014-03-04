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
            
            % Fill the design part in the batch
            for gr1 = 1:5
                % Warning: error in spm12b if files are defined by rows
                % instead of colmuns.
                testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampPairT.scans1.fsubject(gr1).scans = {
                    fullfile(testCase.testDataDir, ['su_control' num2str(gr1, '%02d')], 'cn_sess1', 'con_0001.img,1');...
                    fullfile(testCase.testDataDir, ['su_control' num2str(gr1, '%02d')], 'cn_sess3', 'con_0001.img,1')...
                    };
                testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampPairT.scans1.fsubject(gr1).scindex = [1 2];
            end
            for gr2 = 6:10
                % Warning: error in spm12b if files are defined by rows
                % instead of colmuns.
                testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampPairT.scans2.fsubject(gr2-5).scans = {
                    fullfile(testCase.testDataDir, ['su_control' num2str(gr2, '%02d')], 'cn_sess1', 'con_0001.img,1');...
                    fullfile(testCase.testDataDir, ['su_control' num2str(gr2, '%02d')], 'cn_sess3', 'con_0001.img,1')...
                    };
                testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampPairT.scans2.fsubject(gr2-5).scindex = [1 2];
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
            
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampPairT.scans1.fsubject(3).scans = ...
                testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampPairT.scans1.fsubject(3).scans(2:-1:1);
            
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampPairT.scans1.fsubject(5).scans = ...
                testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampPairT.scans1.fsubject(5).scans(2:-1:1);            
            
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampPairT.scans1.fsubject(3).scindex = [2 1];
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampPairT.scans1.fsubject(5).scindex = [2 1];
        end    
        
         % With 1 covariate
        function test_twosample_twocond_cov(testCase)
            testCase.testName = 'twosample_twocond_cov';
            
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampPairT.cov.c = [1 5 2 21 0 5 4 8 7 1 1 5 2 21 0 5 4 8 7 1];
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampPairT.cov.cname = 'Age';
        end

        % With 3 covariates
        function test_twosample_twocond_cov3(testCase)
            testCase.testName = 'twosample_twocond_cov3';
            
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampPairT.cov(1).c = [1 1 2 3 1 5 4 6 3 1 1 2 3 1 5 4 6 3 1 18];
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampPairT.cov(1).cname = 'Age';
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampPairT.cov(2).c = [0 21 15 18 3 4 22 1 5 4 21 15 18 3 4 22 1 5 4 4];
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampPairT.cov(2).cname = 'Height';
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampPairT.cov(3).c = [-1 -0.5 -1 1 0 2 0.1 1 -1 2 -0.5 -1 1 0 2 0.1 1 -1 2 4];
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampPairT.cov(3).cname = 'Width';
        end

        % With variance smoothing
        function test_twosample_twocond_var(testCase)
            testCase.testName = 'twosample_twocond_var';
            
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampPairT.vFWHM = [9 9 9];
        end

        % With approximate test
        function test_twosample_twocond_approx(testCase)
            testCase.testName = 'twosample_twocond_approx';
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampPairT.nPerm = 100;
            
            rand('seed',200);
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
