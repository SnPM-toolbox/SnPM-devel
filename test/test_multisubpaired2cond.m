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
            
            testCase.matlabbatch{1}.spm.tools.snpm.des.PairT.fsubject(1).scans = {
                  fullfile(testCase.testDataDir, 'su_control01', 'cn_sess1', 'con_0001.img,1')
                  fullfile(testCase.testDataDir, 'su_control01', 'cn_sess2', 'con_0001.img,1')
                  };
            testCase.matlabbatch{1}.spm.tools.snpm.des.PairT.fsubject(1).scindex = [1 2];
            testCase.matlabbatch{1}.spm.tools.snpm.des.PairT.fsubject(2).scans = {
                  fullfile(testCase.testDataDir, 'su_control02', 'cn_sess1', 'con_0001.img,1')
                  fullfile(testCase.testDataDir, 'su_control02', 'cn_sess2', 'con_0001.img,1')
                  };
            testCase.matlabbatch{1}.spm.tools.snpm.des.PairT.fsubject(2).scindex = [1 2];
            testCase.matlabbatch{1}.spm.tools.snpm.des.PairT.fsubject(3).scans = {
                  fullfile(testCase.testDataDir, 'su_control03', 'cn_sess1', 'con_0001.img,1')
                  fullfile(testCase.testDataDir, 'su_control03', 'cn_sess2', 'con_0001.img,1')
                  };      
            testCase.matlabbatch{1}.spm.tools.snpm.des.PairT.fsubject(3).scindex = [1 2];
            testCase.matlabbatch{1}.spm.tools.snpm.des.PairT.fsubject(4).scans = {
                  fullfile(testCase.testDataDir, 'su_control04', 'cn_sess1', 'con_0001.img,1')
                  fullfile(testCase.testDataDir, 'su_control04', 'cn_sess2', 'con_0001.img,1')
                  };      
            testCase.matlabbatch{1}.spm.tools.snpm.des.PairT.fsubject(4).scindex = [1 2];
            testCase.matlabbatch{1}.spm.tools.snpm.des.PairT.fsubject(5).scans = {
                  fullfile(testCase.testDataDir, 'su_control05', 'cn_sess1', 'con_0001.img,1')
                  fullfile(testCase.testDataDir, 'su_control05', 'cn_sess2', 'con_0001.img,1')
                  };      
            testCase.matlabbatch{1}.spm.tools.snpm.des.PairT.fsubject(5).scindex = [1 2];
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
        
            testCase.matlabbatch{1}.spm.tools.snpm.des.PairT.fsubject(2).scans = testCase.matlabbatch{1}.spm.tools.snpm.des.PairT.fsubject(2).scans(2:-1:1);
            testCase.matlabbatch{1}.spm.tools.snpm.des.PairT.fsubject(2).scindex = testCase.matlabbatch{1}.spm.tools.snpm.des.PairT.fsubject(2).scindex(2:-1:1);

            testCase.matlabbatch{1}.spm.tools.snpm.des.PairT.fsubject(4).scans = testCase.matlabbatch{1}.spm.tools.snpm.des.PairT.fsubject(4).scans(2:-1:1);
            testCase.matlabbatch{1}.spm.tools.snpm.des.PairT.fsubject(4).scindex = testCase.matlabbatch{1}.spm.tools.snpm.des.PairT.fsubject(4).scindex(2:-1:1);
        end

        % With variance smoothing
        function test_multisubpaired2cond_var(testCase)
            testCase.testName = 'multisubpaired2cond_var';
            
            testCase.matlabbatch{1}.spm.tools.snpm.des.PairT.vFWHM = [6 6 6];
        end

        % With approximate test
        function test_multisubpaired2cond_approx(testCase)
            testCase.testName = 'multisubpaired2cond_approx';
            rand('seed',200);
            
            testCase.matlabbatch{1}.spm.tools.snpm.des.PairT.nPerm = 24;
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