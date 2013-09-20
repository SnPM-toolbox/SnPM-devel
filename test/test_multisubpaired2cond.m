% Compare results (t-map and positive effects filtered map T>3) of
% interactive and batch multi-subject paired with 2 conditions
classdef test_multisubpaired2cond < matlab.unittest.TestCase
    properties
        testDataDir;
        batchResDir;
        parentDataDir;
        matlabbatch;
        interResDir;
        testName;
    end
    
    methods (TestMethodSetup)
        function setGlobals(testCase)
            global TEST;
            TEST = true;
            
            testCase.parentDataDir = fullfile('/Users/cmaumet/Documents/Data/snpm_test_data');
            testCase.testDataDir = fullfile(testCase.parentDataDir, 'data');
        end
        
        function create_basis_matlabbatch(testCase)
            testCase.matlabbatch{1}.cfg_snpm.Design.PairT.fsubject(1).scans = {
                  fullfile(testCase.testDataDir, 'su_control01', 'cn_sess1', 'con_0001.img,1')
                  fullfile(testCase.testDataDir, 'su_control01', 'cn_sess2', 'con_0001.img,1')
                  };
            testCase.matlabbatch{1}.cfg_snpm.Design.PairT.fsubject(1).scindex = [1 2];
            testCase.matlabbatch{1}.cfg_snpm.Design.PairT.fsubject(2).scans = {
                  fullfile(testCase.testDataDir, 'su_control02', 'cn_sess1', 'con_0001.img,1')
                  fullfile(testCase.testDataDir, 'su_control02', 'cn_sess2', 'con_0001.img,1')
                  };
            testCase.matlabbatch{1}.cfg_snpm.Design.PairT.fsubject(2).scindex = [1 2];
            testCase.matlabbatch{1}.cfg_snpm.Design.PairT.fsubject(3).scans = {
                  fullfile(testCase.testDataDir, 'su_control03', 'cn_sess1', 'con_0001.img,1')
                  fullfile(testCase.testDataDir, 'su_control03', 'cn_sess2', 'con_0001.img,1')
                  };      
            testCase.matlabbatch{1}.cfg_snpm.Design.PairT.fsubject(3).scindex = [1 2];
            testCase.matlabbatch{1}.cfg_snpm.Design.PairT.fsubject(4).scans = {
                  fullfile(testCase.testDataDir, 'su_control04', 'cn_sess1', 'con_0001.img,1')
                  fullfile(testCase.testDataDir, 'su_control04', 'cn_sess2', 'con_0001.img,1')
                  };      
            testCase.matlabbatch{1}.cfg_snpm.Design.PairT.fsubject(4).scindex = [1 2];
            testCase.matlabbatch{1}.cfg_snpm.Design.PairT.fsubject(5).scans = {
                  fullfile(testCase.testDataDir, 'su_control05', 'cn_sess1', 'con_0001.img,1')
                  fullfile(testCase.testDataDir, 'su_control05', 'cn_sess2', 'con_0001.img,1')
                  };      
            testCase.matlabbatch{1}.cfg_snpm.Design.PairT.fsubject(5).scindex = [1 2];

            % Compute
            testCase.matlabbatch{2}.cfg_snpm.snpm_bch_cp.snpmcfg(1) = cfg_dep;
            testCase.matlabbatch{2}.cfg_snpm.snpm_bch_cp.snpmcfg(1).tname = 'SnPMcfg.mat configuration file';
            testCase.matlabbatch{2}.cfg_snpm.snpm_bch_cp.snpmcfg(1).tgt_spec = {};
            testCase.matlabbatch{2}.cfg_snpm.snpm_bch_cp.snpmcfg(1).sname = 'MultiSub: One Sample T test on diffs/contrasts: SnPMcfg.mat configuration file';
            testCase.matlabbatch{2}.cfg_snpm.snpm_bch_cp.snpmcfg(1).src_exbranch = substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
            testCase.matlabbatch{2}.cfg_snpm.snpm_bch_cp.snpmcfg(1).src_output = substruct('.','SnPMcfg');

            % Results   
            testCase.matlabbatch{3}.cfg_snpm.Infer.SnPMmat(1) = cfg_dep;
            testCase.matlabbatch{3}.cfg_snpm.Infer.SnPMmat(1).tname = 'SnPM.mat results file';
            testCase.matlabbatch{3}.cfg_snpm.Infer.SnPMmat(1).tgt_spec = {};
            testCase.matlabbatch{3}.cfg_snpm.Infer.SnPMmat(1).sname = 'Compute: SnPM.mat results file';
            testCase.matlabbatch{3}.cfg_snpm.Infer.SnPMmat(1).src_exbranch = substruct('.','val', '{}',{2}, '.','val', '{}',{1});
            testCase.matlabbatch{3}.cfg_snpm.Infer.SnPMmat(1).src_output = substruct('.','SnPM');
            testCase.matlabbatch{3}.cfg_snpm.Infer.Thr.Vox.VoxSig.Pth = 0.10;
            testCase.matlabbatch{3}.cfg_snpm.Infer.Tsign = 1;
            testCase.matlabbatch{3}.cfg_snpm.Infer.WriteFiltImg.name = 'SnPM_filtered_10none.nii';
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
        
            testCase.matlabbatch{1}.cfg_snpm.Design.PairT.fsubject(2).scans = testCase.matlabbatch{1}.cfg_snpm.Design.PairT.fsubject(2).scans(2:-1:1);
            testCase.matlabbatch{1}.cfg_snpm.Design.PairT.fsubject(2).scindex = testCase.matlabbatch{1}.cfg_snpm.Design.PairT.fsubject(2).scindex(2:-1:1);

            testCase.matlabbatch{1}.cfg_snpm.Design.PairT.fsubject(4).scans = testCase.matlabbatch{1}.cfg_snpm.Design.PairT.fsubject(4).scans(2:-1:1);
            testCase.matlabbatch{1}.cfg_snpm.Design.PairT.fsubject(4).scindex = testCase.matlabbatch{1}.cfg_snpm.Design.PairT.fsubject(4).scindex(2:-1:1);
        end

        % With variance smoothing
        function test_multisubpaired2cond_var(testCase)
            testCase.testName = 'multisubpaired2cond_var';
            
            testCase.matlabbatch{1}.cfg_snpm.Design.PairT.vFWHM = [6 6 6];
        end

        % With approximate test
        function test_multisubpaired2cond_approx(testCase)
            testCase.testName = 'multisubpaired2cond_approx';
            rand('seed',200);
            
            testCase.matlabbatch{1}.cfg_snpm.Design.PairT.nPerm = 24;
        end

    end
    methods (TestMethodTeardown) % Start with last method...
        function assert_check(testCase)
            % Compare t-maps and filtered maps
            
            % Maps obtained with the batch execution
            batch_tmap = spm_select('FPList', testCase.batchResDir, '^snpmT\+\.img');
            batch_filtmap = spm_select('FPList', testCase.batchResDir, '^SnPM_filtered_10none.*\.nii');

            % Maps obtained with the interactive execution
            inter_tmap = spm_select('FPList', testCase.interResDir, '^snpmT\+\.img');
            inter_filtmap = spm_select('FPList', testCase.interResDir, '^SnPMt_filtered_10none\.img');

            testCase.verifyEqual(spm_read_vols(spm_vol(batch_tmap)), spm_read_vols(spm_vol(inter_tmap)), 'AbsTol', 10^-10)
            testCase.verifyEqual(spm_read_vols(spm_vol(batch_filtmap)), spm_read_vols(spm_vol(inter_filtmap)), 'AbsTol', 10^-10)
            
            clear global TEST;
        end
        
        function complete_batch_and_run(testCase)
            % Find the result directory for the batch execution and the
            % corresponding result directory computed manually using the
            % original spm2-like interface
            testCase.batchResDir = fullfile(testCase.parentDataDir, 'results', 'batch', testCase.testName);
            testCase.interResDir = fullfile(spm_str_manip(testCase.batchResDir,'hh'), 'interactive', testCase.testName);
            testCase.matlabbatch{1}.cfg_snpm.Design.PairT.dir = {testCase.batchResDir};
            
            spm_jobman('run', testCase.matlabbatch);
        end
    end
end