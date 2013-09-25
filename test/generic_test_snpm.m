% Compare results (t-map and positive effects filtered map T>3) of
% interactive and batch one-sample t-tests
classdef generic_test_snpm < matlab.unittest.TestCase
    properties
        testDataDir;
        batchResDir;
        parentDataDir;
        matlabbatch;
        interResDir;
        testName;
        stattype;
    end
    
    methods (TestMethodSetup)
        function setGlobals(testCase)
            global TEST;
            TEST = true;

            testCase.parentDataDir = fullfile('/Users/cmaumet/Documents/Data/snpm_test_data');
            testCase.testDataDir = fullfile(testCase.parentDataDir, 'data');
            
            % By default t-test
            testCase.stattype = 't';
        end
        
        function update_basis_matlabbatch(testCase)
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
    
    methods (TestMethodTeardown) % Start with last method...
        
        function assert_check(testCase)
            % Compare t-maps and filtered maps
            
            % Maps obtained with the batch execution
            if strcmp(testCase.stattype, 't')
                statMapRegexp = '^snpmT\+\.img';
            else
                statMapRegexp = '^snpmF\.img';
            end
                
            batch_tmap = spm_select('FPList', testCase.batchResDir, statMapRegexp);
            batch_beta = cellstr(spm_select('FPList', testCase.batchResDir, '^beta_00\d\d\.hdr'));
            batch_filtmap = spm_select('FPList', testCase.batchResDir, '^SnPM_filtered_10none.*\.nii');

            % Maps obtained with the interactive execution
            inter_tmap = spm_select('FPList', testCase.interResDir, statMapRegexp);
            inter_beta = cellstr(spm_select('FPList', testCase.interResDir, '^beta_00\d\d\.hdr'));
            inter_filtmap = spm_select('FPList', testCase.interResDir, '^SnPMt_filtered_10none\.img');

            testCase.verifyEqual(spm_read_vols(spm_vol(batch_tmap)), spm_read_vols(spm_vol(inter_tmap)), 'AbsTol', 10^-10)
            if numel(batch_beta) ~= numel(inter_beta)
                error(['Number of betas are not equal between batch (',...
                        num2str(numel(batch_beta)),...
                        ') and interactive mode ',...
                        num2str(numel(inter_beta))   ]);
            else
                for i = 1:numel(inter_beta)
                    testCase.verifyEqual(spm_read_vols(spm_vol(batch_beta{i})), spm_read_vols(spm_vol(inter_beta{i})), 'AbsTol', 10^-10)
                end
            end          
            testCase.verifyEqual(spm_read_vols(spm_vol(batch_filtmap)), spm_read_vols(spm_vol(inter_filtmap)), 'AbsTol', 10^-10)
            
            clear global TEST;
        end
        
        function complete_batch_and_run(testCase)
            testCase.complete_batch();
            
            spm_jobman('run', testCase.matlabbatch);
        end 
    end
    
    methods
        function complete_batch(testCase)
        end 
        
    end
end