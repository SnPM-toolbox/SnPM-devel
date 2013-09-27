% Compare results (t-map and positive effects filtered map T>3) of
% interactive and batch one-sample t-tests
classdef generic_test_snpm < matlab.unittest.TestCase
    properties
        testDataDir;
        batchResDir;
        parentDataDir;
        matlabbatch;
        spmBatch;
        spmDir;
        interResDir;
        testName;
        stattype;
        compaWithSpm;
    end
    
    methods (TestMethodSetup)
        function setGlobals(testCase)
            global TEST;
            TEST = true;

            testCase.parentDataDir = fullfile('/Users/cmaumet/Documents/Data/snpm_test_data');
            testCase.testDataDir = fullfile(testCase.parentDataDir, 'data');
            
            % By default t-test
            testCase.stattype = 't';
            
            % By default perform comparison of beta maps with SPM results
            testCase.compaWithSpm = true;
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
            
            if testCase.compaWithSpm
                % SPM batch
                testCase.spmBatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep;
                testCase.spmBatch{2}.spm.stats.fmri_est.spmmat(1).tname = 'Select SPM.mat';
                testCase.spmBatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(1).name = 'filter';
                testCase.spmBatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(1).value = 'mat';
                testCase.spmBatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(2).name = 'strtype';
                testCase.spmBatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(2).value = 'e';
                testCase.spmBatch{2}.spm.stats.fmri_est.spmmat(1).sname = 'Factorial design specification: SPM.mat File';
                testCase.spmBatch{2}.spm.stats.fmri_est.spmmat(1).src_exbranch = substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
                testCase.spmBatch{2}.spm.stats.fmri_est.spmmat(1).src_output = substruct('.','spmmat');
            end
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
            
            if testCase.compaWithSpm
                spm_beta = cellstr(spm_select('FPList', testCase.spmDir, '^beta_00\d\d\.hdr'));
            end

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
                    
                    if testCase.compaWithSpm
                        % The beta files must be equal with the one obtained by
                        % SPM (absolute tolerance lowered to 10^-4)

                        % Find corresponding beta in SPM (not necessairily in the same order)
                        corresIndex = 1;
                        minDiff = Inf;
                        for j = 1:numel(spm_beta)
                            currDiff = nansum(nansum(nansum(abs(spm_read_vols(spm_vol(batch_beta{i}))-spm_read_vols(spm_vol(spm_beta{j}))))));
                            if (currDiff < minDiff)
                                corresIndex = j;
                                minDiff = currDiff;
                            end
                        end

                        testCase.verifyEqual(spm_read_vols(spm_vol(batch_beta{i})), spm_read_vols(spm_vol(spm_beta{corresIndex})), 'AbsTol', 10^-1)
                    end
                    
                    testCase.verifyEqual(spm_read_vols(spm_vol(batch_beta{i})), spm_read_vols(spm_vol(inter_beta{i})), 'AbsTol', 10^-10)
                end
            end          
            testCase.verifyEqual(spm_read_vols(spm_vol(batch_filtmap)), spm_read_vols(spm_vol(inter_filtmap)), 'AbsTol', 10^-10)
            
            clear global TEST;
        end
        
        function complete_batch_and_run(testCase)
            testCase.complete_batch();
            spm_jobman('run', testCase.matlabbatch);
            
            if testCase.compaWithSpm
                designName = fieldnames(testCase.matlabbatch{1}.cfg_snpm.Design);
                myBatch = testCase.matlabbatch{1}.cfg_snpm.Design.(designName{1});
                testCase.spmDir = strrep(myBatch.dir{1}, 'batch', 'spm');
                myBatch.dir = {testCase.spmDir};
                testCase.spmBatch{1}.spm.stats.factorial_design = myBatch;

                testCase.create_spm_batch();

                if ~exist(testCase.spmDir, 'dir')
                    mkdir(testCase.spmDir);
                else
                    % Remove SPM.mat to avoid interactive window asking if
                    % model can be overwritten
                    delete(fullfile(testCase.spmDir, 'SPM.mat'));
                end
                % If grand mean scaling then we should calculate mean (otherwise error?)
                if (  isfield(testCase.spmBatch{1}.spm.stats.factorial_design, 'globalm') &&...
                      isfield(testCase.spmBatch{1}.spm.stats.factorial_design.globalm, 'gmsca') &&...
                      isfield(testCase.spmBatch{1}.spm.stats.factorial_design.globalm.gmsca, 'gmsca_yes')) &&...
                   ( ~isfield(testCase.spmBatch{1}.spm.stats.factorial_design, 'globalc') ||...
                      isfield(testCase.spmBatch{1}.spm.stats.factorial_design.globalc, 'g_omit'))
                        testCase.spmBatch{1}.spm.stats.factorial_design.globalc.g_mean = 1;
                        if isfield(testCase.spmBatch{1}.spm.stats.factorial_design.globalc, 'g_omit')
                            testCase.spmBatch{1}.spm.stats.factorial_design.globalc = ...
                                rmfield(testCase.spmBatch{1}.spm.stats.factorial_design.globalc, 'g_omit');
                        end
                end
                spm_jobman('run', testCase.spmBatch);
            end
        end 
    end
    
    methods
        function complete_batch(testCase)
        end 
        
    end
end