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
        numBetas;
        inter_map;
        batch_map;
        tolerance;
        mapName;
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
            testCase.matlabbatch{2}.spm.tools.snpm.cp.snpmcfg(1) = cfg_dep;
            testCase.matlabbatch{2}.spm.tools.snpm.cp.snpmcfg(1).tname = 'SnPMcfg.mat configuration file';
            testCase.matlabbatch{2}.spm.tools.snpm.cp.snpmcfg(1).tgt_spec = {};
            testCase.matlabbatch{2}.spm.tools.snpm.cp.snpmcfg(1).sname = 'MultiSub: One Sample T test on diffs/contrasts: SnPMcfg.mat configuration file';
            testCase.matlabbatch{2}.spm.tools.snpm.cp.snpmcfg(1).src_exbranch = substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
            testCase.matlabbatch{2}.spm.tools.snpm.cp.snpmcfg(1).src_output = substruct('.','SnPMcfg');

            % Results   
            % Uncorrected voxel-wise p<0.1
            testCase.matlabbatch{3}.spm.tools.snpm.inference.SnPMmat(1) = cfg_dep;
            testCase.matlabbatch{3}.spm.tools.snpm.inference.SnPMmat(1).tname = 'SnPM.mat results file';
            testCase.matlabbatch{3}.spm.tools.snpm.inference.SnPMmat(1).tgt_spec = {};
            testCase.matlabbatch{3}.spm.tools.snpm.inference.SnPMmat(1).sname = 'Compute: SnPM.mat results file';
            testCase.matlabbatch{3}.spm.tools.snpm.inference.SnPMmat(1).src_exbranch = substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
            testCase.matlabbatch{3}.spm.tools.snpm.inference.SnPMmat(1).src_output = substruct('.','SnPM');
            testCase.matlabbatch{3}.spm.tools.snpm.inference.Thr.Vox.VoxSig.Pth = 0.10;
            testCase.matlabbatch{3}.spm.tools.snpm.inference.Tsign = 1;
            testCase.matlabbatch{3}.spm.tools.snpm.inference.WriteFiltImg.name = 'SnPM_filtered_10none.nii';
            
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
            batch_ip = cellstr(spm_select('FPList', testCase.batchResDir, '^lP.*.hdr'));
            batch_filtmap = cellstr(spm_select('FPList', testCase.batchResDir, '^SnPMt?_filtered_.*\.nii'));
            
            if testCase.compaWithSpm
                spm_beta = cellstr(spm_select('FPList', testCase.spmDir, '^beta_00\d\d\.hdr'));
            end

            % Maps obtained with the interactive execution
            inter_tmap = spm_select('FPList', testCase.interResDir, statMapRegexp);
            inter_beta = cellstr(spm_select('FPList', testCase.interResDir, '^beta_00\d\d\.hdr'));
            inter_ip = cellstr(spm_select('FPList', testCase.interResDir, '^lP.*.hdr'));
            inter_filtmap = cellstr(spm_select('FPList', testCase.interResDir, '^SnPMt_filtered_.*\.img'));
            
            % Compare t-maps
            testCase.inter_map = inter_tmap;
            testCase.batch_map = batch_tmap;
            testCase.tolerance = 10^-10;
            testCase.mapName = 'tmap';
            testCase.compare_batch_with_inter();
            
            % Compare beta maps
            testCase.inter_map = inter_beta(1:testCase.numBetas);
            testCase.batch_map = batch_beta(1:testCase.numBetas);
            testCase.tolerance = 10^-10;
            testCase.mapName = 'beta';
            testCase.compare_batch_with_inter();
            
            % Compare BATCH beta maps with SPM betas
            % The beta files must be equal with the one obtained by
            % SPM (absolute tolerance lowered to 10^-4)
            if testCase.compaWithSpm
                % Find SPM beta corresponding to each BATCH beta
                spmBetaIndices = NaN*ones(1, testCase.numBetas);
                for i = 1:testCase.numBetas
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
                    spmBetaIndices(i) = corresIndex;
                end
                
                testCase.inter_map = inter_beta(1:testCase.numBetas);
                testCase.batch_map = spm_beta(spmBetaIndices);
                testCase.tolerance = 10^-2;
                testCase.mapName = 'beta, SPM';
                testCase.compare_batch_with_inter();
            end
            
            % Compare lP maps
            testCase.inter_map = inter_ip;
            testCase.batch_map = batch_ip;
            testCase.tolerance = 10^-10;
            testCase.mapName = 'lP';
            testCase.compare_batch_with_inter();
            
            % Compare filtered maps
            testCase.inter_map = inter_filtmap;
            testCase.batch_map = batch_filtmap;
            testCase.tolerance = 10^-10;
            testCase.mapName = 'filtered map';
            testCase.compare_batch_with_inter();
            
            clear global TEST;
        end
        
        function complete_batch_and_run(testCase)
            testCase.complete_batch();
            
            if ~exist(testCase.batchResDir, 'dir')
                mkdir(testCase.batchResDir);
            else
                % Remove filtered maps, betas, lP from previous runs
                prevBetas = spm_select('FPList', testCase.batchResDir, '^beta.*');
                prevlP = spm_select('FPList', testCase.batchResDir, '^lP.*');
                prevtmap = spm_select('FPList', testCase.batchResDir, '^snpmT\+.*');
                prevtmapneg = spm_select('FPList', testCase.batchResDir, '^snpmT-.*');
                prevfilt = spm_select('FPList', testCase.batchResDir, '^SnPM_filtered_.*');
                prevsnpmmat = spm_select('FPList', testCase.batchResDir, '^SnPM.*\.mat');
                prevps = spm_select('FPList', testCase.batchResDir, '^spm.*\.ps');
                prevxyz = spm_select('FPList', testCase.batchResDir, '^XYZ\.mat');
                prevRes = spm_select('FPList', testCase.batchResDir, '^ResMS.*');
                
                filesToDelete = cellstr(strvcat(prevBetas, prevlP, prevtmap, ...
                    prevtmapneg, prevfilt, prevsnpmmat, prevps, prevxyz, prevRes));
                
                if ~isempty(filesToDelete)
                    for i = 1:numel(filesToDelete)
                        delete(filesToDelete{i});
                    end
                end
            end
            spm_jobman('run', testCase.matlabbatch);
            
            if testCase.compaWithSpm
                designName = fieldnames(testCase.matlabbatch{1}.spm.tools.snpm.des);
                myBatch = testCase.matlabbatch{1}.spm.tools.snpm.des.(designName{1});
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
        
        function additional_results(testCase)
            % Rename uncorrected p<0.1
            testCase.matlabbatch{end}.spm.tools.snpm.inference.WriteFiltImg.name = 'SnPMt_filtered_vox_unc_p10.nii';
            
            % Uncorrected voxel-wise TorF > 1.6
            testCase.matlabbatch{end+1}.spm.tools.snpm.inference.SnPMmat(1) = cfg_dep;
            testCase.matlabbatch{end}.spm.tools.snpm.inference.SnPMmat(1).tname = 'SnPM.mat results file';
            testCase.matlabbatch{end}.spm.tools.snpm.inference.SnPMmat(1).tgt_spec = {};
            testCase.matlabbatch{end}.spm.tools.snpm.inference.SnPMmat(1).sname = 'Compute: SnPM.mat results file';
            testCase.matlabbatch{end}.spm.tools.snpm.inference.SnPMmat(1).src_exbranch = substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
            testCase.matlabbatch{end}.spm.tools.snpm.inference.SnPMmat(1).src_output = substruct('.','SnPM');
            testCase.matlabbatch{end}.spm.tools.snpm.inference.Thr.Vox.VoxSig.TFth = 1.6;
            testCase.matlabbatch{end}.spm.tools.snpm.inference.WriteFiltImg.name = 'SnPMt_filtered_vox_unc_t16.nii';
            
            % FWE voxel-wise p<0.5
            testCase.matlabbatch{end+1}.spm.tools.snpm.inference.SnPMmat(1) = cfg_dep;
            testCase.matlabbatch{end}.spm.tools.snpm.inference.SnPMmat(1).tname = 'SnPM.mat results file';
            testCase.matlabbatch{end}.spm.tools.snpm.inference.SnPMmat(1).tgt_spec = {};
            testCase.matlabbatch{end}.spm.tools.snpm.inference.SnPMmat(1).sname = 'Compute: SnPM.mat results file';
            testCase.matlabbatch{end}.spm.tools.snpm.inference.SnPMmat(1).src_exbranch = substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
            testCase.matlabbatch{end}.spm.tools.snpm.inference.SnPMmat(1).src_output = substruct('.','SnPM');
            testCase.matlabbatch{end}.spm.tools.snpm.inference.Thr.Vox.VoxSig.FWEth = 0.5;
            testCase.matlabbatch{end}.spm.tools.snpm.inference.WriteFiltImg.name = 'SnPMt_filtered_vox_fwe_p50.nii'; 
            
            % FDR voxel-wise p<0.5
            testCase.matlabbatch{end+1}.spm.tools.snpm.inference.SnPMmat(1) = cfg_dep;
            testCase.matlabbatch{end}.spm.tools.snpm.inference.SnPMmat(1).tname = 'SnPM.mat results file';
            testCase.matlabbatch{end}.spm.tools.snpm.inference.SnPMmat(1).tgt_spec = {};
            testCase.matlabbatch{end}.spm.tools.snpm.inference.SnPMmat(1).sname = 'Compute: SnPM.mat results file';
            testCase.matlabbatch{end}.spm.tools.snpm.inference.SnPMmat(1).src_exbranch = substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
            testCase.matlabbatch{end}.spm.tools.snpm.inference.SnPMmat(1).src_output = substruct('.','SnPM');
            testCase.matlabbatch{end}.spm.tools.snpm.inference.Thr.Vox.VoxSig.FDRth = 0.5;
            testCase.matlabbatch{end}.spm.tools.snpm.inference.WriteFiltImg.name = 'SnPMt_filtered_vox_fdr_p50.nii';  
        end
        
        function additional_cluster_results(testCase)
            % Uncorrected (cluster forming u>4) cluster-wise p<0.1
            testCase.matlabbatch{end+1}.spm.tools.snpm.inference.SnPMmat(1) = cfg_dep;
            testCase.matlabbatch{end}.spm.tools.snpm.inference.SnPMmat(1).tname = 'SnPM.mat results file';
            testCase.matlabbatch{end}.spm.tools.snpm.inference.SnPMmat(1).tgt_spec = {};
            testCase.matlabbatch{end}.spm.tools.snpm.inference.SnPMmat(1).sname = 'Compute: SnPM.mat results file';
            testCase.matlabbatch{end}.spm.tools.snpm.inference.SnPMmat(1).src_exbranch = substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
            testCase.matlabbatch{end}.spm.tools.snpm.inference.SnPMmat(1).src_output = substruct('.','SnPM');
            testCase.matlabbatch{end}.spm.tools.snpm.inference.Thr.Clus.CFth = 4;
            testCase.matlabbatch{end}.spm.tools.snpm.inference.Thr.Clus.ClusSig.PthC = 0.1;
            testCase.matlabbatch{end}.spm.tools.snpm.inference.WriteFiltImg.name = 'SnPMt_filtered_clus_4_unc_p10.nii';
            
            % Uncorrected (cluster forming u>4) cluster-wise k>6
            testCase.matlabbatch{end+1}.spm.tools.snpm.inference.SnPMmat(1) = cfg_dep;
            testCase.matlabbatch{end}.spm.tools.snpm.inference.SnPMmat(1).tname = 'SnPM.mat results file';
            testCase.matlabbatch{end}.spm.tools.snpm.inference.SnPMmat(1).tgt_spec = {};
            testCase.matlabbatch{end}.spm.tools.snpm.inference.SnPMmat(1).sname = 'Compute: SnPM.mat results file';
            testCase.matlabbatch{end}.spm.tools.snpm.inference.SnPMmat(1).src_exbranch = substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
            testCase.matlabbatch{end}.spm.tools.snpm.inference.SnPMmat(1).src_output = substruct('.','SnPM');
            testCase.matlabbatch{end}.spm.tools.snpm.inference.Thr.Clus.CFth = 4;
            testCase.matlabbatch{end}.spm.tools.snpm.inference.Thr.Clus.ClusSig.Cth = 6;
            testCase.matlabbatch{end}.spm.tools.snpm.inference.WriteFiltImg.name = 'SnPMt_filtered_clus_4_unc_k6.nii';
            
            % FWE (cluster forming u>4) cluster-wise p<0.5
            testCase.matlabbatch{end+1}.spm.tools.snpm.inference.SnPMmat(1) = cfg_dep;
            testCase.matlabbatch{end}.spm.tools.snpm.inference.SnPMmat(1).tname = 'SnPM.mat results file';
            testCase.matlabbatch{end}.spm.tools.snpm.inference.SnPMmat(1).tgt_spec = {};
            testCase.matlabbatch{end}.spm.tools.snpm.inference.SnPMmat(1).sname = 'Compute: SnPM.mat results file';
            testCase.matlabbatch{end}.spm.tools.snpm.inference.SnPMmat(1).src_exbranch = substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
            testCase.matlabbatch{end}.spm.tools.snpm.inference.SnPMmat(1).src_output = substruct('.','SnPM');
            testCase.matlabbatch{end}.spm.tools.snpm.inference.Thr.Clus.CFth = 4;
            testCase.matlabbatch{end}.spm.tools.snpm.inference.Thr.Clus.ClusSig.FWEthC = 0.5;
            testCase.matlabbatch{end}.spm.tools.snpm.inference.WriteFiltImg.name = 'SnPMt_filtered_clus_4_fwe_p50.nii'; 
            
            % FWE (cluster forming u>5) cluster-wise p<0.5
            testCase.matlabbatch{end+1}.spm.tools.snpm.inference.SnPMmat(1) = cfg_dep;
            testCase.matlabbatch{end}.spm.tools.snpm.inference.SnPMmat(1).tname = 'SnPM.mat results file';
            testCase.matlabbatch{end}.spm.tools.snpm.inference.SnPMmat(1).tgt_spec = {};
            testCase.matlabbatch{end}.spm.tools.snpm.inference.SnPMmat(1).sname = 'Compute: SnPM.mat results file';
            testCase.matlabbatch{end}.spm.tools.snpm.inference.SnPMmat(1).src_exbranch = substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
            testCase.matlabbatch{end}.spm.tools.snpm.inference.SnPMmat(1).src_output = substruct('.','SnPM');
            testCase.matlabbatch{end}.spm.tools.snpm.inference.Thr.Clus.CFth = 5;
            testCase.matlabbatch{end}.spm.tools.snpm.inference.Thr.Clus.ClusSig.FWEthC = 0.5;
            testCase.matlabbatch{end}.spm.tools.snpm.inference.WriteFiltImg.name = 'SnPMt_filtered_clus_5_fwe_p50.nii';  
        end
        
        function additional_predifined_cluster_results(testCase)
            % Rename uncorrected p<0.1
            testCase.matlabbatch{end}.spm.tools.snpm.inference.WriteFiltImg.name = 'SnPMt_filtered_vox_unc_p10.nii';
            
            % Uncorrected (cluster forming p<0.1) cluster-wise p<0.1
            testCase.matlabbatch{end+1}.spm.tools.snpm.inference.SnPMmat(1) = cfg_dep;
            testCase.matlabbatch{end}.spm.tools.snpm.inference.SnPMmat(1).tname = 'SnPM.mat results file';
            testCase.matlabbatch{end}.spm.tools.snpm.inference.SnPMmat(1).tgt_spec = {};
            testCase.matlabbatch{end}.spm.tools.snpm.inference.SnPMmat(1).sname = 'Compute: SnPM.mat results file';
            testCase.matlabbatch{end}.spm.tools.snpm.inference.SnPMmat(1).src_exbranch = substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
            testCase.matlabbatch{end}.spm.tools.snpm.inference.SnPMmat(1).src_output = substruct('.','SnPM');
            testCase.matlabbatch{end}.spm.tools.snpm.inference.Thr.Clus.CFth = NaN;
            testCase.matlabbatch{end}.spm.tools.snpm.inference.Thr.Clus.ClusSig.PthC = 0.1;
            testCase.matlabbatch{end}.spm.tools.snpm.inference.WriteFiltImg.name = 'SnPMt_filtered_clus_p10_unc_p10.nii';
            
            % FWE (cluster forming p<0.1) cluster-wise p<0.5
            testCase.matlabbatch{end+1}.spm.tools.snpm.inference.SnPMmat(1) = cfg_dep;
            testCase.matlabbatch{end}.spm.tools.snpm.inference.SnPMmat(1).tname = 'SnPM.mat results file';
            testCase.matlabbatch{end}.spm.tools.snpm.inference.SnPMmat(1).tgt_spec = {};
            testCase.matlabbatch{end}.spm.tools.snpm.inference.SnPMmat(1).sname = 'Compute: SnPM.mat results file';
            testCase.matlabbatch{end}.spm.tools.snpm.inference.SnPMmat(1).src_exbranch = substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
            testCase.matlabbatch{end}.spm.tools.snpm.inference.SnPMmat(1).src_output = substruct('.','SnPM');
            testCase.matlabbatch{end}.spm.tools.snpm.inference.Thr.Clus.CFth = NaN;
            testCase.matlabbatch{end}.spm.tools.snpm.inference.Thr.Clus.ClusSig.FWEthC = 0.5;
            testCase.matlabbatch{end}.spm.tools.snpm.inference.WriteFiltImg.name = 'SnPMt_filtered_clus_p10_fwe_p50.nii'; 
        end
        
        function compare_batch_with_inter(testCase)
            if ~iscell(testCase.inter_map)
                testCase.inter_map = {testCase.inter_map};
            end
            if ~iscell(testCase.batch_map)
                testCase.batch_map = {testCase.batch_map};
            end
            
            if numel(testCase.inter_map) ~= numel(testCase.batch_map)
                error(['Number of ' testCase.mapName ' maps are not equal between batch (',...
                        num2str(numel(testCase.batch_map)),...
                        ') and interactive mode ',...
                        num2str(numel(testCase.inter_map))   ]);
            else
                for i = 1:numel(testCase.inter_map)
                    data1 = spm_read_vols(spm_vol(testCase.batch_map{i}));
                    data2 = spm_read_vols(spm_vol(testCase.inter_map{i}));
                    testCase.verifyEqual(data1, data2, 'AbsTol', testCase.tolerance, [testCase.batch_map{i}])
                end
            end
        end
    end
end
