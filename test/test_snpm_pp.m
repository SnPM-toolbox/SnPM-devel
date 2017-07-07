% Perform tests on error and warnings thrown by snpm_pp
%_______________________________________________________________________
% Copyright (C) 2016 The University of Warwick
% Camille Maumet
classdef test_snpm_pp < matlab.unittest.TestCase
    properties
        testDataDir;
        batchResDir;
        parentDataDir;
        matlabbatch;
        testName;
        SnPMrefVersion;
        warningId;
    end
    
    methods (TestMethodSetup)
        
        function setGlobals(testCase)           
            % Random number generator should not be initialised with a
            % shuffled seed
            global SnPMdefs
            SnPMdefs.shuffle_seed = false;
            
            % Run the tests in command line mode (no display)
            global defaults;
            defaults.cmdline = true;
            
            % Disable warning on very small number of permutations
            warning('off','SnPM:VeryFewPermsCoarseExactPValues')

            snpm_test_config;
            cd(spm_str_manip(which('snpm_test_config'), 'h'))
            global testDataDir;
            global SnPMrefVersion;
            testCase.SnPMrefVersion = SnPMrefVersion;

            if isempty(testDataDir)
              error('SnPM:NotTestDataDir', 'Test data directory not set, please update snpm_test_config');
            end
            
            testCase.parentDataDir = testDataDir;
            testCase.testDataDir = fullfile(testDataDir, 'data');
            
            testCase.warningId = '';
        end
    end
    
    methods (Test)
        % Test error: No covariate, no variance smoothing and cluster stat
        % and huge SnPM_ST file leading to Matlab error
        function test_onesample_cluster_err_bigfile(testCase)
            testCase.testName = 'onesample_cluster_err_bigfile';
            batch_dir = fullfile(testCase.parentDataDir, 'results', 'batch');
            testCase.batchResDir = fullfile(batch_dir, testCase.testName);
            
            % Assumes that the onesample_cluster test has been computed          
            cluster_dir = fullfile(batch_dir, 'onesample_cluster');
            if ~isdir(cluster_dir)
                error('SnPM:test:MissingOneSampleClusterTest', ...
                    'One-sample cluster test is missing: Test of the snpm_pp cannot be computed')
            end

            copyfile(fullfile(cluster_dir),testCase.batchResDir)            
            
            % Create very big SnPM.mat
            snpmmatfile = fullfile(testCase.batchResDir, 'SnPM_ST.mat');
            load(snpmmatfile);
            init_ST = SnPM_ST;
            
            SnPM_ST = init_ST;
            SnPM_ST = [SnPM_ST; SnPM_ST];
%                 save(fullfile(testCase.batchResDir, 'SnPM_ST.mat'), 'SnPM_ST')

            fid   = fopen(snpmmatfile,'a');
            fwrite(fid,SnPM_ST,'double');
            fclose(fid);    
            
%             MATLAB:load:unableToReadMatFile
            
            testCase.matlabbatch{1}.spm.tools.snpm.inference.SnPMmat = {fullfile(testCase.batchResDir, 'SnPM.mat')};            
            testCase.matlabbatch{1}.spm.tools.snpm.inference.Thr.Clus.ClusSize.CFth = 4;
            testCase.matlabbatch{1}.spm.tools.snpm.inference.Thr.Clus.ClusSize.ClusSig.PthC = 0.1;
            testCase.matlabbatch{1}.spm.tools.snpm.inference.WriteFiltImg.name = 'SnPMt_filtered_clus_4_unc_p10.nii';
            testCase.matlabbatch{1}.spm.tools.snpm.inference.export.export_no = 0;
            
            
            testCase.warningId = 'SnPM:SnPMSTFileNotLOaded';
        end
    end
    
    methods (TestMethodTeardown) % Start with last method...
        function run_batch_and_test(testCase)
            verifyError(testCase, @()(spm_jobman('run', testCase.matlabbatch)),'matlabbatch:run:jobfailederr')
            % Check manually for warning as is overwritten by error and 
            % therefore cannot be checked directly with verifyWarning
            [~, msgid] = lastwarn; 
            if ~msgid == testCase.warningId
                error(['SnPM:test:MissingWarwning:' testCase.warningId], ['Missing warning: ' testCase.warningId]);
            end
        end
    end
end