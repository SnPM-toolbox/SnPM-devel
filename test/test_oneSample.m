% Perform non-regression tests on one sample tests in SnPM. 
% Check that results obtained using the batch version are identical to the 
% results computed manually (using the interactive GUI).
%_______________________________________________________________________
% Copyright (C) 2013 The University of Warwick
% Id: test_oneSample.m  SnPM13 2013/10/12
% Camille Maumet
classdef test_oneSample < generic_test_snpm
    properties
    end
    
    methods (TestMethodSetup)
        function create_basis_matlabbatch(testCase)
            testCase.numBetas = 1;
            
            testCase.matlabbatch{1}.spm.tools.snpm.des.OneSampT.P = {
                     fullfile(testCase.testDataDir, 'su_control01', 'cn_sess1', 'con_0001.img,1')
                     fullfile(testCase.testDataDir, 'su_control02', 'cn_sess1', 'con_0001.img,1')
                     fullfile(testCase.testDataDir, 'su_control03', 'cn_sess1', 'con_0001.img,1')
                     fullfile(testCase.testDataDir, 'su_control04', 'cn_sess1', 'con_0001.img,1')
                     fullfile(testCase.testDataDir, 'su_control05', 'cn_sess1', 'con_0001.img,1')
                     };
        end
    end

    methods (Test)
        % No covariate, no variance smoothing, no cluster stat
        function test_onesample_1(testCase)
            testCase.testName = 'onesample_1';
            
            % Test FDR, FWE et uncorrected T thresh as well
            additional_results(testCase);
        end
        
        % No covariate, no variance smoothing and cluster stat
        function test_onesample_cluster(testCase)
            testCase.testName = 'onesample_cluster';
            
            testCase.matlabbatch{1}.spm.tools.snpm.des.OneSampT.bVolm = 1;
            testCase.matlabbatch{1}.spm.tools.snpm.des.OneSampT.ST.ST_later = -1;
            
            % Test FDR, FWE and uncorrected T thresh as well
            additional_results(testCase);
            % Test cluster size statistic
            additional_cluster_results(testCase);
            % Test cluster mass statistic
            additional_cluster_mass_results(testCase);
        end

        % No covariate, no variance smoothing and cluster stat with
        % pre-defined height threshold
        function test_onesample_cluster_predefined(testCase)
            testCase.testName = 'onesample_cluster_predefined';
            
            testCase.matlabbatch{1}.spm.tools.snpm.des.OneSampT.bVolm = 1;
            testCase.matlabbatch{1}.spm.tools.snpm.des.OneSampT.ST.ST_U = 0.1;
            
            additional_predifined_cluster_results(testCase);
        end
        
        % With 1 covariate
        function test_onesample_cov(testCase)
            testCase.testName = 'onesample_cov';
            
            testCase.matlabbatch{1}.spm.tools.snpm.des.OneSampT.cov.c = [1 5 2 21 0];
            testCase.matlabbatch{1}.spm.tools.snpm.des.OneSampT.cov.cname = 'age';
        end

        % With 3 covariates
        function test_onesample_cov3(testCase)
            testCase.testName = 'onesample_cov3';
            
            testCase.matlabbatch{1}.spm.tools.snpm.des.OneSampT.cov(1).c = [1 1 2 3 1];
            testCase.matlabbatch{1}.spm.tools.snpm.des.OneSampT.cov(1).cname = 'age';
            testCase.matlabbatch{1}.spm.tools.snpm.des.OneSampT.cov(2).c = [0 21 15 18 3];
            testCase.matlabbatch{1}.spm.tools.snpm.des.OneSampT.cov(2).cname = 'height';
            testCase.matlabbatch{1}.spm.tools.snpm.des.OneSampT.cov(3).c = [-1 -0.5 -1 1 0];
            testCase.matlabbatch{1}.spm.tools.snpm.des.OneSampT.cov(3).cname = 'width';
        end

        % With variance smoothing
        function test_onesample_var(testCase)
            testCase.testName = 'onesample_var';
            
            testCase.matlabbatch{1}.spm.tools.snpm.des.OneSampT.vFWHM = [6 6 6];
        end

        % With approximate test
        function test_onesample_approx(testCase)
            testCase.testName = 'onesample_approx';
            
            rand('seed',200);
            testCase.matlabbatch{1}.spm.tools.snpm.des.OneSampT.P(end+1:end+8) = {
                 fullfile(testCase.testDataDir, 'su_control06', 'cn_sess1', 'con_0001.img,1')
                 fullfile(testCase.testDataDir, 'su_control07', 'cn_sess1', 'con_0001.img,1')
                 fullfile(testCase.testDataDir, 'su_control08', 'cn_sess1', 'con_0001.img,1')
                 fullfile(testCase.testDataDir, 'su_control09', 'cn_sess1', 'con_0001.img,1')
                 fullfile(testCase.testDataDir, 'su_control10', 'cn_sess1', 'con_0001.img,1')
                 fullfile(testCase.testDataDir, 'su_control11', 'cn_sess1', 'con_0001.img,1')
                 fullfile(testCase.testDataDir, 'su_control12', 'cn_sess1', 'con_0001.img,1')
                 fullfile(testCase.testDataDir, 'su_control13', 'cn_sess1', 'con_0001.img,1')
                 };
            testCase.matlabbatch{1}.spm.tools.snpm.des.OneSampT.nPerm = 100;
        end
        
        % Global normalisation, normalisation: Proportional scaling scaled 
        % to default value (50)
        function test_onesample_propscaling(testCase)
            testCase.testName = 'onesample_propscaling';
            
            testCase.matlabbatch{1}.spm.tools.snpm.des.OneSampT.globalm.glonorm = 2;
        end
        
        % Global normalisation, normalisation: Proportional scaling scale 
        % to user-defined value
        function test_onesample_propscaling_to_user(testCase)
            testCase.testName = 'onesample_propscaling_to_user';
            
            testCase.matlabbatch{1}.spm.tools.snpm.des.OneSampT.globalm.glonorm = 2;
            testCase.matlabbatch{1}.spm.tools.snpm.des.OneSampT.globalm.gmsca.gmsca_yes.gmscv = 145;
        end

        % Global normalisation: overall grand mean scaling to 145
        function test_onesample_grandmean_145(testCase)
            testCase.testName = 'onesample_grandmean_145';
            
            testCase.matlabbatch{1}.spm.tools.snpm.des.OneSampT.globalm.gmsca.gmsca_yes.gmscv = 145;
        end
        
        % Global normalisation: overall grand mean scaling to 50
        function test_onesample_grandmean_50(testCase)
            testCase.testName = 'onesample_grandmean_50';
            
            testCase.matlabbatch{1}.spm.tools.snpm.des.OneSampT.globalm.gmsca.gmsca_yes.gmscv = 50;
        end
        
        % Global normalisation, normalisation: ANCOVA
        function test_onesample_ancova(testCase)
            testCase.testName = 'onesample_ancova';
            
            testCase.matlabbatch{1}.spm.tools.snpm.des.OneSampT.globalm.glonorm = 3;
        end
        
        % Work slice by slice
        function test_onesample_slice(testCase)
            testCase.compaWithSpm = false;
            
            testCase.testName = 'onesample_slice';
            
            rand('seed',200);
            testCase.matlabbatch{1}.spm.tools.snpm.des.OneSampT.P(end+1:end+12) = {
                 fullfile(testCase.testDataDir, 'su_control06', 'cn_sess1', 'con_0001.img,1')
                 fullfile(testCase.testDataDir, 'su_control07', 'cn_sess1', 'con_0001.img,1')
                 fullfile(testCase.testDataDir, 'su_control08', 'cn_sess1', 'con_0001.img,1')
                 fullfile(testCase.testDataDir, 'su_control09', 'cn_sess1', 'con_0001.img,1')
                 fullfile(testCase.testDataDir, 'su_control10', 'cn_sess1', 'con_0001.img,1')
                 fullfile(testCase.testDataDir, 'su_control11', 'cn_sess1', 'con_0001.img,1')
                 fullfile(testCase.testDataDir, 'su_control12', 'cn_sess1', 'con_0001.img,1')
                 fullfile(testCase.testDataDir, 'su_control13', 'cn_sess1', 'con_0001.img,1')
                 fullfile(testCase.testDataDir, 'su_control01', 'cn_sess2', 'con_0001.img,1')
                 fullfile(testCase.testDataDir, 'su_control01', 'cn_sess3', 'con_0001.img,1')
                 fullfile(testCase.testDataDir, 'su_control01', 'cn_sess4', 'con_0001.img,1')
                 fullfile(testCase.testDataDir, 'su_control02', 'cn_sess2', 'con_0001.img,1')
                 };
            testCase.matlabbatch{1}.spm.tools.snpm.des.OneSampT.nPerm = 100;
            testCase.matlabbatch{1}.spm.tools.snpm.des.OneSampT.bVolm = 0;
        end
        
    end
    
    methods
        function complete_batch(testCase)
            % Find the result directory for the batch execution and the
            % corresponding result directory computed manually using the
            % original spm2-like interface
            testCase.batchResDir = fullfile(testCase.parentDataDir, 'results', 'batch', testCase.testName);
            testCase.interResDir = fullfile(spm_str_manip(testCase.batchResDir,'hh'), 'GT', testCase.testName);
            testCase.matlabbatch{1}.spm.tools.snpm.des.OneSampT.dir = {testCase.batchResDir};
        end
        
        function create_spm_batch(testCase)
            factoDesign = testCase.spmBatch{1}.spm.stats.factorial_design;
            
            factoDesign.des.t1.scans = factoDesign.P;
            factoDesign = rmfield(factoDesign, 'P');
            
            testCase.spmBatch{1}.spm.stats.factorial_design = factoDesign;
        end
    end
end