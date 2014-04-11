% Perform non-regression tests on two-sample tests in SnPM. 
% Check that results obtained using the batch version are identical to the 
% results computed manually (using the interactive GUI).
%_______________________________________________________________________
% Copyright (C) 2013 The University of Warwick
% Id: test_twoSample.m  SnPM13 2013/10/12
% Camille Maumet
classdef test_twoSample < generic_test_snpm
    properties
    end
    
    methods (TestMethodSetup)
        function create_basis_matlabbatch(testCase)
            testCase.numBetas = 2;
            
            for k = 1:3
                testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampT.scans1{k,1} = ...
                     fullfile(testCase.testDataDir, ['test_data_', num2str(k, '%02.0f'), '.nii']);
            end
            for k = 18:20
                testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampT.scans2{k-17,1} = ...
                     fullfile(testCase.testDataDir, ['test_data_gr2_', num2str(k, '%02.0f'), '.nii']);
            end
        end
    end
    
    methods (Test)
        % No covariate, no variance smoothing
        function test_twosample_1(testCase)
            testCase.testName = 'twosample_1';
        end
        
        % No covariate, no variance smoothing (1 vs. group)
        function test_twosample_unbalanced(testCase)
            testCase.testName = 'twosample_unbalanced';
            
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampT.scans1{end+1,1} = ...
                     fullfile(testCase.testDataDir, ['test_data_04.nii']);
            
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampT.scans2(end,:) = [];
        end
        
        % No covariate, no variance smoothing and cluster stat
        function test_twosample_cluster(testCase)
            testCase.testName = 'twosample_cluster';
            
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampT.ST.ST_later = -1;
            
            % Test FDR, FWE et uncorrected T thresh as well
            additional_results(testCase);
            additional_cluster_results(testCase);
            additional_cluster_mass_results(testCase);
        end

        % No covariate, no variance smoothing and cluster stat with
        % pre-defined height threshold
        function test_twosample_cluster_predefined(testCase)
            testCase.testName = 'twosample_cluster_predefined';
            
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampT.ST.ST_U = 0.1;
            
            additional_predefined_cluster_results(testCase);
        end
        
        
        % No covariate, no variance smoothing and cluster stat with
        % pre-defined height threshold at default value
        function test_twosample_cluster_predef_stat(testCase)
            testCase.testName = 'twosample_cluster_predef_stat';
            
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampT.ST.ST_U = 2.03;
            
            additional_predefined_cluster_results(testCase);
        end

%         % With covariate with wrong number of values
%         function test_twosample_wrong_cov(testCase)
%             testCase.testName = 'twosample_wrong_cov';
%             testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampT.cov.c = [1 5 2 21];
%             testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampT.cov.cname = 'Wrong covar';
%         end
        
        % With 1 covariate
        function test_twosample_cov(testCase)
            testCase.compaWithSpm = false; % Not matching SPM orthog?
            
            testCase.numBetas = 3;
            testCase.testName = 'twosample_cov';
            
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampT.cov.c = [1 5 2 21 0 3];
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampT.cov.cname = 'Age';
        end
        
        % With 3 covariates
        function test_twosample_cov3(testCase)
            testCase.compaWithSpm = false; % Not matching SPM orthog?
            
            testCase.testName = 'twosample_cov3';
            testCase.numBetas = 5;
            
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampT.cov(1).c = [1 5 2 21 0 3];
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampT.cov(1).cname = 'Age';
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampT.cov(2).c = [1 3 5 7 3 5];
%             u = [1 5 2 21 0 3 6 14 8 5];
%             x=u(:).'/norm(u);
%             yz=null(x).';
%             xyz=[x;yz]*norm(u);
%             testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampT.cov(2).c = xyz(2,:);
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampT.cov(2).cname = 'Height';
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampT.cov(3).c = [-1 0.5 0.6 -0.1 2 1];
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampT.cov(3).cname = 'Width';
        end

        % With variance smoothing
        function test_twosample_var(testCase)
            testCase.testName = 'twosample_var';
            
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampT.vFWHM = [6 6 6];
        end

        % With approximate test
        function test_twosample_approx(testCase)
            testCase.testName = 'twosample_approx';
            
            rand('seed',200);            
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampT.nPerm = 15;
        end
        
        % Global normalisation, normalisation: Proportional scaling scaled 
        % to default value (50)
        function test_twosample_propscaling(testCase)
            testCase.testName = 'twosample_propscaling';
            
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampT.globalm.glonorm = 2;
        end
        
        % Global normalisation, normalisation: Proportional scaling scale 
        % to user-defined value
        function test_twosample_propscaling_to_user(testCase)
            testCase.testName = 'twosample_propscaling_to_user';
            
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampT.globalm.glonorm = 2;
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampT.globalm.gmsca.gmsca_yes.gmscv = 145;
        end

        % Global normalisation: overall grand mean scaling to 145
        function test_twosample_grandmean_145(testCase)
            testCase.testName = 'twosample_grandmean_145';
            
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampT.globalm.gmsca.gmsca_yes.gmscv = 145;
        end
        
        % Global calculation -> "User". 
        % Global normalisation -> "Proportional" 
        % Overall grand mean scaling: No
        function test_twosample_proportional_global_user(testCase)
            
            testCase.testName = 'twosample_proportional_global_user';
            
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampT.globalc.g_user.global_uval = [1 3 2 2 3 1];
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampT.globalm.gmsca.gmsca_no = 1;
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampT.globalm.glonorm = 2;
        end
        
        % Global normalisation: overall grand mean scaling to 50
        function test_twosample_grandmean_50(testCase)
            testCase.testName = 'twosample_grandmean_50';
            
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampT.globalm.gmsca.gmsca_yes.gmscv = 50;
        end
        
        % Global normalisation, normalisation: ANCOVA
        function test_twosample_ancova(testCase)
            testCase.compaWithSpm = false; % Not matching SPM orthog?
            testCase.testName = 'twosample_ancova';
            
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampT.globalm.glonorm = 3;
        end
    end
    
    methods % Start with last method...
        function complete_batch(testCase)
            % Find the result directory for the batch execution and the
            % corresponding result directory computed manually using the
            % original spm2-like interface
            testCase.batchResDir = fullfile(testCase.parentDataDir, 'results', 'batch', testCase.testName);
            testCase.interResDir = fullfile(spm_str_manip(testCase.batchResDir,'hh'), 'GT', testCase.testName);
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampT.dir = {testCase.batchResDir};
        end
        
        function create_spm_batch(testCase)
            
            factoDesign = testCase.spmBatch{1}.spm.stats.factorial_design;
            
            factoDesign.des.t2.scans1 = factoDesign.scans1;
            
            % Very important as otherwise scans are mixed in SPM
            if size(factoDesign.des.t2.scans1, 2) > 1
                factoDesign.des.t2.scans1 = factoDesign.des.t2.scans1';
            end
            factoDesign.des.t2.scans2 = factoDesign.scans2;
            if size(factoDesign.des.t2.scans2, 2) > 1
                factoDesign.des.t2.scans2 = factoDesign.des.t2.scans2';
            end
                
            factoDesign = rmfield(factoDesign, 'scans1');
            factoDesign = rmfield(factoDesign, 'scans2');
            if isfield(factoDesign, 'ST')
                factoDesign = rmfield(factoDesign, 'ST');
            end
            
            % % Equal variance
            %factoDesign.des.t2.variance = 0;
            
            testCase.spmBatch{1}.spm.stats.factorial_design = factoDesign;
        end
    end
end

