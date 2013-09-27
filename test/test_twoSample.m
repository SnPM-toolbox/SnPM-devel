% Compare results (t-map and positive effects filtered map T>3) of
% interactive and batch two-sample t-tests
classdef test_twoSample < matlab.unittest.TestCase & generic_test_snpm
    properties
    end
    
    methods (TestMethodSetup)
        function create_basis_matlabbatch(testCase)
            testCase.matlabbatch{1}.cfg_snpm.Design.TwoSampT.scans1 = {...
                 fullfile(testCase.testDataDir, 'su_control01', 'cn_sess1', 'con_0001.img,1'),...
                 fullfile(testCase.testDataDir, 'su_control02', 'cn_sess1', 'con_0001.img,1'),...
                 fullfile(testCase.testDataDir, 'su_control03', 'cn_sess1', 'con_0001.img,1'),...
                 fullfile(testCase.testDataDir, 'su_control04', 'cn_sess1', 'con_0001.img,1'),...
                 fullfile(testCase.testDataDir, 'su_control05', 'cn_sess1', 'con_0001.img,1')};
             testCase.matlabbatch{1}.cfg_snpm.Design.TwoSampT.scans2 = {...
                 fullfile(testCase.testDataDir, 'su_control06', 'cn_sess1', 'con_0001.img,1'),...
                 fullfile(testCase.testDataDir, 'su_control07', 'cn_sess1', 'con_0001.img,1'),...
                 fullfile(testCase.testDataDir, 'su_control08', 'cn_sess1', 'con_0001.img,1'),...
                 fullfile(testCase.testDataDir, 'su_control09', 'cn_sess1', 'con_0001.img,1'),...
                 fullfile(testCase.testDataDir, 'su_control10', 'cn_sess1', 'con_0001.img,1')};
            %testCase.matlabbatch{1}.cfg_snpm.Design.TwoSampT.group_memb = 'A A A A A B B B B B';
        end
    end
    
    methods (Test)
        % No covariate, no variance smoothing
        function test_twosample_1(testCase)
            testCase.testName = 'twosample_1';
        end

        % With 1 covariate
        function test_twosample_cov(testCase)
            testCase.compaWithSpm = false; % Not matching SPM orthog?
            
            testCase.testName = 'twosample_cov';
            
            testCase.matlabbatch{1}.cfg_snpm.Design.TwoSampT.cov.c = [1 5 2 21 0 3 6 14 8 5];
            testCase.matlabbatch{1}.cfg_snpm.Design.TwoSampT.cov.cname = 'Age';
        end

        % With 3 covariates
        function test_twosample_cov3(testCase)
            testCase.compaWithSpm = false; % Not matching SPM orthog?
            
            testCase.testName = 'twosample_cov3';
            
            testCase.matlabbatch{1}.cfg_snpm.Design.TwoSampT.cov(1).c = [1 5 2 21 0 3 6 14 8 5];
            testCase.matlabbatch{1}.cfg_snpm.Design.TwoSampT.cov(1).cname = 'Age';
            testCase.matlabbatch{1}.cfg_snpm.Design.TwoSampT.cov(2).c = [1 3 5 7 3 5 11 7 8 4];
%             u = [1 5 2 21 0 3 6 14 8 5];
%             x=u(:).'/norm(u);
%             yz=null(x).';
%             xyz=[x;yz]*norm(u);
%             testCase.matlabbatch{1}.cfg_snpm.Design.TwoSampT.cov(2).c = xyz(2,:);
            testCase.matlabbatch{1}.cfg_snpm.Design.TwoSampT.cov(2).cname = 'Height';
            testCase.matlabbatch{1}.cfg_snpm.Design.TwoSampT.cov(3).c = [-1 0.5 0.6 -0.1 2 1 1.5 0.5 1 -1];
            testCase.matlabbatch{1}.cfg_snpm.Design.TwoSampT.cov(3).cname = 'Width';
        end

        % With variance smoothing
        function test_twosample_var(testCase)
            testCase.testName = 'twosample_var';
            
            testCase.matlabbatch{1}.cfg_snpm.Design.TwoSampT.vFWHM = [6 6 6];
        end

        % With approximate test
        function test_twosample_approx(testCase)
            testCase.testName = 'twosample_approx';
            
            rand('seed',200);
            
            testCase.matlabbatch{1}.cfg_snpm.Design.TwoSampT.scans1(end+1:end+2) = {...
                 fullfile(testCase.testDataDir, 'su_control11', 'cn_sess1', 'con_0001.img,1'),...
                 fullfile(testCase.testDataDir, 'su_control12', 'cn_sess1', 'con_0001.img,1')...
                 };
            testCase.matlabbatch{1}.cfg_snpm.Design.TwoSampT.scans2(end+1) = {...
                 fullfile(testCase.testDataDir, 'su_control13', 'cn_sess1', 'con_0001.img,1')...
                 };
             
            testCase.matlabbatch{1}.cfg_snpm.Design.TwoSampT.nPerm = 100;
        end
    end
    
    methods % Start with last method...
        function complete_batch(testCase)
            % Find the result directory for the batch execution and the
            % corresponding result directory computed manually using the
            % original spm2-like interface
            testCase.batchResDir = fullfile(testCase.parentDataDir, 'results', 'batch', testCase.testName);
            testCase.interResDir = fullfile(spm_str_manip(testCase.batchResDir,'hh'), 'interactive', testCase.testName);
            testCase.matlabbatch{1}.cfg_snpm.Design.TwoSampT.dir = {testCase.batchResDir};
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
            
            % % Equal variance
            %factoDesign.des.t2.variance = 0;
            
            testCase.spmBatch{1}.spm.stats.factorial_design = factoDesign;
        end
    end
end

