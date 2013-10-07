% Compare results (t-map and positive effects filtered map T>3) of
% interactive and batch one subject, 2 conditions with replication
classdef test_onesub_twocondrepl < matlab.unittest.TestCase & generic_test_snpm
    properties
    end
    
    methods (TestMethodSetup)
        function create_basis_matlabbatch(testCase)
            testCase.compaWithSpm = false;
            
            % Fill the design part in the batch
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampTss.P = {
                 fullfile(testCase.testDataDir, 'PET_motor', 's8np01160em01R.img,1')
                 fullfile(testCase.testDataDir, 'PET_motor', 's8np01160em02R.img,1')
                 fullfile(testCase.testDataDir, 'PET_motor', 's8np01160em03R.img,1')
                 fullfile(testCase.testDataDir, 'PET_motor', 's8np01160em04R.img,1')
                 fullfile(testCase.testDataDir, 'PET_motor', 's8np01160em05R.img,1')
                 fullfile(testCase.testDataDir, 'PET_motor', 's8np01160em06R.img,1')
                 fullfile(testCase.testDataDir, 'PET_motor', 's8np01160em07R.img,1')
                 fullfile(testCase.testDataDir, 'PET_motor', 's8np01160em08R.img,1')
                 fullfile(testCase.testDataDir, 'PET_motor', 's8np01160em09R.img,1')
                 fullfile(testCase.testDataDir, 'PET_motor', 's8np01160em10R.img,1')
                 fullfile(testCase.testDataDir, 'PET_motor', 's8np01160em11R.img,1')
                 fullfile(testCase.testDataDir, 'PET_motor', 's8np01160em12R.img,1')
                 };
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampTss.condidx = [1 2 1 2 1 2 1 2 1 2 1 2];
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampTss.Tss_repc = 6;
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampTss.TwosampTss_Block = 4;
        end
    end
    
    methods (Test)
        % No covariate, no variance smoothing
        function test_onesub_twocondrepl_1(testCase)
            % Nominal test
            testCase.testName = 'onesub_twocondrepl_1';
        end

        % With variance smoothing
        function test_onesub_twocondrepl_var(testCase)
            testCase.testName = 'onesub_twocondrepl_var';
            
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampTss.vFWHM = [12 12 12];
        end
    end
    
    methods (TestMethodTeardown)
        function complete_batch(testCase)
            % Find the result directory for the batch execution and the
            % corresponding result directory computed manually using the
            % original spm2-like interface
            testCase.batchResDir = fullfile(testCase.parentDataDir, 'results', 'batch', testCase.testName);
            testCase.interResDir = fullfile(spm_str_manip(testCase.batchResDir,'hh'), 'interactive', testCase.testName);
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampTss.dir = {testCase.batchResDir};
        end
    end
end