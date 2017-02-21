% Perform non-regression tests on RapidPT in SnPM. 
% Check that results obtained using the batch version are identical to the 
% results computed manually (using the interactive GUI).
%_______________________________________________________________________
% Copyright (C) 2013 The University of Warwick
% Id: test_rapidpt.m  SnPM13 2017/21/12
% Camille Maumet
classdef test_rapidpt < generic_test_snpm
    properties
    end
    
    methods (TestMethodSetup)
        function setGlobals(testCase)           
            % Random number generator should not be initialised with a
            % shuffled seed
            global SnPMdefs
            SnPMdefs.RapidPT = 2;
            rng(1);
            % Run the tests in command line mode (no display)
            global defaults;
            defaults.cmdline = true;
            
        end
        
        function create_basis_matlabbatch(testCase)
            
            for k = 1:7
                testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampT.scans1{k,1} = ...
                     fullfile(testCase.testDataDir, ['test_data_', num2str(k, '%02.0f'), '.nii']);
            end
            for k = 18:24
                testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampT.scans2{k-17,1} = ...
                     fullfile(testCase.testDataDir, ['test_data_gr2_', num2str(k, '%02.0f'), '.nii']);
            end
        end
        
    end
    
    methods (Test)
        % Test maxnull kldiv
        function test_rapidpt_maxt(testCase)
            testCase.testName = 'test_rapidpt_maxt';
            disp('Testting rapidpt...');
            
            testCase.matlabbatch{1}.spm.tools.snpm.des.TwoSampT.nPerm = 3000;

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
        
%         function create_spm_batch(testCase)
%             
%             factoDesign = testCase.spmBatch{1}.spm.stats.factorial_design;
%             
%             factoDesign.des.t2.scans1 = factoDesign.scans1;
%             
%             % Very important as otherwise scans are mixed in SPM
%             if size(factoDesign.des.t2.scans1, 2) > 1
%                 factoDesign.des.t2.scans1 = factoDesign.des.t2.scans1';
%             end
%             factoDesign.des.t2.scans2 = factoDesign.scans2;
%             if size(factoDesign.des.t2.scans2, 2) > 1
%                 factoDesign.des.t2.scans2 = factoDesign.des.t2.scans2';
%             end
%                 
%             factoDesign = rmfield(factoDesign, 'scans1');
%             factoDesign = rmfield(factoDesign, 'scans2');
%             if isfield(factoDesign, 'ST')
%                 factoDesign = rmfield(factoDesign, 'ST');
%             end
%             
%             % % Equal variance
%             %factoDesign.des.t2.variance = 0;
%             
%             testCase.spmBatch{1}.spm.stats.factorial_design = factoDesign;
%         end
    end
end

