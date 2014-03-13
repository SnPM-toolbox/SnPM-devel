% Compute ground trouth for testing purposes
%
% Update snpm_test_config before running this function.

%_______________________________________________________________________
% Copyright (C) 2013-14 The University of Warwick
% Id: snpm_test_ground_truth.m  SnPM8 2014/01/31
% Camille Maumet
function snpm_test_ground_truth(redo)

if nargin == 0
    redo = false;
end

snpm_test_config;
global testDataDir;

disp('Recomputing ground thruth for SnPM testing...')

if isempty(testDataDir)
    error('Test data directory not set, please update snpm_test_config');
end

disp(['Current version of SnPM: ' snpm('ver')])
disp(['Path: ' which('snpm')])

buttonName = questdlg(['Current version of SnPM: ' snpm('ver') '. Is that ok?'],...
    'Check SnPM version','yes','no','no');

switch buttonName,
    case 'no',
        error('Re-computation of ground truth stopped!');
    case 'yes',
        disp(sprintf('\nStart re-computation of ground truth\n'))
        
        resultDir = fullfile(spm_str_manip(testDataDir, 'h'), 'results');
        gtDirName = ['GT_' strrep(snpm('ver'), '.', '')];
        
        cwd = pwd;
        
        testOneSample = {} % { 'onesample_slice', 'onesample_ancova', ...
%             'onesample_grandmean_50', 'onesample_grandmean_145', ...
%             'onesample_propscaling_to_user', 'onesample_propscaling', ...
%             'onesample_cluster', 'onesample_cluster_predefined', ...
%             'onesample_1', 'onesample_propscaling', 'onesample_approx', ...
%             'onesample_var', 'onesample_cov3', 'onesample_cov'
%             } ;
         
        testTwoSample = {'twosample_var', ...
                        'twosample_approx', 'twosample_propscaling', ...
                        'twosample_propscaling_to_user', ...
                        'twosample_grandmean_145', 'twosample_grandmean_50',...
                        'twosample_ancova'}%'twosample_cov', 'twosample_cov3' 'twosample_cluster_predef_stat'}; %'twosample_cluster_predefined', 'twosample_cluster', 'twosample_1'
        testOneSubTwoSample = {}% {'onesub_twocondrepl_1_other_design', ...
%              'onesub_twocondrepl_1', 'onesub_twocondrepl_var'};
        allTests = [testOneSample testTwoSample testOneSubTwoSample];
        
        for i = 1:numel(allTests)
            currTest = allTests{i};
            
            % ---- one-sample tests ----
            resDir = fullfile(resultDir, gtDirName, allTests{i});
            if ~isdir(resDir)
                mkdir(resDir)
                cd(resDir)
            end
            
            cfgFile = spm_select('FPList', resDir, '^SnPMcfg\.mat$');           
            
            switch(allTests{i})
                % *** One-sample test ***
                case {'onesample_1'}
                    if isempty(cfgFile) || redo
                        design_one_sample_test(testDataDir, resDir, '0', {}, '0')
                    end
                    
                case {'onesample_cov'}
                    if isempty(cfgFile) || redo
                        design_one_sample_test(testDataDir, resDir, '1', {'1 5 2 21 0'}, '0')
                    end
                    
                case {'onesample_cov3'}
                    if isempty(cfgFile) || redo
                        design_one_sample_test(testDataDir, resDir, '3', ...
                            {'1 1 2 3 1', '0 21 15 18 3', '-1 -0.5 -1 1 0'}, '0')
                    end
                    
                case {'onesample_var'}
                    if isempty(cfgFile) || redo
                        design_one_sample_test(testDataDir, resDir, '0', {}, '6')
                    end
                    
                case {'onesample_approx'}
                    if isempty(cfgFile) || redo
                        rand('seed',200);
                        design_one_sample_test(testDataDir, resDir, '0', {}, '0', 5, '15')
                    end
                
                case {'onesample_propscaling'}
                    if isempty(cfgFile) || redo
                        design_one_sample_test(testDataDir, resDir, ...
                            '0', {}, '0', 5, '', 'proportional scaling')
                    end
                    
                case {'onesample_propscaling_to_user'}
                    if isempty(cfgFile) || redo
                        design_one_sample_test(testDataDir, resDir, ...
                            '0', {}, '0', 5, '', 'proportional scaling', '145')
                    end
                    
                case {'onesample_grandmean_145'}
                    if isempty(cfgFile) || redo
                        design_one_sample_test(testDataDir, resDir, ...
                            '0', {}, '0', 5, '', '', '', 'scaling to overall grand mean', '145')
                    end
                    
                case {'onesample_grandmean_50'}
                    if isempty(cfgFile) || redo
                        design_one_sample_test(testDataDir, resDir, ...
                            '0', {}, '0', 5, '', '', '', 'scaling to overall grand mean', '50')
                    end
                    
                case {'onesample_ancova'}
                    if isempty(cfgFile) || redo
                        design_one_sample_test(testDataDir, resDir, ...
                            '0', {}, '0', 5, '', 'AnCova', '')
                    end
                    
                case {'onesample_slice'}
                    rand('seed',200);
                    if isempty(cfgFile) || redo
                        design_one_sample_test(testDataDir, resDir, ...
                            '0', {}, '0', 17, '15')
                    end
                    
                case {'onesample_cluster', 'onesample_cluster_predefined'}
                    nominalCfg = spm_select('FPList', fullfile(spm_str_manip(resDir, 'h'), 'onesample_1'), '^SnPMcfg\.mat$');
                    if isempty(nominalCfg)
                        error('No nominal config file for SnPM');
                    end
                    copyfile(nominalCfg, resDir);
                    cfgFile = spm_select('FPList', resDir, '^SnPMcfg\.mat$');
                    configSnPM = load(cfgFile);
                    
                    switch(allTests{i})
                        case {'onesample_cluster'}
                            configSnPM.bST = 1;
                            configSnPM.pU_ST_Ut = -1;
                        case {'onesample_cluster_predefined'}
                            configSnPM.bST = 1;
                            configSnPM.pU_ST_Ut = 0.1;     
                    end                    
                   save(cfgFile, '-struct', 'configSnPM')
                   
                % *** Two-sample test ***
                case {'twosample_1'}
                    if isempty(cfgFile) || redo
                        design_two_sample_test(testDataDir, resDir, '0', {}, '0')
                    end   
                    
                case {'twosample_cluster', 'twosample_cluster_predefined', ...
                        'twosample_cluster_predef_stat'}
                    nominalCfg = spm_select('FPList', fullfile(spm_str_manip(resDir, 'h'), 'twosample_1'), '^SnPMcfg\.mat$');
                    if isempty(nominalCfg)
                        error('No nominal config file for SnPM');
                    end
                    copyfile(nominalCfg, resDir);
                    cfgFile = spm_select('FPList', resDir, '^SnPMcfg\.mat$');
                    configSnPM = load(cfgFile);
                    
                    switch(allTests{i})
                        case {'twosample_cluster'}
                            configSnPM.bST = 1;
                            configSnPM.pU_ST_Ut = -1;
                        case {'twosample_cluster_predefined'}
                            configSnPM.bST = 1;
                            configSnPM.pU_ST_Ut = 0.1;  
                        case {'twosample_cluster_predef_stat'}
                            configSnPM.bST = 1;
                            configSnPM.pU_ST_Ut = 2.03;  
                    end                    
                   save(cfgFile, '-struct', 'configSnPM')
                   
                case {'twosample_cov'}
                    if isempty(cfgFile) || redo
                        design_two_sample_test(testDataDir, resDir, '1', {'1 5 2 21 0 3'}, '0')
                    end
                    
                case {'twosample_cov3'}
                    if isempty(cfgFile) || redo
                        design_two_sample_test(testDataDir, resDir, '3', ...
                            {'1 5 2 21 0 3', '1 3 5 7 3 5', '-1 0.5 0.6 -0.1 2 1'}, '0')
                    end
                    
                case {'twosample_var'}
                    if isempty(cfgFile) || redo
                        design_two_sample_test(testDataDir, resDir, '0', {}, '6')
                    end
                    
                case {'twosample_approx'}
                    if isempty(cfgFile) || redo
                        rand('seed',200);
                        design_two_sample_test(testDataDir, resDir, '0', {}, '0', 5, '15')
                    end
                
                case {'twosample_propscaling'}
                    if isempty(cfgFile) || redo
                        design_two_sample_test(testDataDir, resDir, ...
                            '0', {}, '0', 5, '', 'proportional scaling')
                    end
                    
                case {'twosample_propscaling_to_user'}
                    if isempty(cfgFile) || redo
                        design_two_sample_test(testDataDir, resDir, ...
                            '0', {}, '0', 5, '', 'proportional scaling', '145')
                    end
                    
                case {'twosample_grandmean_145'}
                    if isempty(cfgFile) || redo
                        design_two_sample_test(testDataDir, resDir, ...
                            '0', {}, '0', 5, '', '', '', 'scaling to overall grand mean', '145')
                    end
                    
                case {'twosample_grandmean_50'}
                    if isempty(cfgFile) || redo
                        design_two_sample_test(testDataDir, resDir, ...
                            '0', {}, '0', 5, '', '', '', 'scaling to overall grand mean', '50')
                    end
                    
                case {'twosample_ancova'}
                    if isempty(cfgFile) || redo
                        design_two_sample_test(testDataDir, resDir, ...
                            '0', {}, '0', 5, '', 'AnCova', '')
                    end
                
                % *** One-subject two-sample test ***
                    case {'onesub_twocondrepl_1'}
                    if isempty(cfgFile) || redo
                        design_one_sub_two_sample_test(testDataDir, resDir, '0')
                    end
                    
                    case {'onesub_twocondrepl_1_other_design'}
                    if isempty(cfgFile) || redo
                        design_one_sub_two_sample_test(testDataDir, resDir, '0', 10)
                    end
                    
                    case {'onesub_twocondrepl_var'}
                    if isempty(cfgFile) || redo
                        design_one_sub_two_sample_test(testDataDir, resDir, '12')
                    end
                
                otherwise
                    error('undefined test')
            end
            
            imageFiles = spm_select('FPList', resDir, '.*\.(img|hdr|nii)$');
            SnPMcomputeFiles = spm_select('FPList', resDir, '^(SnPM|SnPMt|SnPMucp|XYZ|SnPM_combo|SnPM_ST)\.mat$');
            psFile = spm_select('FPList', resDir, '^spm_.*\.ps$');
            toDel = cellstr(strvcat(imageFiles, SnPMcomputeFiles, psFile));
            if ~isempty(toDel{1})
                for j = 1:numel(toDel)
                    delete(toDel{j});
                end
            end
            
            snpm_cp(resDir);
            
            % Results have to be computed interactively...
            switch(allTests{i})
                case {'onesample_1'}
                    additional_interactive_results(resDir)
                    
                case {'onesample_cluster', 'twosample_cluster'}
                    additional_interactive_results(resDir, 'Voxelwise')
                    additional_interactive_cluster_results(resDir)
                    interactive_cluster_mass_results(resDir)
                    
                case {'onesample_cluster_predefined', 'twosample_cluster_predefined',...
                        'twosample_cluster_predef_stat'}
                    additional_interactive_predefined_cluster_results(resDir)
                    
                case {'onesample_cov', 'onesample_cov3', 'onesample_var', ...
                        'onesample_approx', 'onesample_propscaling', ...
                        'onesample_propscaling_to_user', ...
                        'onesample_grandmean_145', 'onesample_grandmean_50',...
                        'onesample_ancova', 'onesample_slice',...
                        'onesub_twocondrepl_1', 'onesub_twocondrepl_var', ...
                        'onesub_twocondrepl_1_other_design', ...
                        'twosample_1', 'twosample_cov', 'twosample_cov3', 'twosample_var', ...
                        'twosample_approx', 'twosample_propscaling', ...
                        'twosample_propscaling_to_user', ...
                        'twosample_grandmean_145', 'twosample_grandmean_50',...
                        'twosample_ancova'}
                    interactive_results(resDir, 'SnPM_filtered_10none', 'P', 'None', '0.1');
                    
                otherwise
                    error('undefined test')
            end
        end
        
end % switch
end

function additional_interactive_results(resDir, clusterInference)
if nargin == 1
    clusterInference = '';
end
interactive_results(resDir, 'SnPMt_filtered_vox_unc_p10', 'P', 'None', '0.1');
interactive_results(resDir, 'SnPMt_filtered_vox_unc_t16', 'T', 'None', '1.6', clusterInference);
interactive_results(resDir, 'SnPMt_filtered_vox_fwe_p10', 'T', 'FWE', '0.1', clusterInference);
interactive_results(resDir, 'SnPMt_filtered_vox_fdr_p70', 'P', 'FDR', '0.7');
end

function additional_interactive_predefined_cluster_results(resDir)
interactive_results(resDir, 'SnPMt_filtered_vox_unc_p10', 'P', 'None', '0.1');
interactive_results(resDir, 'SnPMt_filtered_clus_p10_unc_p10', 'T', 'Uncorr', '0.1', 'Clusterwise', '', 'P-value');
interactive_results(resDir, 'SnPMt_filtered_clus_p10_fwe_p50', 'T', 'FWE', '0.5', 'Clusterwise', '', 'P-value');
end

function additional_interactive_cluster_results(resDir)
interactive_results(resDir, 'SnPMt_filtered_clus_4_unc_p10', 'T', 'Uncorr', '0.1', 'Clusterwise', '4', 'P-value');
interactive_results(resDir, 'SnPMt_filtered_clus_4_unc_k6', 'T', 'Uncorr', '6', 'Clusterwise', '4', 'ClusterSize');
interactive_results(resDir, 'SnPMt_filtered_clus_4_fwe_p50', 'T', 'FWE', '0.5', 'Clusterwise', '4', 'P-value');
interactive_results(resDir, 'SnPMt_filtered_clus_5_fwe_p50', 'T', 'FWE', '0.5', 'Clusterwise', '5', 'P-value');
end

% Instructions for interactive one-subject two-sample test
function design_one_sub_two_sample_test(testDataDir, resDir, varSmoothing, nSubjects)
    if nargin < 4
        nSubjects = 12;
    end

    cwd = pwd;
    cd(testDataDir)
    % There is no snpmcfg.mat start snpm_ui and create it
    % interactively (with instructions for user)
    disp('* Select design type: SingleSub: Two Sample T test; 2 conditions');
    disp(['* # replications per condition: ' num2str(nSubjects/2)]);
    if nSubjects == 12
        disp('* Size of exchangeability block: 4');
    else
        disp('* Size of exchangeability block: 10');
    end
    disp('* Select scans in time order:')
    for i = 1:nSubjects
        disp(sprintf(['\t' fullfile(testDataDir, ['test_data_' num2str(i, '%02.0f') '.nii'])]));
    end
    disp(['* Enter conditions index (B/A) [' num2str(nSubjects) ']: ' repmat('AB', 1, nSubjects/2)]);
    common_choices(1, varSmoothing, '', '', '', '')
    snpm_ui_and_copy_config(cwd, resDir);
end

% Instructions for interactive two-sample tests
function design_two_sample_test(testDataDir, resDir, numCovariates, ...
                valueCov, varSmoothing, nSubjects, nPerm, propScaling, ...
                userPropScaling, grandMeanScaling, userGrandMean)
    if nargin < 11
        userGrandMean = '';
        if nargin < 10
            grandMeanScaling = '';
            if nargin < 9
                userPropScaling = '50';
                if nargin < 8
                    propScaling = '';
                    if nargin < 7
                        nPerm = '';
                        if nargin < 6
                            nSubjects = 5;
                            if nargin < 5
                                varSmoothing = '0';
                            end
                        end
                    end
                end
            end
        end
    end
    
    cwd = pwd;
    cd(testDataDir)
    % There is no snpmcfg.mat start snpm_ui and create it
    % interactively (with instructions for user)
    disp('* Select design type: 2 Groups: Two Sample T test; 1 scan/subject');
    disp('* Select all scans:');
    for i = 1:3
        disp(sprintf(['\t' ...
                fullfile(testDataDir, ['test_data_' num2str(i, '%02.0f') '.nii'])]));
    end
    for i = 18:20
        disp(sprintf(['\t' ...
                fullfile(testDataDir, ['test_data_gr2_' num2str(i, '%02.0f') '.nii'])]));
    end
    disp(['Enter subject index (A/B) [6]: AAABBB'])
    disp(['* # of confounding covariates: ' numCovariates])
    for i = 1:str2double(numCovariates)
        disp(['* [5] - Covariate ' num2str(i) ': ' valueCov{i}])
    end
    if isempty(nPerm)
        approx = 'No';
    else
        approx = 'Yes';
    end
    disp(['* ' num2str(2^nSubjects) ' Perms. Use approx. test?: ' approx])
    if ~isempty(nPerm)
        disp(['* # perms. to use? (Max ' num2str(2^nSubjects) '): ' nPerm])
    end
    common_choices(5, varSmoothing, propScaling, grandMeanScaling, ...
        userGrandMean, userPropScaling);
    snpm_ui_and_copy_config(cwd, resDir)
end

% Instructions for interactive one-sample tests
function design_one_sample_test(testDataDir, resDir, numCovariates, ...
                valueCov, varSmoothing, nSubjects, nPerm, propScaling, ...
                userPropScaling, grandMeanScaling, userGrandMean)
    if nargin < 11
        userGrandMean = '';
        if nargin < 10
            grandMeanScaling = '';
            if nargin < 9
                userPropScaling = '50';
                if nargin < 8
                    propScaling = '';
                    if nargin < 7
                        nPerm = '';
                        if nargin < 6
                            nSubjects = 5;
                        end
                    end
                end
            end
        end
    end
    
    cwd = pwd;
    cd(testDataDir)
    % There is no snpmcfg.mat start snpm_ui and create it
    % interactively (with instructions for user)
    disp('* Select design type: MultiSub: One Sample T test on differences; 1 condition');
    disp('* Select all scans:');
    for i = 1:nSubjects
        disp(sprintf(['\t' fullfile(testDataDir, ['test_data' num2str(i, '%02.0f') '.nii'])]))
    end
    disp(['* # of confounding covariates: ' numCovariates])
    for i = 1:str2double(numCovariates)
        disp(['* [5] - Covariate ' num2str(i) ': ' valueCov{i}])
    end
    if isempty(nPerm)
        approx = 'No';
    else
        approx = 'Yes';
    end
    disp(['* ' num2str(2^nSubjects) ' Perms. Use approx. test?: ' approx])
    if ~isempty(nPerm)
        disp(['* # perms. to use? (Max ' num2str(2^nSubjects) '): ' nPerm])
    end
    common_choices(nSubjects, varSmoothing, propScaling, grandMeanScaling, ...
        userGrandMean, userPropScaling);
    snpm_ui_and_copy_config(cwd, resDir)
end

function snpm_ui_and_copy_config(cwd, resDir)
    snpm_ui
    cfgFile = spm_select('FPList', pwd, '^SnPMcfg\.mat$');
    if isempty(cfgFile)
        error(['No SnPM configuration file in' pwd]);
    end
    movefile(cfgFile, resDir);
    cd(cwd);
end

function common_choices(nSubjects, varSmoothing, propScaling, ...
    grandMeanScaling, userGrandMean, userPropScaling)
    disp(['* FWHM(mm) for Variance smooth: ' varSmoothing])
    if nSubjects > 5
        disp(['* 17 scans: Work volumetrically?: no'])
    end
    disp('* Collect Supra-Threshold stats?: No')
    if isempty(propScaling)
        disp('* Select global normalisation: <no Global normalisation>')
    else
        disp(['* Select global normalisation: ' propScaling])
        if strcmp(propScaling, 'proportional scaling')
            disp(['* Propsca global mean to: ' userPropScaling])
        end
    end
    if isempty(propScaling)
        if isempty(grandMeanScaling)
            disp('* grand mean scaling: <no grand Mean scaling>')
        else
            disp(['* grand mean scaling: ' grandMeanScaling])
            disp(['* scale overall grand mean to...: ' userGrandMean])
        end
    end
    if ~isempty(propScaling) || ~isempty(grandMeanScaling)
        disp('* Select global calculation...: mean voxel value (within per image fullmean/8 mask)')
    end
    disp('* Threshold masking: none')
    disp('* Analysis mask?: No')
end

function interactive_cluster_mass_results(resDir)
    disp('* Positive or negative effects?: +ve');
    disp('* Write filtered statistic img?: yes');
    disp('* Filename?: SnPMt_filtered_cluss_mass');
    disp('* Corrected p value for filtering?: 0.05' );
    disp('* Ut (0.00<p<0.01 | xx<t<xx): 3.8');
    disp('* Theta value for voxel-cluster combining?: 0.5');
    disp('* Choose combining function: Mass');
    
	snpm_combo_pp(resDir);
end

function interactive_results(resDir, fileName, PorT, correctedThresh, ...
    threshValue, clusterInference, formingThresh, threshType)
if nargin < 6
    clusterInference = '';
end
if nargin < 8
    threshType = '';
end

disp(sprintf('\n--- Please parametrise the inference according to:'));

disp('* Positive or negative effects?: +ve');
disp('* Write filtered statistic img?: yes');
disp(['* Filename?: ' fileName]);
disp(['* Results for which img?: ' PorT]);
if ~isempty(clusterInference)
    disp(['* Inference method?: ' clusterInference]);
    if strcmp(clusterInference, 'Clusterwise') && ~isempty(formingThresh)
        disp(['* Clus-def thresh(pseudo t>xxx)?: ' formingThresh]);
    end
end
disp(['* Use corrected threshold?: ' correctedThresh]);
if strcmp(correctedThresh, 'Uncorr')
    disp(['* Define uncorrected: ' threshType]);
end
if strcmp(threshType, 'ClusterSize')
    disp(['* Uncorr cluster size threshold: ' threshValue]);
else
    disp(['* Uncorrected p value threshold: ' threshValue]);
end
snpm_pp(resDir);
end

