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

disp('Current version of SnPM')
disp(which('snpm'))

buttonName = questdlg(['Current version of SnPM: ' snpm('ver') '. Is that ok?'],...
    'Check SnPM version','yes','no','no');

switch buttonName,
    case 'no',
        error('Re-computation of ground truth stopped!');
    case 'yes',
        disp('Start re-computation of ground truth')
        
        resultDir = fullfile(spm_str_manip(testDataDir, 'h'), 'results');
        gtDirName = ['GT_' strrep(snpm('ver'), '.', '')];
        
        cwd = pwd;
        
        testOneSample = {'onesample_propscaling_to_user'}%'onesample_propscaling','onesample_cluster' 'onesample_cluster_predefined'} %{'onesample_1', 'onesample_propscaling', 'onesample_approx', 'onesample_var', 'onesample_cov3', 'onesample_cov', , } % };
        allTests = testOneSample;
        
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
                    
                case {'onesample_cluster'}
                    additional_interactive_results(resDir, 'Voxelwise')
                    additional_interactive_cluster_results(resDir)
                    interactive_cluster_mass_results(resDir)
                    
                case {'onesample_cluster_predefined'}
                    additional_interactive_predefined_cluster_results(resDir)
                    
                case {'onesample_cov', 'onesample_cov3', 'onesample_var', ...
                        'onesample_approx', 'onesample_propscaling', ...
                        'onesample_propscaling_to_user'}
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
interactive_results(resDir, 'SnPMt_filtered_vox_fwe_p50', 'T', 'FWE', '0.5', clusterInference);
interactive_results(resDir, 'SnPMt_filtered_vox_fdr_p50', 'P', 'FDR', '0.5');
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

function design_one_sample_test(testDataDir, resDir, numCovariates, ...
                valueCov, varSmoothing, nSubjects, nPerm, propScaling, ...
                userPropScaling)
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
    
    cwd = pwd;
    cd(resDir)
    % There is no snpmcfg.mat start snpm_ui and create it
    % interactively (with instructions for user)
    disp('* Select design type: MultiSub: One Sample T test on differences; 1 condition');
    disp('* Select all scans:');
    for i = 1:nSubjects
        disp(fullfile(testDataDir, ['su_control' num2str(i, '%02.0f')], ...
            'cn_sess1', 'con_0001.img,1'))
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
    disp(['* FWHM(mm) for Variance smooth: ' varSmoothing])
    disp('* Collect Supra-Threshold stats?: No')
    if isempty(propScaling)
        disp('* Select global normalisation: <no Global normalisation>')
    elseif strcmp(propScaling, 'proportional scaling')
        disp(['* Select global normalisation: ' propScaling])
        disp(['* Propsca global mean to: ' userPropScaling])
        disp('* Select global calculation...: mean voxel value (within per image fullmean/8 mask)')
    end
    disp('* grand mean scaling: <no grand Mean scaling>')
    disp('* Threshold masking: none')
    disp('* Analysis mask?: No')
    snpm_ui
    cd(cwd);
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

