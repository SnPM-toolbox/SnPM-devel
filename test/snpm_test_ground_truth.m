% Compute ground trouth for testing purposes
%
% Update snpm_test_config before running this function.
%_______________________________________________________________________
% Copyright (C) 2013-14 The University of Warwick
% Id: test_oneSample.m  SnPM13.01 2014/01/31
% Camille Maumet
function snpm_test_ground_truth()
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
      error('Re-computation of groung truth stopped!');
    case 'yes',
      disp('Start re-computation of ground truth')
      
      resultDir = fullfile(spm_str_manip(testDataDir, 'h'), 'results');

      % ---- one-sample tests ----

      % test 1: onesample_1
      matlabbatch{1}.spm.tools.snpm.des.OneSampT.P = {
          fullfile(testDataDir, 'su_control01', 'cn_sess1', 'con_0001.img,1')
          fullfile(testDataDir, 'su_control02', 'cn_sess1', 'con_0001.img,1')
          fullfile(testDataDir, 'su_control03', 'cn_sess1', 'con_0001.img,1')
          fullfile(testDataDir, 'su_control04', 'cn_sess1', 'con_0001.img,1')
          fullfile(testDataDir, 'su_control05', 'cn_sess1', 'con_0001.img,1')
      };
      matlabbatch{1}.spm.tools.snpm.des.OneSampT.dir = {fullfile(resultDir, 'GT', 'onesample_1')}; 
      SnPMmatFile = fullfile(resultDir, 'GT', 'onesample_1', 'SPM.mat');
      matlabbatch{2}.spm.tools.snpm.cp.snpmcfg = {SnPMmatFile};

      % Results   
      % Uncorrected voxel-wise p<0.1
      matlabbatch{3}.spm.tools.snpm.inference.SnPMmat = {SnPMmatFile};
      matlabbatch{3}.spm.tools.snpm.inference.Thr.Vox.VoxSig.Pth = 0.10;
      matlabbatch{3}.spm.tools.snpm.inference.Tsign = 1;
      matlabbatch{3}.spm.tools.snpm.inference.WriteFiltImg.name = 'SnPM_filtered_10none.nii';

      % Rename uncorrected p<0.1
      matlabbatch{end}.spm.tools.snpm.inference.WriteFiltImg.name = 'SnPMt_filtered_vox_unc_p10.nii';

      % Uncorrected voxel-wise TorF > 1.6
      matlabbatch{end+1}.spm.tools.snpm.inference.SnPMmat = {SnPMmatFile};
      matlabbatch{end}.spm.tools.snpm.inference.Thr.Vox.VoxSig.TFth = 1.6;
      matlabbatch{end}.spm.tools.snpm.inference.WriteFiltImg.name = 'SnPMt_filtered_vox_unc_t16.nii';

      % FWE voxel-wise p<0.5
      matlabbatch{end+1}.spm.tools.snpm.inference.SnPMmat = {SnPMmatFile};
      matlabbatch{end}.spm.tools.snpm.inference.Thr.Vox.VoxSig.FWEth = 0.5;
      matlabbatch{end}.spm.tools.snpm.inference.WriteFiltImg.name = 'SnPMt_filtered_vox_fwe_p50.nii'; 

      % FDR voxel-wise p<0.5
      matlabbatch{end+1}.spm.tools.snpm.inference.SnPMmat = {SnPMmatFile};
      matlabbatch{end}.spm.tools.snpm.inference.Thr.Vox.VoxSig.FDRth = 0.5;
      matlabbatch{end}.spm.tools.snpm.inference.WriteFiltImg.name = 'SnPMt_filtered_vox_fdr_p50.nii'; 

      spm_jobman('run', matlabbatch);

  end % switch
end

