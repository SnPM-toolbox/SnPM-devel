function snpmBatch = tbx_cfg_snpm()
%  Function called by SPM to create the SnPM submenu in SPM -> Tools in the
%  batch Window. This function must be called "tbx_cfg_snpm" in order to be
%  automatically detected and called by SPM.
% 
% OUTPUT
%  snpmBatch      - Design menu structure
%
%_______________________________________________________________________
% Camille Maumet

addpath(fullfile(spm('dir'),'toolbox','snpm'));
addpath(fullfile(spm('dir'),'toolbox','snpm', 'config'));

snpmBatch = snpm_cfg_master;
end