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

toolboxDir = spm_str_manip(mfilename('fullpath'), 'h');

addpath(toolboxDir);
addpath(fullfile(toolboxDir, 'config'));

snpmBatch = snpm_cfg_master;
end