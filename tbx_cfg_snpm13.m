function snpmBatch = tbx_cfg_snpm13()
%  Function called by SPM to create the SnPM submenu in SPM -> Tools in the
%  batch Window. This function must be called "tbx_cfg_snpm" in order to be
%  automatically detected and called by SPM.
% 
% OUTPUT
%  snpmBatch      - Design menu structure
%
%_______________________________________________________________________
% Copyright (C) 2013 The University of Warwick
% Id: tbx_cfg_snpm13.m  SnPM13 2013/10/12
% Camille Maumet

toolboxDir = spm_str_manip(mfilename('fullpath'), 'h');

addpath(toolboxDir);
addpath(fullfile(toolboxDir, 'config'));
addpath(fullfile(toolboxDir, 'test'));
addpath(fullfile(toolboxDir, 'test', 'common'));

snpmBatch = snpm_cfg_master;
end