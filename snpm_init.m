function snpm_init
% FORMAT snpm_init
% Initialise SnPM toolbox and SPM matlabbatch.
%_______________________________________________________________________
% Copyright (C) 2013 The University of Warwick
% Id: snpm_init.m  SnPM13 2013/10/12
% Thomas Nichols
% Based on toy_example.m 1716 2008-05-23 08:18:45Z volkmar

rev = '$Rev: 1716 $'; %#ok

if exist('snpm_defaults')~=2
  error('Cannot initialize - SnPM toolbox not found in path')
end
if isempty(whos('global','SnPMdefs'))
  snpm_defaults
end

% the application config resides in the current directory, therefore add
% this to the path
p = fileparts(mfilename('fullpath'));
addpath(fullfile(p,'config'));

% cfg_util initialisation
cfg_util('initcfg');
