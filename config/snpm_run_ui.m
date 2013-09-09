function out = snpm_run_ui(job)
% 
% This collects menu items specified by the design, checks their
% validity, and calls snpm_ui in batch mode, which will then save
% a suitable SnPMcfg.mat configuration file, 
%_______________________________________________________________________
% Thomas Nichols, Emma Thomas
% $Id$

% Based on Volkmar Glauche's MatlabBatch example code (Rev 1716)

snpm_ui(job)
out.SnPMcfg = {fullfile(job.dir{1},'SnPMcfg.mat')};

% function vout = cfg_example_vout_add1(job)
% % Determine what outputs will be present if this job is run. In this case,
% % the structure of the inputs is fixed, and the output is always a single
% % number. Note that input items may not be numbers, they can also be
% % dependencies.

% vout = cfg_dep;                        % The dependency object
% vout.sname      = 'SnPM: SnPMcfg config file';       % Displayed dependency name
% vout.src_output = substruct('()',{1}); % The output subscript reference. This could be any reference into the output variable created during computation

