function out = snpm_run_cp(job)
% Run a SnPM analysis
%
%_______________________________________________________________________
% Thomas Nichols, Emma Thomas
% $Id$

% Based on Volkmar Glauche's MatlabBatch example code (Rev 1716)

rev = '$Rev: 1716 $'; %#ok

snpmdir=fileparts(job.snpmcfg{1});
snpm_cp(snpmdir);
out.SnPM = {fullfile(snpmdir,'SnPM.mat')};
