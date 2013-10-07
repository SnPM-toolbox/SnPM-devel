function snpmpp = snpm_bch_pp_FWEreport
% Example script that creates an cfg_exbranch to sum two numbers. The
% inputs are entered as two single numbers, the output is just a single
% number.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_example_add1.m 1716 2008-05-23 08:18:45Z volkmar $

rev = '$Rev: 1716 $'; %#ok

%% Input Items for our SnPM Compute Run
% Input SnPM results mat file
snpmres         = cfg_files;
snpmres.tag     = 'snpmreslt';
snpmres.name    = 'SnPM.mat results file';
snpmres.help    = {'Select a SnPM.mat results file.'};
snpmres.filter = 'any';
snpmres.ufilter = 'SnPM.mat';
snpmres.num     = [1 1];

% Input SnPM results mat file
action         = cfg_files;
action.tag     = 'Action';
action.name    = 'SnPM.mat results file';
action.help    = {'Select a SnPM.mat results file.'};
action.filter = 'any';
action.ufilter = 'SnPM.mat';
action.num     = [1 1];


%% Executable Branch
snpmpp      = cfg_exbranch;       % This is the branch that has information about how to run this module
snpmpp.name = 'Results';             % The display name
snpmpp.tag  = 'pp'; % The name appearing in the harvested job structure. This name must be unique among all items in the val field of the superior node
snpmpp.val  = {snpmres action};    % The items that belong to this branch. All items must be filled before this branch can run or produce virtual outputs
snpmpp.prog = @snpm_run_pp;  % A function handle that will be called with the harvested job to run the computation
%snpmpp.vout = @cfg_example_vout_snpmpp; % A function handle that will be called with the harvested job to determine virtual outputs
snpmpp.help = {'Runs a configured SnPM.'};

function vout = snpm_bch_pp_vout(job)
% Determine what outputs will be present if this job is run. In this case,
% the structure of the inputs is fixed, and the output is always a single
% number. Note that input items may not be numbers, they can also be
% dependencies.

%%%%% LATER, add the output to this... i.e. the filter image
vout = cfg_dep;                        % The dependency object
vout.sname      = 'Add1: a + b';       % Displayed dependency name
vout.src_output = substruct('()',{1}); % The output subscript reference. This could be any reference into the output variable created during computation
