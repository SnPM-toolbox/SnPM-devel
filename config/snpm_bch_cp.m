function snpmcp = snpm_bch_cp
% Setup menus to run SnPM analysis
%_______________________________________________________________________
% Copyright (C) 2013 The University of Warwick
% Id: snpm_bch_cp.m  SnPM13 2013/10/12
% Thomas Nichols, Emma Thomas
%
% Based on Volkmar Glauche's MatlabBatch example code (Rev 1716)

%% Input Items for our SnPM Compute Run
% Input config file
spmmat         = cfg_files;
spmmat.tag     = 'spmmat';
% spmmat.name    = 'spmmat.mat configuration file';
spmmat.name    = 'SPM.mat configuration file';
spmmat.help    = {'Computes the permutation analysis specified in the selected SPM.mat file.'};
spmmat.filter = 'mat';
spmmat.ufilter = 'SPM.mat';
spmmat.num     = [1 1];

%% Executable Branch
snpmcp      = cfg_exbranch;       % This is the branch that has information about how to run this module
snpmcp.name = 'Compute';             % The display name
snpmcp.tag  = 'cp'; % The name appearing in the harvested job structure. This name must be unique among all items in the val field of the superior node
snpmcp.val  = {spmmat};    % The items that belong to this branch. All items must be filled before this branch can run or produce virtual outputs
%snpmcp.prog = @snpm_run_cp;  % A function handle that will be called with the harvested job to run the computation
snpmcp.prog = @cg_snpm_estimate;
snpmcp.vout = @snpm_bch_cp_vout;
snpmcp.help = {'Runs a configured SnPM.'};

function vout = snpm_bch_cp_vout(job)
% Determine what outputs will be present if this job is run. In this case,
% the structure of the inputs is fixed, and the output is always a single
% number. Note that input items may not be numbers, they can also be
% dependencies.

%%%% Later, add the following:
% lP{+,-} images
% lP_FDR{+,-} images
% lP_FWE{+,-} images
% snpmT{+,-} or snpmF images  - determine from job (set by explict flag?) whether T or F
% 
vout = cfg_dep;                        % The dependency object
vout.sname      = 'SnPM.mat results file';       % Displayed dependency name
vout.src_output = substruct('.','SnPM'); % The output subscript reference. This could be any reference into the output variable created during computation

