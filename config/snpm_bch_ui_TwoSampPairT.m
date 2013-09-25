function snpmui = snpm_bch_ui_TwoSampPairT
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

snpm_defaults 

rev = '$Rev: 1716 $'; %#ok

DesNm = '2 Groups: Test diff of response; 2 conditions, 1 scan per condition';
DesFile = mfilename;
DesHelp = {'Create design and permutation matrix appropriate for two group, two condition interaction analyses where there is just *one* scan per condition.  Only the interaction between activation and group is examined. Use snpm_MSA2x to examine intragroup (main) effects.',...
            '',...
            'Number of permutations. There are nTot-choose-nGrp possible permutations, where nTot is the total number of subjects (and scans) and nGrp is the size of one of the groups. prod(1:nTot)/prod(1:nGrp)^2',...
	  };

%% Questions
% Reuse flexible factorial component of SPM to get the design files and
% associated conditions
factorial_design = spm_cfg_factorial_design();
sub_files = factorial_design.val{2}.values{end}.val{2}.val{1};

% We expect exactly two files per subject
sub_files.values{1}.val{1}.num = [2 2];
sub_files.values{1}.val{1}.help = {'Select a pair of scans (2 conditions)'};
% We need to define the order of the two conditions
sub_files.values{1}.val{2}.num = [1 2];
sub_files.values{1}.val{2}.name = 'Scan index';
sub_files.values{1}.val{2}.help = {'Enter scan index: (1 2|2 1)'};
sub_files.values{1}.val{2}.tag = 'scindex';

% ---------------------------------------------------------------------
% scans1 Group 1 scans
% ---------------------------------------------------------------------
scans1         = cfg_branch;
scans1.tag     = 'scans1';
scans1.name    = 'Group 1 scans';
scans1.help    = {'Select the images from sample 1.  They must all have the same image dimensions, orientation, voxel size etc.'};
scans1.val     = {sub_files}; 

% ---------------------------------------------------------------------
% scans2 Group 2 scans
% ---------------------------------------------------------------------
scans2         = cfg_branch;
scans2.tag     = 'scans2';
scans2.name    = 'Group 2 scans';
scans2.help    = {'Select the images from sample 2.  They must all have the same image dimensions, orientation, voxel size etc.'};
scans2.val     = {sub_files};


% ---------------------------------------------------------------------
% c Vector
% ---------------------------------------------------------------------
c         = cfg_entry;
c.tag     = 'c';
c.name    = 'Vector';
c.help    = {'Vector of covariate values'};
c.strtype = 'e';
c.num     = [Inf 1];
% ---------------------------------------------------------------------
% cname Name
% ---------------------------------------------------------------------
cname         = cfg_entry;
cname.tag     = 'cname';
cname.name    = 'Name';
cname.help    = {'Name of covariate'};
cname.strtype = 's';
cname.num     = [1 Inf];

% ---------------------------------------------------------------------
% mcov Covariate
% ---------------------------------------------------------------------
mcov         = cfg_branch;
mcov.tag     = 'mcov';
mcov.name    = 'Covariate';
mcov.val     = {c cname };
mcov.help    = {'Add a new covariate to your experimental design'};
% ---------------------------------------------------------------------
% generic Covariates
% ---------------------------------------------------------------------
generic         = cfg_repeat;
generic.tag     = 'generic';
generic.name    = 'Covariates';
generic.help    = {'Covariates'};
generic.values  = {mcov };
generic.num     = [0 Inf];


%% Executable Branch
snpmui = snpm_bch_ui(DesNm,DesFile,DesHelp,{scans1, scans2, generic}, true);
