function snpmui = snpm_bch_ui_ANOVAwithinS
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

DesNm = 'Within Subject ANOVA, k diffs/contrasts per subject';
DesFile = mfilename;
DesHelp = {'',...
	  'stuff.',...
	  '',...
	  'stuff ',...
	  '    stuff',...
	  'stuff.',...
	  };

%% Questions


% Scans per Subject
scans_sub            = cfg_entry;
scans_sub.name       = 'Scans per subject';
scans_sub.tag        = 'Anovawithin_scans';
scans_sub.strtype    = 'i';
scans_sub.val        = {};
scans_sub.num        = [1 1];
scans_sub.help       = {'This is the number of subjects'};

% Covariate Value
cv_none         = cfg_const;
cv_none.tag     = 'cv_none';
cv_none.name    = 'None';
cv_none.val     = {1};
cv_none.help    = {'Covariate value = none'};

cov_Val         = cfg_entry;
cov_Val.tag     = 'cov_Val';
cov_Val.name    = 'Covariate';%arbitary name
cov_Val.help    = {'Help'};
cov_Val.strtype = 'e';
cov_Val.num     = [1 Inf];

cv_one         = cfg_branch;
cv_one.tag     = 'cv_one';
cv_one.name    = 'Enter Different Covariate Value';
cv_one.val     = {cov_Val};
cv_one.help    = {'Help'};

covariate         = cfg_choice;
covariate.tag     = 'covariate';
covariate.name    = 'Covariate Value'; %arbitary name
covariate.val     = {cv_none };
covariate.help    = {'Help'};
covariate.values  = {cv_none cov_Val };


%% Executable Branch
snpmui = snpm_bch_ui(DesNm,DesFile,DesHelp,{covariate});
