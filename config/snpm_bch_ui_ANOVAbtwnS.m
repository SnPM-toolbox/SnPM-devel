function snpmui = snpm_bch_ui_ANOVAbtwnS
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

DesNm = 'Between group ANOVA; 1 scan per subject';
DesFile = mfilename;
DesHelp = {'',...
	  'stuff.',...
	  '',...
	  'stuff ',...
	  '    stuff',...
	  'stuff.',...
	  };


%% Questions
% ---------------------------------------------------------------------
% scans Scans [1,2]
% ---------------------------------------------------------------------
scans         = cfg_files;
scans.tag     = 'scans';
scans.name    = 'Scans';
scans.help    = {'Select images for this group. '};
scans.filter = 'image';
scans.ufilter = '.*';
scans.num     = [1 Inf];
% ---------------------------------------------------------------------
% group group
% ---------------------------------------------------------------------
group         = cfg_branch;
group.tag     = 'group';
group.name    = 'group';
group.val     = {scans };
group.help    = {'Add a new group of scans to your experimental design'};
% ---------------------------------------------------------------------
% generic groups
% ---------------------------------------------------------------------
generic         = cfg_repeat;
generic.tag     = 'generic';
generic.name    = 'groups';
generic.help    = {''};
generic.values  = {group};
generic.val     = {group group group};
generic.num     = [3 Inf];

% ---------------------------------------------------------------------
% nullHypAllZero - Null Hypothesis: Groups are all equal? Or all zero?
% ---------------------------------------------------------------------
nullHypAllZero         = cfg_menu;
nullHypAllZero.tag     = 'nullHypAllZero';
nullHypAllZero.name    = 'Null Hypothesis';
nullHypAllZero.labels  = {'Groups are all equal' 'Groups are all zero'};
nullHypAllZero.values  = {true false};
nullHypAllZero.help    = {'','Null Hypothesis: Groups are all zero|all equal.'};


% % F Contrasts
% contr            = cfg_entry;
% contr.name       = 'F Contrasts';
% contr.tag        = 'Abtwns_contrast';
% contr.strtype    = 'i';
% contr.val        = {};
% contr.num        = [1 Inf];
% contr.help       = {'This is the number of Contrasts'};
% 
% % Covariate Values
% cv_none         = cfg_const;
% cv_none.tag     = 'cv_none';
% cv_none.name    = 'None';
% cv_none.val     = {1};
% cv_none.help    = {'Covariate value = none'};
% 
% cov_Val         = cfg_entry;
% cov_Val.tag     = 'cov_Val';
% cov_Val.name    = 'Covariate';%arbitary name
% cov_Val.help    = {'Help'};
% cov_Val.strtype = 'e';
% cov_Val.num     = [1 Inf];
% 
% cv_one         = cfg_branch;
% cv_one.tag     = 'cv_one';
% cv_one.name    = 'Enter Different Covariate Value';
% cv_one.val     = {cov_Val};
% cv_one.help    = {'Help'};
% 
% covariate         = cfg_choice;
% covariate.tag     = 'covariate';
% covariate.name    = 'Covariate'; %arbitary name
% covariate.val     = {cv_none };
% covariate.help    = {'Help'};
% covariate.values  = {cv_none cov_Val }; 



%% Executable Branch
snpmui = snpm_bch_ui(DesNm,DesFile,DesHelp,{generic nullHypAllZero}, true);
