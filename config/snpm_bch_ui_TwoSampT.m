function snpmui = snpm_bch_ui_TwoSampT
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

DesNm = '2 Groups: Two Sample T test; 1 scan per subject';
DesFile = mfilename;
DesHelp = {'Create design and permutation matrix appropriate for two group analyses where there is just *one* scan per subject.', ...
            '', ...
            'Keep in mind that when only 1 scan per subject is used there is no way to control for anatomical differences, hence the differences identified will be attributable to both functional and anatomical differences between the groups.  ', ...
            '', ...
            'A common source of between group anatomical differences is age; older subjects tend to have larger ventricles and thinner gray matter relative to younger subjects.  One approach to address this difference is to include a linear confounding covariate of age; alternatively, a dichotomous covariate (consisting of just 0''s and 1''s) indicating young/old can be used.  Including such covariates will ensure that group differences are not atributable to linear (or constant, for 0/1 covariate) effects of age.  Hence, with both of these approaches if ages are not equally distributed between groups then including such covariates can reduce the signal attributable to group differences, since some of the signal could just be due to age.',...
            '',...
	  };

% Group memberships
group_memb         = cfg_entry;
group_memb.tag     = 'group_memb';
group_memb.name    = 'Group membership';%Name displayed in the batch tree
group_memb.help    = {'Enter subject index: (A/B)'};
group_memb.strtype = 's'; % string
group_memb.num     = [1 Inf]; % Expected format  
  
% % Number of Covariates
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
% covariate.name    = 'Covariate Value'; %arbitary name
% covariate.val     = {cv_none };
% covariate.help    = {'Help'};
% covariate.values  = {cv_none cov_Val };

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
snpmui = snpm_bch_ui(DesNm,DesFile,DesHelp,{group_memb generic});
