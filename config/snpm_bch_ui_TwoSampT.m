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
DesHelp = { '',...
            'Create design and permutation matrix appropriate for two group analyses where there is just *one* scan per subject.', ...
            '', ...
            'Keep in mind that when only 1 scan per subject is used there is no way to control for anatomical differences, hence the differences identified will be attributable to both functional and anatomical differences between the groups.  ', ...
            '', ...
            'A common source of between group anatomical differences is age; older subjects tend to have larger ventricles and thinner gray matter relative to younger subjects.  One approach to address this difference is to include a linear confounding covariate of age; alternatively, a dichotomous covariate (consisting of just 0''s and 1''s) indicating young/old can be used.  Including such covariates will ensure that group differences are not atributable to linear (or constant, for 0/1 covariate) effects of age.  Hence, with both of these approaches if ages are not equally distributed between groups then including such covariates can reduce the signal attributable to group differences, since some of the signal could just be due to age.',...
            '',...
            'Number of permutations. There are nTot-choose-nGrp possible permutations, where nTot is the total number of subjects (and scans) and nGrp is the size of one of the groups.'
	  };

  % Similar to SPM two-sample t-test
% ---------------------------------------------------------------------
% scans1 Group 1 scans
% ---------------------------------------------------------------------
scans1         = cfg_files;
scans1.tag     = 'scans1';
scans1.name    = 'Group 1 scans';
scans1.help    = {'Select the images from sample 1.  They must all have the same image dimensions, orientation, voxel size etc.'};
scans1.filter = 'image';
scans1.ufilter = '.*';
scans1.num     = [1 Inf];
% ---------------------------------------------------------------------
% scans2 Group 2 scans
% ---------------------------------------------------------------------
scans2         = cfg_files;
scans2.tag     = 'scans2';
scans2.name    = 'Group 2 scans';
scans2.help    = {'Select the images from sample 2.  They must all have the same image dimensions, orientation, voxel size etc.'};
scans2.filter = 'image';
scans2.ufilter = '.*';
scans2.num     = [1 Inf];  
  
% % Group memberships
% group_memb         = cfg_entry;
% group_memb.tag     = 'group_memb';
% group_memb.name    = 'Group membership';%Name displayed in the batch tree
% group_memb.help    = {'Enter subject index: (A/B)'};
% group_memb.strtype = 's'; % string
% group_memb.num     = [1 Inf]; % Expected format  

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
% cov Covariate
% ---------------------------------------------------------------------
cov         = cfg_branch;
cov.tag     = 'cov';
cov.name    = 'Covariate';
cov.val     = {c cname };
cov.help    = {'Add a new covariate to your experimental design'};
% ---------------------------------------------------------------------
% generic Covariates
% ---------------------------------------------------------------------
generic         = cfg_repeat;
generic.tag     = 'generic';
generic.name    = 'Covariates';
generic.help    = {'Covariates'};
generic.values  = {cov };
generic.num     = [0 Inf];

%% Executable Branch
snpmui = snpm_bch_ui(DesNm,DesFile,DesHelp,{scans1, scans2, generic}, true);
