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
    'Create design and permutation matrix appropriate for one group analyses where there are multiple scans per subject, and where each scan is itself a difference image or contrast image.  This plug in effects a within subject ANOVA.',...
    '',...
    'A common use of this PlugIn is for an F test for a set of contrasts. For each subject, we have k contrasts jointly expressing some effect of interest.  Under the null hypothesis we assume that the data for each subject is unpeturbed by multiplication by -1.  That is, under the null hypothesis the multivariate measurements are all mean zero and symmetrically distributed.  We assume exchangeability between subjects (just as we usually assume independent subjects) but *do* *not* assume that the k values for each subject are independent.',...
    '',...
    'The PlugIn tests for the presence of *any* effect among the k contrasts.  That is, it tests the null hypothesis that all of the effects are mean zero.',...
    '',...
    'Number of permutations. There are 2^(nSubj-1) possible permutations, where nSubj is the total number of subjects. Intuitively, each subject can be assigned to +1 or -1, so we should have 2^nSubj possible permutations. However, since we are doing an F test and all +1''s and all -1''s would give us the same F statistic. To avoid the redundance, therefore we explicitly assign the first subject to +1 group.',...
    '',...
    'It is recommended that at least 7 or 8 subjects are used; with only 6 subjects, the permutation distribution will only have 2^5 = 32 elements and the smallest p-value will be 1/32=0.03125.',...
	  };

%% Questions
% ---------------------------------------------------------------------
% scans Scans
% ---------------------------------------------------------------------
scans         = cfg_files;
scans.tag     = 'scans';
scans.name    = 'Scans';
scans.help    = {'Select the images to be analysed.  They must all have the same image dimensions, orientation, voxel size etc.'};
scans.filter = 'image';
scans.ufilter = '.*';
scans.num     = [1 Inf];
% ---------------------------------------------------------------------
% fsubject Subject
% ---------------------------------------------------------------------
fsubject         = cfg_branch;
fsubject.tag     = 'fsubject';
fsubject.name    = 'Subject';
fsubject.val     = {scans};
fsubject.help    = {'Enter data and conditions for a new subject'};
% ---------------------------------------------------------------------
% generic Subjects
% ---------------------------------------------------------------------
generic         = cfg_repeat;
generic.tag     = 'generic';
generic.name    = 'Subjects';
generic.help    = {''};
generic.values  = {fsubject };
generic.num     = [1 Inf];


% % Scans per Subject
% scans_sub            = cfg_entry;
% scans_sub.name       = 'Scans per subject';
% scans_sub.tag        = 'Anovawithin_scans';
% scans_sub.strtype    = 'i';
% scans_sub.val        = {};
% scans_sub.num        = [1 1];
% scans_sub.help       = {'This is the number of subjects'};

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
cov         = cfg_branch;
cov.tag     = 'cov';
cov.name    = 'Covariate';
cov.val     = {c cname };
cov.help    = {'Add a new covariate to your experimental design'};
% ---------------------------------------------------------------------
% generic Covariates
% ---------------------------------------------------------------------
generic_cov         = cfg_repeat;
generic_cov.tag     = 'generic_cov';
generic_cov.name    = 'Covariates of no interest';
generic_cov.help    = {'Covariates of no interest'};
generic_cov.values  = {cov };
generic_cov.num     = [0 Inf];

%% Executable Branch
snpmui = snpm_bch_ui(DesNm,DesFile,DesHelp,{generic generic_cov}, true);
