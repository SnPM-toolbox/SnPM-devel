function snpmui = snpm_bch_ui_OneSampT
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
% $Id: Thomas Nichols, Emma Thomas, Camille Maumet $


snpm_defaults

rev = '$Rev: 1716 $';

DesNm = 'MultiSub: One Sample T test on diffs/contrasts';
DesFile = mfilename;
DesHelp = {'',...
    'Create design and permutation matrix appropriate for one group analyses where there is just *one* scan per subject.  This plug in effects a one-sample t-test.',...
    '',...
    'A common use of this plug is for random effects analysis of contrast images.  For this analysis we only need to assume, under the null hypothesis,  that each of the images are exchangeble and the contrast images have mean zero, symmetrically distributed data at each voxel. (Exchangeability follows from independence of different subjects.)',...
    '',...
    'Number of permutations. There are 2^nSubj possible permutations, where nScan is the total number of scans. It is recommended that at least 6 or 7 subjects are used; with only 5 subjects, the permutation distribution will only have 2^5 = 32 elements and the smallest p-value will be 1/32=0.03125.',...
	  };
  
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
generic         = cfg_repeat;
generic.tag     = 'generic';
generic.name    = 'Covariates';
generic.help    = {'Covariates'};
generic.values  = {cov };
generic.num     = [0 Inf];


%% Executable Branch
snpmui = snpm_bch_ui(DesNm,DesFile,DesHelp,{generic});

