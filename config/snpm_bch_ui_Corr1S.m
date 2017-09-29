function snpmui = snpm_bch_ui_Corr1S
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

 
DesNm = 'SingleSub: Simple Regression (correlation); single covariate of interest';
DesFile = mfilename;
DesHelp = {'',...
    'Create design and permutation matrix appropriate for single-subject, correlation design.',...
    '',...
    'Number of permutations. There are nScan! (nScan factoral) possible permutations, where nScan is the number of scans of the subject. You can compute this using the gamma function in Matlab: nScan! is gamma(nScan+1); or by direct computation as prod(1:nScan)',...
    };


%% Questions

%Enter Covariate Value
CovInt            = cfg_entry;
CovInt.name       = 'Covariate';
CovInt.tag        = 'CovInt';
CovInt.strtype    = 'e';
CovInt.val        = {};
CovInt.num        = [1 Inf];
CovInt.help       = {'This is the variable to correlate with the imaging data.'}; 

%Size of echangability block
xblock            = cfg_entry;
xblock.name       = 'Size of Exchangability Block';
xblock.tag        = 'xblock';
xblock.strtype    = 'i';
xblock.val        = {};
xblock.num        = [1 1];
xblock.help       = {'This is the number of exhangability blocks'};

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
generic.name    = 'Covariates of no interest';
generic.help    = {'Covariates of no interest'};
generic.values  = {cov };
generic.num     = [0 Inf];

%% Executable Branch
snpmui = snpm_bch_ui(DesNm,DesFile,DesHelp,{CovInt generic xblock});
