function snpmui = snpm_bch_ui_Corr1s
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
DesHelp = {'Create design and permutation matrix appropriate for single-subject, correlation design.'};


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

%Covariate Values
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
snpmui = snpm_bch_ui(DesNm,DesFile,DesHelp,{CovInt xblock covariate});
