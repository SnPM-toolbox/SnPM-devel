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

%% Executable Branch
snpmui = snpm_bch_ui(DesNm,DesFile,DesHelp,{CovInt xblock});
