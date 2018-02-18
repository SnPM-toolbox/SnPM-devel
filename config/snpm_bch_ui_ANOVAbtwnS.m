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
    'Create design and permutation matrix appropriate for k groups analyses where there is just *one* scan per subject.',...
    '',...
    'This PlugIn can be regarded as a generalization of a two-sample t-test to k>2 groups.  Instead of producing a t-statistic, it produces a F statistic.',...
    '',...
    'The PlugIn can test for two different types of effects.  It can test for the presence of any *non-zero* effect among the k groups; that is, it tests the null hypothesis that all of the groups are mean zero.  Or it can test for the presence of any differences *between* the k groups; that is, it tests the null hypothesis that all of the groups have some (possibly non-zero) common mean.',...
    '',...
    'Number of permutations. There are (nScan)!/(GrpCnt[1]!*GrpCnt[2]!*...*GrpCnt[k]!) possible permutations, where nScan is the total number of scans and GrpCnt(i) is the size of the ith group. round(exp(gammaln(nScan+1)-sum(gammaln(GrpCnt+1))))',...
	  };

% ---------------------------------------------------------------------
% scans 
% ---------------------------------------------------------------------
scans         = cfg_files;
scans.tag     = 'scans';
scans.name    = 'Scans';
scans.help    = {'Select images for this group. '};
scans.filter = 'image';
scans.ufilter = '.*';
scans.num     = [1 Inf];
% ---------------------------------------------------------------------
% group
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
nullHypAllZero.values  = {false true};
nullHypAllZero.help    = {'','Null Hypothesis: Groups are all zero|all equal.'};


%% Executable Branch
snpmui = snpm_bch_ui(DesNm,DesFile,DesHelp,{generic nullHypAllZero}, true);
