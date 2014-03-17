function snpmui = snpm_bch_ui_PairT
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

DesNm = 'MultiSub: Paired T test; 2 conditions, 1 scan per condition';
DesFile = mfilename;
DesHelp = {'',...
    'Create design and permutation matrix appropriate for multi-subject, two condition no replication design, where the condition labels have been applied to the subjects in a nonrandomized fashion.',...
    '',...
    'With only two scans per subject, there are only two possible sets of labels: AB and BA.  No restriction is placed on how many subjects received A first.  That is, unlike snpm_MSA2x, this PlugIn handles designs where all subjects receive condition A first.',...
    '',...
    'Since there is no randomization we must justify exchangibility with assumptions under the null hypothesis.  In particular, we must assume that at each voxel, the distribution of the data is the same for A and B scans and for all subjects.  This is equivalent to assuming that distribution of A-B is symmetric and the same for all subjects. Note that that an unmodeled temporal effect would violated this assumption (no such assumption is needed when a randomization design is used; if randomization was used snpm_MSA2x should be used).',...
    '',...
    '(If there are replications, it is recommended that a "first level") (model is fit, reducing each condition to a single summary image. )',...
    '',...
    'Number of permutations. There are 2^nSubj possible labelings, where nSubj is the number of subjects. This is because each subject can have two possible states,flipped or unflipped. For example, 2^7 = 128 and 2^8 = 256.  Hence at least eight subjects are needed to characterize the permutation distribution well, and 10 or more are best.',...
  };
  
%% Questions

% Reuse flexible factorial component of SPM to get the design files and
% associated conditions
factorial_design = spm_cfg_factorial_design();
sub_files = factorial_design.val{2}.values{end}.val{2}.val{1};

% We expect exactly two files per subject
sub_files.values{1}.val{1}.num = [2 2];
sub_files.values{1}.val{1}.help = {'Select scans'};
% We need to define the order of the two conditions
sub_files.values{1}.val{2}.num = [1 2];
sub_files.values{1}.val{2}.name = 'Scan index';
sub_files.values{1}.val{2}.help = {'Enter scan index: (1 2|2 1)'};
sub_files.values{1}.val{2}.tag = 'scindex';
sub_files.val{1} = sub_files.values{1};

%% Executable Branch
snpmui = snpm_bch_ui(DesNm,DesFile,DesHelp,{ sub_files }, true);
