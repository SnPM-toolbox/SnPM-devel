function snpmui = snpm_bch_ui_TwoSampTss
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

DesNm = 'SingleSub: Two Sample T test; 2 conditions, replications';
DesFile = mfilename;
DesHelp = {'',...
        'Create design and permutation matrix appropriate for single subject, two condition activation (with replication) experiments.',...
        '',...
        'Number of permutations. There are nScan-choose-nRepl possible permutations, where nScan is the number of scans and nRepl is the number of replications (nScan = 2*nRepl). In Matlab, use the nchoosek command to check the number of permutations, as in nchoosek(nScan,nRepl).',...
        ...%Matlab doesn''t have a choose function but you can use this expression prod(1:nScan)/prod(1:nRepl)^2',...
      };

%% Questions

% Number of Replications 
ReplC            = cfg_entry;
ReplC.name       = 'Replications of Conditions';
ReplC.tag        = 'Tss_repc';
ReplC.strtype    = 'e';
ReplC.val        = {};
ReplC.num        = [1 1];
ReplC.help       = {'This is the Number of replications'}; 

% We need to define the order of the conditions
cond            = cfg_entry;
cond.name       = 'Condition index';
cond.tag        = 'condidx';
cond.strtype    = 'i';
cond.val        = {};
cond.num        = [1 Inf];
cond.help = {'Enter conditions indices: (1/2)'};

% Exchangability Blocks
ex_blk            = cfg_entry;
ex_blk.name       = 'Size of Exchangability Block';
ex_blk.tag        = 'TwosampTss_Block';
ex_blk.strtype    = 'i';
ex_blk.val        = {};
ex_blk.num        = [1 1];
ex_blk.help       = {'This is the number of subjects'};

%% Executable Branch
snpmui = snpm_bch_ui(DesNm,DesFile,DesHelp,{ReplC cond ex_blk});
