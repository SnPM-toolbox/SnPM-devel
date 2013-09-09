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
	  'stuff.',...
	  '',...
	  'stuff ',...
	  '    stuff',...
	  'stuff.',...
	  };
  
%% Questions

%Number of Subjects
N_sub            = cfg_entry;
N_sub.name       = 'Enter Number of Subjects';
N_sub.tag        = 'N_sub';
N_sub.strtype    = 'e';
N_sub.val        = {};
N_sub.num        = [1 1];
N_sub.help       = {'Number of subjects'}; 

%Replications of Conditions
rep_cond            = cfg_entry;
rep_cond.name       = 'Replications of Conditions';
rep_cond.tag        = 'rep_cond';
rep_cond.strtype    = 'e';
rep_cond.val        = {};
rep_cond.num        = [1 1];
rep_cond.help       = {'This is the covariate value'}; 

%% Executable Branch
snpmui = snpm_bch_ui(DesNm,DesFile,DesHelp,{N_sub rep_cond});
