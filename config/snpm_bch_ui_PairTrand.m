function snpmui = snpm_bch_ui_PairTrand
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

DesNm = 'MultiSub: 2 conditions, replications - permutation test';
DesFile = mfilename;
DesHelp = {'',...
	  'stuff.',...
	  '',...
	  'stuff ',...
	  '    stuff',...
	  'stuff.',...
	  };
  
%% Questions


% Number of Subjects
Nsub            = cfg_entry;
Nsub.name       = 'Enter Number of Subjects';
Nsub.tag        = 'PairTrand_Subjects';
Nsub.strtype    = 'e';
Nsub.val        = {};
Nsub.num        = [1 1];
Nsub.help       = {'This is the covariate value'}; 

% Number of Replications
Repcond            = cfg_entry;
Repcond.name       = 'Replications of Conditions';
Repcond.tag        = 'PairTrand_Rep';
Repcond.strtype    = 'e';
Repcond.val        = {};
Repcond.num        = [1 1];
Repcond.help       = {'This is the Number of replications'}; 


%% Executable Branch
snpmui = snpm_bch_ui(DesNm,DesFile,DesHelp,{Nsub Repcond});
