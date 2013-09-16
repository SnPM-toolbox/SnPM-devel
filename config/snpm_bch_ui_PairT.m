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

%% Executable Branch
snpmui = snpm_bch_ui(DesNm,DesFile,DesHelp,{ sub_files }, true);
