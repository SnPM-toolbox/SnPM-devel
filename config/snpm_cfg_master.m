function cfg = cfg_example_master
% Master file that collects the cfg_exbranches in conceptually similar
% groups.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_example_master.m 1716 2008-05-23 08:18:45Z volkmar $

rev = '$Rev: 1716 $'; %#ok


%%%%
%%%%  Later... maybe re-write this to be like (e.g.) spm_bch_realign,
%%%%  such that common snpm_bch_ui components are setup up *here*, along
%%%%  with the design-specific aspects.  The "pre" aspects I can put
%%%%  here, and the "post" aspects can stay in snpm_bch_ui (or I can
%%%%  create two functions).  This will simplify the snpm_bch_ui_* files
%%%%  

des        = cfg_choice;
des.name   = 'Specify';
des.tag    = 'Design';
des.values = {snpm_bch_ui_TwoSampTss...
	      snpm_bch_ui_Corr1S  ...
	      snpm_bch_ui_OneSampT  ...
	      snpm_bch_ui_Corr ...
	      snpm_bch_ui_PairT...
          ... % 	      snpm_bch_ui_PairTrand  ... 
	      snpm_bch_ui_ANOVAwithinS  ...
	      snpm_bch_ui_TwoSampPairT  ...
	      snpm_bch_ui_TwoSampT  ...
	      snpm_bch_ui_ANOVAbtwnS};
des.help   = {['These are the different models that are available to run with SnPM']};

%% Collect above Collections
cfg        = cfg_choice;
cfg.name   = 'SnPM';
cfg.tag    = 'cfg_snpm';
cfg.values = {des snpm_bch_cp snpm_bch_pp}; % Values in a cfg_repeat can be any cfg_item objects
cfg.help   = {'SnPM - In batch mode'};
