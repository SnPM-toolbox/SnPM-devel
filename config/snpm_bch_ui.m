function [snpmui] = snpm_bch_ui(DesignName,DesignFile,DesignHelp,DesignItems,removeScansNode)
%  Set up common parameters to be set for all analyses.
% INPUT
%  DesignName  - Short, one-line description of design (string)
%  DesignFile  - Filename of the design configuration file (string)
%  DesignHelp  - Help to be displayed with this design (cell array of strings)
%  DesignItems - Design-specific menu items (cell array of Matlabbatch menu items)
% 
% OUTPUT
%  snpmui      - Design menu structure, with common menu items set before
%                and after design-specific ones.
%
%_______________________________________________________________________
% Thomas Nichols, Emma Thomas
% $Id$

% Based on Volkmar Glauche's MatlabBatch example code (Rev 1716)

% ---------------------------------------------------------------------
% Record internal information on the design
% ---------------------------------------------------------------------
DesNm         = cfg_const;
DesNm.tag     = 'DesignName';
DesNm.name    = 'DesignName';
DesNm.val     = {DesignName};
DesNm.help    = {['One-line descritption of the design.']};
DesNm.hidden  = true;

DesFile        = cfg_const;
DesFile.tag   = 'DesignFile';
DesFile.name   = 'DesignFile';
DesFile.val    = {DesignFile};
DesFile.help   = {['Filename of design configuration used.']};
DesFile.hidden = true;

DesignTag = spm_str_manip(DesignFile,'rt');
DesignTag = strrep(DesignTag,'snpm_bch_ui_','');


% ---------------------------------------------------------------------
% dir - Directory for analysis results
% ---------------------------------------------------------------------
dir         = cfg_files; % This is the generic data entry item
dir.name    = 'Analysis Directory'; % The displayed name
dir.tag     = 'dir';       % The name appearing in the harvested job structure. This name must be unique among all items in the val field of the superior node
dir.filter  = 'dir';
dir.ufilter = '.*';
dir.num     = [1 1];     % Number of inputs required (2D-array with exactly one row and one column)
dir.help    = {'','This sets the SnPM anlaysis directory.','All results will appear in this directory.'}; % help text displayed

% ---------------------------------------------------------------------
% scans Scans
% ---------------------------------------------------------------------
scans         = cfg_files; % This is the generic data entry item
scans.name    = 'Images to analyze'; % The displayed name
scans.tag     = 'P';       % The name appearing in the harvested job structure. This name must be unique among all items in the val field of the superior node
scans.filter  = 'image';      
scans.num     = [1 Inf];     % Number of inputs required (2D-array with exactly one row and one column)
scans.help    = {'','Model will be fit to these images.'}; % help text displayed

% SnPM specific options
% ---------------------------------------------------------------------

%Number of Permutations
nPerm            = cfg_entry;
nPerm.name       = 'Number of Permutations';
nPerm.tag        = 'nPerm';
nPerm.strtype    = 'e';
nPerm.def        = @(val)snpm_get_defaults('nPerm', val{:});
nPerm.num        = [1 1];
nPerm.help       = {'',...
		    'This sets the maximum number of permutations to use.  Actual number of permutations used may be smaller depending on the nubmer of scans and the design.',...
		    'The recommended number of permutations is 10,000.  If the number of possible permutations is smaller than the maximum, all possible permutations will be computed and the test will be "exact".  If there are more possible permutations than the maximum set, a "Monte Carlo Permutation Test" will be performed, and a random subset of permutations will be used.',...
		    'Stability of Monte Carlo P-values.  If the number of possible permutations is very large, the Monte Carlo Permutation Test P-values will have a standard deviation (over multiple runs of SnPM with different subsets of permutations) of sqrt(p*(1-p)/n), where p is the exact P-value (obtained from computing every possible permutation) and n is the number of permutations used.',...
		    'Thus, if the true P-value is 0.05 and only, say, 100 permutations are computed, the Monte Carlo standard deviation is 0.0218 and the 95% CI is approximately +/- 0.0436, almost equal to the true P-value itself.  To reduce the Monte Carlo confidence interval to one tenth of P=0.05, 7,500 permutations are required; 0,000 gives a 95% CI of +/- 0.0044.',...
		    'The use of a tiny number of permutations (e.g. 50 or 100) is recommended for testing large anlayses and to get an idea of how long an analysis will take to run.  Just realize the results will likely change with more permutations.'};

% ---------------------------------------------------------------------
% vFWHM - Variance FWHM smoothing
% ---------------------------------------------------------------------
vFWHM         = cfg_entry; % This is the generic data entry item
vFWHM.name    = 'Variance smoothing'; % The displayed name
vFWHM.tag     = 'vFWHM';       % The name appearing in the harvested job structure. This name must be unique among all items in the val field of the superior node
vFWHM.num     = [1 3];     % Number of inputs required (3-vector)
vFWHM.def     = @(val)snpm_get_defaults('vFWHM', val{:});
vFWHM.help    = {'','FWHM in mm smoothing applied to the variance image; if non-zero, this produces a "Pesudo T" image.',...
		 'With fewer than 20 degrees of freedom power can be improved with some variance smoothing.',...
		 'Try setting FWHM equal to the original smoothing applied to the data.'}; % help text displayed

% ---------------------------------------------------------------------
% bVolm - Load all data at once? Or work in 2D?
% ---------------------------------------------------------------------
bVolm         = cfg_menu;
bVolm.tag     = 'bVolm';
bVolm.name    = 'Memory usage';
bVolm.help    = {''};
bVolm.labels  = {'High' 'Low'};
bVolm.values  = {true false};
bVolm.def     = @(val)snpm_get_defaults('bVolm', val{:});
bVolm.help    = {'','Choosing "High" will result in SnPM loading the entire 4D dataset into memory, which will produce the fastest possible results.',...
		 'Choosing "Low" will cause SnPM to only load up a single plane of data at a time, running somewhat slower.',...
		 'Working plane-by-plane allows you to work with larger datasets than in "High" mode, but prevents you from doing variance smoothing in the z-dimension and cluster inference.'};

% ---------------------------------------------------------------------
% tm_none None
% ---------------------------------------------------------------------


% ---------------------------------------------------------------------
% bST - Collect Suprathreshold Statistics
% ---------------------------------------------------------------------
ST_none         = cfg_const;
ST_none.tag     = 'ST_none';
ST_none.name    = 'None';
ST_none.val     = {0};
ST_none.help    = {'No cluster inference - only voxel-wise.'};

ST_later         = cfg_const;
ST_later.tag     = 'ST_later';
ST_later.name    = 'Yes (slow, may create huge SnPM_ST.mat file)';
ST_later.val     = {-1};
ST_later.help    = {
    'Cluster inference, choose cluster-forming threshold *post-analysis*.'
    ''
    'Choosing this option will cause SnPM to save all of the "Mountain Tops" of the statistic image, creating a possibly very large SnPM_ST.mat file.  This will also slow down computation somewhat.  However, it allows you to chose the cluster-forming threshold after the computation has completed'
		   };

ST_now         = cfg_entry;
ST_now.tag     = 'ST_U';
ST_now.name    = 'Yes, set cluster-forming threshold now (fast)';
ST_now.strtype = 'e';
ST_now.num     = [1 1];
ST_now.help    = {
    'Cluster inference, choose cluster-forming threshold *now*.'
    ''
    'To speed computation and avoid a large SnPM_ST.mat file, specify the cluster-forming threshold now (as an uncorrected P-value or statistic (T or F) value).  If you want change your cluster-forming threshold, however, you will need to re-compute the SnPM analysis.'
		 };
ST_now.def     = @(val)snpm_get_defaults('ST_U');

ST         = cfg_choice;
ST.tag     = 'ST';
ST.name    = 'Cluster inference';
ST.val     = {ST_none};
ST.values  = {ST_none ST_later ST_now};
ST.help    = {'','Cluster inference in SnPM requires a substantial addtional computational and (possible) disk space burden, and hence is not done by default.',...
	      'If performing cluster size inference, you can choose to set the cluster-forming threshold later, in which case a large amount of information is written to disk on each permutation.',...
	      'Alternatively, you can choose to set the cluster-forming threshold now (pre-computation); this requires no additional disk space.  However, if you decide you want a different cluster-forming threshold you have to completely re-run the analysis.'};



% Tradtional SPM inputs...
% ---------------------------------------------------------------------


% ---------------------------------------------------------------------
% gmsca Grand mean scaling
% ---------------------------------------------------------------------
gmsca         = cfg_menu;
gmsca.tag     = 'gmsca';
gmsca.name    = 'Grand mean scaling';
gmsca.help    = {
                 'This option is only used for PET data.'
                 ''
                 'Selecting YES will specify ''grand mean scaling by factor'' which could be eg. ''grand mean scaling by subject'' if the factor is ''subject''. '
                 ''
                 'Since differences between subjects may be due to gain and sensitivity effects, AnCova by subject could be combined with "grand mean scaling by subject" to obtain a combination of between subject proportional scaling and within subject AnCova. '
                 ''
}';
gmsca.labels = {
                'No'
                'Yes'
}';
gmsca.values = {0 1};
gmsca.val    = {0};
% ---------------------------------------------------------------------
% ancova ANCOVA
% ---------------------------------------------------------------------
ancova         = cfg_menu;
ancova.tag     = 'ancova';
ancova.name    = 'ANCOVA';
ancova.help    = {
                  'This option is only used for PET data.'
                  ''
                  'Selecting YES will specify ''ANCOVA-by-factor'' regressors. This includes eg. ''Ancova by subject'' or ''Ancova by effect''. These options allow eg. different subjects to have different relationships between local and global measurements. '
                  ''
}';
ancova.labels = {
                 'No'
                 'Yes'
}';
ancova.values = {0 1};
ancova.val    = {0};



% 
% Masking options
% ---------------------------------------------------------------------


% ---------------------------------------------------------------------
% tm_none None
% ---------------------------------------------------------------------
tm_none         = cfg_const;
tm_none.tag     = 'tm_none';
tm_none.name    = 'None';
tm_none.val     = {1};
tm_none.help    = {'No threshold masking'};
% ---------------------------------------------------------------------
% athresh Threshold
% ---------------------------------------------------------------------
athresh         = cfg_entry;
athresh.tag     = 'athresh';
athresh.name    = 'Threshold';
athresh.help    = {
                   'Enter the absolute value of the threshold.'
                   ''
}';
athresh.strtype = 'e';
athresh.num     = [1 1];
athresh.val     = {100};
% ---------------------------------------------------------------------
% tma Absolute
% ---------------------------------------------------------------------
tma         = cfg_branch;
tma.tag     = 'tma';
tma.name    = 'Absolute';
tma.val     = {athresh };
tma.help    = {
               'Images are thresholded at a given value and only voxels at which all images exceed the threshold are included. '
               ''
               'This option allows you to specify the absolute value of the threshold.'
               ''
}';
% ---------------------------------------------------------------------
% rthresh Threshold
% ---------------------------------------------------------------------
rthresh         = cfg_entry;
rthresh.tag     = 'rthresh';
rthresh.name    = 'Threshold';
rthresh.help    = {
                   'Enter the threshold as a proportion of the global value'
                   ''
}';
rthresh.strtype = 'e';
rthresh.num     = [1 1];
rthresh.val     = {0.8};
% ---------------------------------------------------------------------
% tmr Relative
% ---------------------------------------------------------------------
tmr         = cfg_branch;
tmr.tag     = 'tmr';
tmr.name    = 'Relative';
tmr.val     = {rthresh };
tmr.help    = {
               'Images are thresholded at a given value and only voxels at which all images exceed the threshold are included. '
               ''
               'This option allows you to specify the value of the threshold as a proportion of the global value. '
               ''
}';
% ---------------------------------------------------------------------
% tm Threshold masking
% ---------------------------------------------------------------------
tm         = cfg_choice;
tm.tag     = 'tm';
tm.name    = 'Threshold masking';
tm.val     = {tm_none };
tm.help    = {
              'Images are thresholded at a given value and only voxels at which all images exceed the threshold are included. '
              ''
}';
tm.values  = {tm_none tma tmr };
% ---------------------------------------------------------------------
% im Implicit Mask
% ---------------------------------------------------------------------
im         = cfg_menu;
im.tag     = 'im';
im.name    = 'Implicit Mask';
im.help    = {
              'An "implicit mask" is a mask implied by a particular voxel value. Voxels with this mask value are excluded from the analysis. '
              ''
              'For image data-types with a representation of NaN (see spm_type.m), NaN''s is the implicit mask value, (and NaN''s are always masked out). '
              ''
              'For image data-types without a representation of NaN, zero is the mask value, and the user can choose whether zero voxels should be masked out or not.'
              ''
              'By default, an implicit mask is used. '
              ''
}';
im.labels = {
             'Yes'
             'No'
}';
im.values = {1 0};
im.val    = {1};

% ---------------------------------------------------------------------
% em Explicit Mask
% ---------------------------------------------------------------------
em         = cfg_files;
em.tag     = 'em';
em.name    = 'Explicit Mask';
em.val     = {{''}};
em.help    = {
              'Explicit masks are other images containing (implicit) masks that are to be applied to the current analysis.'
              ''
              'All voxels with value NaN (for image data-types with a representation of NaN), or zero (for other data types) are excluded from the analysis. '
              ''
              'Explicit mask images can have any orientation and voxel/image size. Nearest neighbour interpolation of a mask image is used if the voxel centers of the input images do not coincide with that of the mask image.'
              ''
}';
em.filter = 'image';
em.ufilter = '.*';
em.num     = [0 1];
% ---------------------------------------------------------------------
% masking Masking
% ---------------------------------------------------------------------
masking         = cfg_branch;
masking.tag     = 'masking';
masking.name    = 'Masking';
masking.val     = {tm im em };
masking.help    = {
                   'The mask specifies the voxels within the image volume which are to be assessed. SPM supports three methods of masking (1) Threshold, (2) Implicit and (3) Explicit. The volume analysed is the intersection of all masks.'
                   ''
}';


% Global Computation
% ---------------------------------------------------------------------

% ---------------------------------------------------------------------
% g_omit Omit
% ---------------------------------------------------------------------
g_omit         = cfg_const;
g_omit.tag     = 'g_omit';
g_omit.name    = 'Omit';
g_omit.val     = {1};
g_omit.help    = {'Omit'};
% ---------------------------------------------------------------------
% global_uval Global values
% ---------------------------------------------------------------------
global_uval         = cfg_entry;
global_uval.tag     = 'global_uval';
global_uval.name    = 'Global values';
global_uval.help    = {
                       'Enter the vector of global values'
                       ''
}';
global_uval.strtype = 'e';
global_uval.num     = [Inf 1];
% ---------------------------------------------------------------------
% g_user User
% ---------------------------------------------------------------------
g_user         = cfg_branch;
g_user.tag     = 'g_user';
g_user.name    = 'User';
g_user.val     = {global_uval };
g_user.help    = {
                  'User defined  global effects (enter your own '
                  'vector of global values)'
}';
% ---------------------------------------------------------------------
% g_mean Mean
% ---------------------------------------------------------------------
g_mean         = cfg_const;
g_mean.tag     = 'g_mean';
g_mean.name    = 'Mean';
g_mean.val     = {1};
g_mean.help    = {
                  'SPM standard mean voxel value'
                  ''
                  'This defines the global mean via a two-step process. Firstly, the overall mean is computed. Voxels with values less than 1/8 of this value are then deemed extra-cranial and get masked out. The mean is then recomputed on the remaining voxels.'
                  ''
}';
% ---------------------------------------------------------------------
% globalc Global calculation
% ---------------------------------------------------------------------
globalc         = cfg_choice;
globalc.tag     = 'globalc';
globalc.name    = 'Global calculation';
globalc.val     = {g_omit };
globalc.help    = {
                   'This option is only used for PET data.'
                   ''
                   'There are three methods for estimating global effects (1) Omit (assumming no other options requiring the global value chosen) (2) User defined (enter your own vector of global values) (3) Mean: SPM standard mean voxel value (within per image fullmean/8 mask) '
                   ''
}';
%%%
%%%  Stupid SnPM problem demands that a global is always calculated - To fix later
%%%  For now disallow "omit" option
%%%
%globalc.values  = {g_omit g_user g_mean };
globalc.values  = {g_user g_mean };


% Grand Mean scaling
% ---------------------------------------------------------------------

% ---------------------------------------------------------------------
% gmsca_no No
% ---------------------------------------------------------------------
gmsca_no         = cfg_const;
gmsca_no.tag     = 'gmsca_no';
gmsca_no.name    = 'No';
gmsca_no.val     = {1};
gmsca_no.help    = {'No overall grand mean scaling'};
% ---------------------------------------------------------------------
% gmscv Grand mean scaled value
% ---------------------------------------------------------------------
gmscv         = cfg_entry;
gmscv.tag     = 'gmscv';
gmscv.name    = 'Grand mean scaled value';
gmscv.help    = {
                 'The default value of 50, scales the global flow to a physiologically realistic value of 50ml/dl/min.'
                 ''
}';
gmscv.strtype = 'e';
gmscv.num     = [Inf 1];
gmscv.val     = {50};
% ---------------------------------------------------------------------
% gmsca_yes Yes
% ---------------------------------------------------------------------
gmsca_yes         = cfg_branch;
gmsca_yes.tag     = 'gmsca_yes';
gmsca_yes.name    = 'Yes';
gmsca_yes.val     = {gmscv };
gmsca_yes.help    = {
                     'Scaling of the overall grand mean simply scales all the data by a common factor such that the mean of all the global values is the value specified. For qualitative data, this puts the data into an intuitively accessible scale without altering the statistics. '
                     ''
}';
% ---------------------------------------------------------------------
% gmsca Overall grand mean scaling
% ---------------------------------------------------------------------
gmsca         = cfg_choice;
gmsca.tag     = 'gmsca';
gmsca.name    = 'Overall grand mean scaling';
gmsca.val     = {gmsca_no };
gmsca.help    = {
                 'Scaling of the overall grand mean simply scales all the data by a common factor such that the mean of all the global values is the value specified. For qualitative data, this puts the data into an intuitively accessible scale without altering the statistics. '
                 ''
                 'When proportional scaling global normalisation is used each image is separately scaled such that it''s global value is that specified (in which case the grand mean is also implicitly scaled to that value). So, to proportionally scale each image so that its global value is eg. 20, select <Yes> then type in 20 for the grand mean scaled value.'
                 ''
                 'When using AnCova or no global normalisation, with data from different subjects or sessions, an intermediate situation may be appropriate, and you may be given the option to scale group, session or subject grand means separately. '
                 ''
}';
gmsca.values  = {gmsca_no gmsca_yes };

% Global Normalisation
% ---------------------------------------------------------------------

% ---------------------------------------------------------------------
% glonorm Normalisation
% ---------------------------------------------------------------------
glonorm         = cfg_menu;
glonorm.tag     = 'glonorm';
glonorm.name    = 'Normalisation';
glonorm.help    = {
                   'Global nuisance effects are usually accounted for either by scaling the images so that they all have the same global value (proportional scaling), or by including the global covariate as a nuisance effect in the general linear model (AnCova). Much has been written on which to use, and when. Basically, since proportional scaling also scales the variance term, it is appropriate for situations where the global measurement predominantly reflects gain or sensitivity. Where variance is constant across the range of global values, linear modelling in an AnCova approach has more flexibility, since the model is not restricted to a simple proportional regression. '
                   ''
                   '''Ancova by subject'' or ''Ancova by effect'' options are implemented using the ANCOVA options provided where each experimental factor (eg. subject or effect), is defined. These allow eg. different subjects to have different relationships between local and global measurements. '
                   ''
                   'Since differences between subjects may be due to gain and sensitivity effects, AnCova by subject could be combined with "grand mean scaling by subject" (an option also provided where each experimental factor is originally defined) to obtain a combination of between subject proportional scaling and within subject AnCova. '
                   ''
}';
glonorm.labels = {
                  'None'
                  'Proportional'
                  'ANCOVA'
}';
glonorm.values = {1 2 3};
glonorm.val    = {1};
% ---------------------------------------------------------------------
% globalm Global normalisation
% ---------------------------------------------------------------------
globalm         = cfg_branch;
globalm.tag     = 'globalm';
globalm.name    = 'Global normalisation';
globalm.val     = {gmsca glonorm };
globalm.help    = {
                   'This option is only used for PET data.'
                   ''
                   'Global nuisance effects are usually accounted for either by scaling the images so that they all have the same global value (proportional scaling), or by including the global covariate as a nuisance effect in the general linear model (AnCova). Much has been written on which to use, and when. Basically, since proportional scaling also scales the variance term, it is appropriate for situations where the global measurement predominantly reflects gain or sensitivity. Where variance is constant across the range of global values, linear modelling in an AnCova approach has more flexibility, since the model is not restricted to a simple proportional regression. '
                   ''
                   '''Ancova by subject'' or ''Ancova by effect'' options are implemented using the ANCOVA options provided where each experimental factor (eg. subject or effect), is defined. These allow eg. different subjects to have different relationships between local and global measurements. '
                   ''
                   'Since differences between subjects may be due to gain and sensitivity effects, AnCova by subject could be combined with "grand mean scaling by subject" (an option also provided where each experimental factor is originally defined) to obtain a combination of between subject proportional scaling and within subject AnCova. '
                   ''
}';

if ~exist('removeScansNode', 'var') || ~removeScansNode
    snpmui_pre = {DesNm DesFile dir scans};
else
    snpmui_pre = {DesNm DesFile dir};
end
%snpmui_pre = {dir scans};
snpmui_des = DesignItems; 
snpmui_post = {nPerm vFWHM bVolm ST masking globalc globalm };

%% Executable Branch
snpmui      = cfg_exbranch;  % This is the branch that has information about how to run this module
snpmui.name = DesignName;   % The display name
snpmui.tag  = DesignTag; % The name appearing in the harvested job structure. This name must be unique among all items in the val field of the superior node
snpmui.val  = {snpmui_pre{:} snpmui_des{:} snpmui_post{:}};    % The items that belong to this branch. All items must be filled before this branch can run or produce virtual outputs
snpmui.prog = @snpm_run_ui;
snpmui.help = DesignHelp;
snpmui.vout = @snpm_bch_ui_vout;

function vout = snpm_bch_ui_vout(job)
% Determine what outputs will be present if this job is run. In this case,
% the structure of the inputs is fixed, and the output is always a single
% number. Note that input items may not be numbers, they can also be
% dependencies.

vout = cfg_dep;                        % The dependency object
vout.sname      = 'SnPMcfg.mat configuration file';       % Displayed dependency name
vout.src_output = substruct('.','SnPMcfg'); % The output subscript reference. This could be any reference into the output variable created during computation

