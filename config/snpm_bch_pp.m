function snpmpp = snpm_bch_pp
% Setup menus to draw inference in an SnPM analysis
%_______________________________________________________________________
% Copyright (C) 2013 The University of Warwick
% Id: snpm_bch_cp.m  SnPM13 2013/10/12
% Thomas Nichols, Camille Maumet
%
% Based on Volkmar Glauche's MatlabBatch example code (Rev 1716)

% Common items for all Results actions
%
% Input SnPM results mat file
snpmres         = cfg_files;
snpmres.tag     = 'SnPMmat';
snpmres.name    = 'SnPM.mat results file';
snpmres.help    = {'Select a SnPM.mat results file.'};
snpmres.filter = 'any';
snpmres.ufilter = 'SnPM.mat';
snpmres.num     = [1 1];

%
% "Results"
%

Tth        = cfg_entry;
Tth.tag    = 'Tth';
Tth.name   = 'Uncorrected T or F (statistic threshold)';
Tth.strtype= 'e';
Tth.num    = [1 1];
Tth.help   = {'Enter a statistic value (uncorrected) threshold.  The threshold will be interpreted as T or F depending on the design used.'};

CFth        = cfg_entry;
CFth.tag    = 'CFth';
CFth.name   = 'Cluster-Forming Threshold';
CFth.strtype= 'e';
CFth.num    = [1 1];
CFth.val    = {NaN};
CFth.help   = {
    'Enter an uncorrected P-value or statistic value used to define clusters (values strictly less than 1 are taken to be P-values).'
    'If you configured "Cluster Inference: Yes ... (fast)" option, leave this as NaN, as the threshold has already been set.  Otherwise ("Cluster Inference: Yes (slow...)", this *must* be set in order to perform cluster inference.'
    'Note, with variance smoothing (''pseduo t'') there are no uncorrected P-values available, and hence the cluster forming threshold must be specfied as greater than 1.'};
 
Cth         = cfg_entry;
Cth.tag     = 'Cth';
Cth.name    = 'Uncorrected k';
Cth.strtype = 'e';
Cth.num     = [1 1];
Cth.help    = {'Enter a cluster size threshold (uncorrected; all clusters with strictly fewer voxels will be discarded).'};

PthC         = cfg_entry;
PthC.tag     = 'PthC';
PthC.name    = 'Uncorrected Nonparametric P(k)';
PthC.strtype = 'e';
PthC.num     = [1 1];
PthC.help    = {'Enter an P-value threshold to obtain uncorrected cluster inference.'
		'WARNING 1.  This produces P-values *only* for an a priori selected cluster (compare to SPM''s uncorrected P(k)).'
		'WARNING 2.  These P-values are obtained with an assumption of stationarity which is likely to be inappropriate for VBM, M/EEG and other types of data.  Specfically, the permutation P-value is computed by  pooling all cluster sizes over space and permutations, which is only valid if cluster sizes arising over the brain are homogeneous.'};

Pth         = cfg_entry;
Pth.tag     = 'Pth';
Pth.name    = 'Uncorrected Nonparametric P';
Pth.strtype = 'e';
Pth.num     = [1 1];
Pth.help    = {'Enter an uncorrected P-value threshold to applied to the nonparametric P-value image.  Note that, because SnPM conducts a permutation test at each voxel, there is not a 1-to-1 mapping between T/F and P.  Consider choosing a voxel-wise threshold on the T (or F) image that marks 100 voxels as significant; now, choose a voxel-wise threshold on the nonparametric P image that also marks 100 voxels significiant.  Because of the voxel-spefic permutaiton P-values, the 100 voxels identified in each case will be different.  '
	       ''
	       'The nonparametric P image is "safer" because the null distribution is determined individually at every voxel, but it may be noisier if there are a small number of permutations.  Inference based on the T (or F) image implicitly assumes that the statistic values mean the same thing over the image; FWE-corrected inferences are always valid even if the true null distribution varies over the image, but specificity may vary voxel-by-voxel.'};

TFth         = cfg_entry;
TFth.tag     = 'TFth';
TFth.name    = 'Uncorrected T or F';
TFth.strtype = 'e';
TFth.num     = [1 1];
TFth.help    = {'Enter an uncorrected T or F threshold'};

FWEth        = cfg_entry;
FWEth.tag    = 'FWEth';
FWEth.name   = 'FWE Corrected';
FWEth.strtype= 'e';
FWEth.num    = [1 1];
FWEth.def    = @(val)snpm_get_defaults('FWElevel', val{:});
FWEth.help   = {'Enter a FWE significance level threshold.  This threshold is applied to a T (or F) image.'};

FWEthC        = cfg_entry;
FWEthC.tag    = 'FWEthC';
FWEthC.name   = 'FWE Corrected';
FWEthC.strtype= 'e';
FWEthC.num    = [1 1];
FWEthC.def    = @(val)snpm_get_defaults('FWElevel', val{:});
FWEthC.help   = {'Enter a FWE significance level threshold.  This threshold is applied to clusters formed on a T (or F) image.'};

FDRth        = cfg_entry;
FDRth.tag    = 'FDRth';
FDRth.name   = 'FDR Corrected';
FDRth.strtype= 'e';
FDRth.num    = [1 1];
FDRth.def    =  @(val)snpm_get_defaults('FDRlevel', val{:});
FDRth.help   = {'Enter a FDR significance level threshold.  This threshold is applied to a nonparametric P-value image.'};

VoxSig         = cfg_choice;
VoxSig.name    = 'Significance Level';
VoxSig.tag     = 'VoxSig';
VoxSig.values  = {Pth TFth FDRth FWEth};
VoxSig.val     = {FWEth};
VoxSig.help    = {'Select voxel-wise significance thresholding method.  You can choose between uncorrected and corrected inference, but note the different type of images they apply to.'
    '"Uncorrected Nonparametric P" thresholds to the result of permutation tests applied at each voxel.' 
    '"Uncorrected T or F" thresholds the T (or F) statistic image.'
    '"FDR Corrected" thresholds the nonparametric P image to control the False Discovery Rate, the expected proportion of false positive voxels among detected voxels.'
    '"FWE Corrected" thresholds the T (or F) statistic image to control the Familywise Error Rate, the chance of one or more false positive voxels'
    ''
    'Note, in standard parametric statistics, there is a fixed 1-to-1 mapping between T (or F) values and P-values, and thus it is irrelevant whether you make inferences base on an image of statistic values or P-values.  In SnPM, a nonpararmetric permutation test is conducted at each and every voxel, making the mapping between T (or F) values and P-values voxel-specific.  In particular, there is no equivalent threshold on the T and nonparametric P images that will produce the same set of significant voxels.'
		 };

ClusSig         = cfg_choice;
ClusSig.name    = 'Significance Level';
ClusSig.tag     = 'ClusSig';
ClusSig.values  = {Cth PthC FWEthC};
ClusSig.val     = {FWEthC};
ClusSig.help    = {'Select cluster-wise significance thresholding method.  You can choose between uncorrected and corrected inference.'
    '"Uncorrected k" removes all clusters smaller than the specified threshold, expressed in units of voxels.' 
    '"Uncorrected Nonparametric P(k)" determines a cluster size threshold with an uncorrected P-value; see the help for this options for important limitations on its use.'
    '"FWE Corrected" removes all clusters smaller than a threshold that will control the Familywise Error Rate, the chance of one or more false positive clusters.'
    ''
    'Note, clusters are defined on T (or F) images, and not on the nonparametric P-value image.'
		 };
     
Vox = cfg_branch;
Vox.name = 'Voxel-Level Inference';
Vox.tag = 'Vox';
Vox.val = {VoxSig};
Vox.help = {'Specify voxel-level inference'};

ClusSize = cfg_branch;
ClusSize.name = 'Cluster size';
ClusSize.tag = 'ClusSize';
ClusSize.val = {CFth ClusSig};
ClusSize.help = {'Specify cluster-size inference'};

PFilt         = cfg_entry;
PFilt.tag     = 'PFilt';
PFilt.name    = 'Corrected p-value for filtering';
PFilt.strtype = 'e';
PFilt.num     = [1 1];
PFilt.val    = {0.05};
PFilt.help    = {'Select corrected threshold'};

PrimThresh         = cfg_entry;
PrimThresh.tag     = 'PrimThresh';
PrimThresh.name    = 'Primary threshold';
PrimThresh.strtype = 'e';
PrimThresh.num     = [1 1];
PrimThresh.help    = {'Select primary threshold for STC analysis.'};

theta         = cfg_entry;
theta.tag     = 'Theta';
theta.name    = 'Theta';
theta.strtype = 'e';
theta.num     = [1 1];
theta.val    = {0.5};
theta.help    = {'Theta value for voxel-cluster combining'};

ClusMass = cfg_branch;
ClusMass.name = 'Cluster mass';
ClusMass.tag = 'ClusMass';
ClusMass.val = {PFilt PrimThresh theta};
ClusMass.help = {'Specify cluster-mass inference'};



Clus = cfg_choice;
Clus.name = 'Cluster-Level Inference';
Clus.tag = 'Clus';
Clus.values = {ClusSize ClusMass};
Clus.val = {ClusSize};
Clus.help = {'Select cluster-level statistic'};

ThrType        = cfg_choice;
ThrType.name   = 'Type of Thresholding';
ThrType.tag    = 'Thr';
ThrType.values = {Vox Clus};
ThrType.val    = {Vox};
ThrType.help   = {
    'Choose between voxel-wise and cluster-wise inference.'
    };


%
% Menus
%
posneg        = cfg_menu;
posneg.name   = 'Display positive or negative effects?';
posneg.tag    = 'Tsign';
posneg.labels = {'Positives','Negatives'};
posneg.values = {1,-1};
posneg.val    = {1};
posneg.help   = {
    'For models that produce T images, select the direction of the effect to display.'
    ''
    'Ignored for F statistic images.'		};

WF_no      = cfg_const;
WF_no.tag  = 'WF_no';
WF_no.name = 'No';
WF_no.val  = {0};
WF_no.help = {'Do not write any filtered statistic image.'};

WF_name      = cfg_entry;
WF_name.tag  = 'name';
WF_name.name = 'Image name to write';
WF_name.strtype  = 's';
WF_name.val  = {'SnPM_filtered'};
WF_name.help = {'Write filtered statistic image with this filename.'};


WrtFilt      = cfg_choice;
WrtFilt.tag  = 'WriteFiltImg';
WrtFilt.name = 'Write thresholded/filtered statistic image?';
WrtFilt.values  = {WF_no,WF_name};
WrtFilt.val  = {WF_no};
WrtFilt.help = {
    'Choose whether to write out the statistic image after filtering with the selected voxel-wise or cluster-wise thresholding.'
    '(The unfiltered statistic image and P-value images are always written out by default.)'
};

Report       = cfg_menu;
Report.tag   = 'Report';
Report.name  = 'Results displayed:';
Report.labels= {'Standard (MIP & table)','FWE report','FDR report'};
Report.values= {'MIPtable','FWEreport','FDRreport'};
Report.val   = {'MIPtable'};
Report.help  = {
    'Select type of results report.'
    '"Standard" shows usual maximum intensity project (MIP) and the table of P-values and x,y,z peak locations.'
    '"FWE report" shows the distribution of the maximum statistic value (and, if computed, the maximum cluster size), used to determine the FWE-corrected P-values.'
    '"FDR report" shows the histogram of voxel-wise uncorrected P-values, and the log-log plot of observed P-values versus (null hypothesis) expected P-values. The latter is the basis of the voxel-wise FDR-corrected P-values.'
    };

% ---- Export to NIDM, adapted from SPM's spm_cfg_esults
%--------------------------------------------------------------------------
% nsubj Number of subjects
%--------------------------------------------------------------------------
nsubj         = cfg_entry;
nsubj.tag     = 'nsubj';
nsubj.name    = 'Number of subjects';
nsubj.help    = {'Number of subjects.'};
nsubj.strtype = 'r';
nsubj.num     = [1 1];

%--------------------------------------------------------------------------
% label Label
%--------------------------------------------------------------------------
grplabel         = cfg_entry;
grplabel.tag     = 'label';
grplabel.name    = 'Label';
grplabel.help    = {'Group label.'};
grplabel.strtype = 's';
grplabel.num     = [0 Inf];

%--------------------------------------------------------------------------
% group 
%--------------------------------------------------------------------------
group      = cfg_branch;
group.tag  = 'group';
group.name = 'Group';
group.val  = {nsubj grplabel};
group.help = {['Number of subjects and labels per group. ', ...
    'For a single subject analysis, enter "1" and "single subject".']};

%--------------------------------------------------------------------------
% groups
%--------------------------------------------------------------------------
groups        = cfg_repeat;
groups.tag    = 'groups';
groups.name   = 'Groups';
groups.help   = {['Number of groups. ', ...
    'For a single subject analysis, specify one group.']};
groups.values = {group};
groups.num    = [1 Inf];

%--------------------------------------------------------------------------
% modality Modality
%--------------------------------------------------------------------------
modality        = cfg_menu;
modality.tag    = 'modality';
modality.name   = 'Modality';
modality.help   = {'Modality.'};
modality.labels = {'Anatomical MRI',...
                   'Functional MRI',...
                   'Diffusion MRI',...
                   'PET',...
                   'SPECT',...
                   'EEG',...
                   'MEG'
}';
modality.values = {'AMRI','FMRI','DMRI','PET','SPECT','EEG','MEG'};

%--------------------------------------------------------------------------
% refspace Reference space
%--------------------------------------------------------------------------
refspace        = cfg_menu;
refspace.tag    = 'refspace';
refspace.name   = 'Reference space';
refspace.help   = {['Reference space. For an experiment completed only ',...
    'within SPM, choose one of the first four options.']};
refspace.labels = {'Subject space (no normalisation)',...
                   'Normalised space (using segment)',...
                   'Normalised space (using old segment)',...
                   'Customised space',...
                   'Other normalised MNI space',...
                   'Other normalised Talairach space',...
}';
refspace.values = {'subject','ixi','icbm','custom','mni','talairach'};

%--------------------------------------------------------------------------
% export Export results to NIDM
%--------------------------------------------------------------------------
export_no      = cfg_const;
export_no.tag  = 'export_no';
export_no.name = 'No';
export_no.val  = {0};
export_no.help = {'Do not export to NIDM.'};


export_nidm    = cfg_branch;
export_nidm.tag  = 'nidm';
export_nidm.name = 'Export to NIDM';
export_nidm.val  = {modality refspace groups};
export_nidm.help = {'NIDM (Neuroimaging Data Model)'};

export        = cfg_choice;
export.tag    = 'export';
export.name   = 'Export results to NIDM';
export.help   = {['Export your results to NIDM and share them easily with your collaborators.']};
export.values = {export_nidm, export_no};
export.val  = {export_nidm};

% ----

snpmpp        = cfg_exbranch;
snpmpp.name   = 'Inference';
snpmpp.tag    = 'inference';
snpmpp.val    = {snpmres ThrType posneg WrtFilt Report export};
snpmpp.prog   = @snpm_run_pp;
%snpmpp.vout
snpmpp.help   = {'Examine the results of the SnPM computation.'};
		  

% snpmpp         = cfg_choice;
% snpmpp.name    = 'Results';
% snpmpp.tag     = 'Results';
% snpmpp.values  = {report infer};
% snpmpp.help    = {'Select the type of results you would like to examine.'};

function vout = snpm_bch_pp_vout(job)
% Determine what outputs will be present if this job is run. In this case,
% the structure of the inputs is fixed, and the output is always a single
% number. Note that input items may not be numbers, they can also be
% dependencies.

%%%%% LATER, add the output to this... i.e. the filter image
vout = cfg_dep;                        % The dependency object
vout.sname      = 'Add1: a + b';       % Displayed dependency name
vout.src_output = substruct('()',{1}); % The output subscript reference. This could be any reference into the output variable created during computation


% function t = Clus_check(job)
% t={};
% try, load(job.snpmres{1})
% catch
%   t={'Cannot open SnPM mat file!  Check that it is valid.'};
% end
% snpmcwd=fileparts(job.snpmres{1});
% load(fullfile(snpmcwd,'SnPMcfg.mat'))
% if isempty(t)
%   if ~bST
%     t={'No cluster size information saved; cannot perform cluster-wise inference.'};
%   end
% end
% return

% function t = CFth_check(job)
% t={};
% try, load(job.snpmres{1})
% catch
%   t={'Cannot open SnPM mat file!  Check that it is valid.'};
% end
% snpmcwd=fileparts(job.snpmres{1});
% load(fullfile(snpmcwd,'SnPMcfg.mat'))
% if isempty(t)
%   if bST & pU_St_Ut==-1 & isnan(job.CFth)
%     t={'Cluster inference performed in "slow" mode, so cluster-forming threshold must be set (cannot be NaN).'};
%   end
% end
% return


