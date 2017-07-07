function snpm_cp(CWD)
% Nonparametric Perm/Rando statistical analysis with General linear model
% FORMAT snpm_cp(CWD)
%
% CWD -	Directory containing SnPMcfg.mat configuration file
%_______________________________________________________________________
% 
% snpm_cp is the engine of the SnPM toolbox and implements the general
% linear model for a set of design matrices, each design matrix 
% constituting one permutation.  First the "correct" permutation
% is calculated in its entirety, then all subsequent permutations are 
% calculated, possibly on a plane-by-plane basis.  
%
% The output of snpm_cp parallels spm_spm: for the correct permutation
% image files containing parameter estimates, statistic values, and F
% values are saved (this is in distinction from SnPM2 and previous
% versions, where this information was saved in .mat files); the permutation 
% distribution of the statistic interest and (optionally) suprathreshold
% stats are also saved.  All results are written to the directory
% that CfgFile resides in.  IMPORTANT: Existing results are overwritten
% without prompting.
%
% Unlike spm_spm, voxels are not discarded on the basis of the F statistic.
% All gray matter voxels (as defined by the gray matter threshold) are
% retained for analysis; note that this will increase the size of all .mat
% files.
%
%
%-----------------------------------------------------------------------
%
% Output File Descriptions:
%
%   XYZ.mat contains a 3 x N matrix of the x,y and z location of the
% voxels in SPMF in mm (usually referring the the standard anatomical
% space (Talairach and Tournoux 1988)} (0,0,0) corresponds to the
% centre of the voxel specified by ORIGIN in the *.hdr of the original
% and related data.
%
%   BETA.mat contains a p x S matrix of the p parameter estimates at
% each of the S voxels for the correct permutation.  These parameters 
% include all effects specified by the design matrix.
%
%   SnPMt.mat contains a 1 x S matrix of the statistic of interest (either
% t or pseudo-t if variance smoothing is used) supplied for all S voxels at
% locations XYZ.
%
%   SnPMucp.mat contains a 1 x S matrix of the nonparametric P values of
% the statistic of interest supplied for all S voxels at locations XYZ.
%
%   SnPM.mat contains a collection of strings and matrices that pertain 
% to the analysis.  In contrast to spm_spm's SPM.mat, most of the essential
% matrices are in the any of the matrices stored here in the CfgFile
% and hence are not duplicated here.   Included are the number of voxels
% analyzed (S) and the image and voxel dimensions [V].  See below
% for complete listing.
%
% snpm_cp writes out the following image files (for each image, there are
% two files: .img and .hdr files)  
%  
% beta_**** (from 0001 to p): p images of p parameter estimates at each
% voxel for the correct permutation. These p parameters include all
% effects specified by the design matrix. 
%
% ResMS: One image of residual mean square errors at each voxel. 
% 
% (SnPM, like SPM, only implements single tailed tests. In the following
% files, '+' or '-' correspond to 'positive' or 'negative' effects (as in
% snpm_pp.m). Here, '+' images are the images for large values,
% indicating evidence against the null hypothesis in favour of a positive
% alternative (activation, or positive slope in a covariate analysis))
% 
% snpmT+ & snpmT-: Images of the statistic of interest (either t or
% pseduo-t if variance smoothing is used), positive or negative. 
% The numbers (i.e. not NaN) saved in snpmT+ images are also saved in the
% SnPMt.mat file. 
% 
% lP+ & lP-: Images of -log10(uncorrected non-parametric P-values,
% positive or negative).
% 
% lP_FWE+ & lP_FWE-: Images of -log10(FWE-corrected non-parametric
% P-values, positive or negative). Here, FWE-corrected non-parametric
% P-values are the proportion of the permutation distribution for the
% maximal statistic which exceeds the statistic image at the voxel. 
%
% lP_FDR+ & lP_FDR-: Images of -log10(FDR-corrected non-parametric
% P-values, positive or negative). 
%
% The following is an example of matlab codes for reading in an image file. 
% P='.../.../beta_0001.img';
% V=spm_vol(P);
% Y=spm_read_vols(V);
% Y(~isfinite(Y))=[]; %delete NaN values from vector Y.
%
%
%-----------------------------------------------------------------------
%
% As an "engine", snpm_cp does not produce any graphics; if the SPM windows
% are open, a progress thermometer bar will be displayed.  
%
% If out-of-memory problems are encountered, the first line of defense is to
% run snpm_cp in a virgin matlab session with out first starting SPM.
%
%
% Variables saved in SnPM.mat
%=======================================================================
%
% S              Volume analyzed (in voxels)
% V              Volume handles (see spm_vol)
% df             Residual degrees of freedom of raw t-statistic
% MaxT           nPerm x 2 matrix of [max;min] t-statistics per perm
% ST_Ut          Threshold above which suprathreshold info was collected.
%                Voxel locations, t and perm are saved in SnPM_ST.mat for
%                t's greater than ST_Ut. ST_Ut=Inf if not saving STCdata
%
% s_SnPM_save    List of variables saved in SnPM.mat file
% CfgFile        SnPM config sile used (full pathname)
% s_SnPMcfg_save List of variables saved in SnPMcfg.mat file
% 
% Data structure of SnPM_ST.mat: suprathreshold stats (if collected)
%-----------------------------------------------------------------------
% 5xn matrix, each column containing:
%       [x, y, z, abs(T), perm]'
%       perm is negative if T was negative
%
%_______________________________________________________________________
% Copyright (C) 2013-2014 The University of Warwick
% Id: snpm_cp.m  SnPM13.01 2014/01/31
% Thomas Nichols, Andrew Holmes



%-----------------------------functions-called------------------------
% spm_append_96
% spm_conv
% spm_figure
% spm_select
% spm_hread
% spm_invTcdf
% spm_matrix
% spm_progress_bar
% spm_smooth
% spm_str_manip
%-----------------------------functions-called------------------------


% Programmers / Hackers help / Code notes...
%=======================================================================
% snpm_cp is modeled after SPM95's spm_spm, which is different from the
% current spm_spm (which is similar to SPM99).  While the current
% (SPM99-SPM5) spm_spm reads data in 'planks' (a partial or whole slice),
% snpm_cp reads either a plane at a time or reads the entire dataset into
% memory ('Volumetric mode', bVolm=1).
%
% If bVolm is true, this function will load the entire dataset (all
% planes, all subjects) into memory.  If bVolm is false the dataset
% will be loaded a plane at a time (all subjects), but smoothing the
% variance in the z-direction is not permitted.  A version that 
% does not load the whole volume but allows volumetric smoothing
% is under development.
% 
% For designs with large number of permutations AND bVolm true, a possible
% approach would be have a stopping feature where the user had decided
% "enough" permutations had run.  Not sure if this is useful.

%-Variable "decoder"
%-----------------------------------------------------------------------
% bWin    - Do we have windows?
% bVarSm  - Variance Smoothing?
% bVolm   - Work on whole volume at once?
% q       - number of observations
% p       - number of predictors
% r       - Model degrees of freedom
% df      - Residual degrees of freedom
% nPerm   - number of permutations
% WorkDim - Number of voxels read in (either a plane's worth or the whole lot)
% MaxT    - Permutation Distribution of intensity maximum
% nP      - WorkDim vector of nonparametric P-values
%
%-Supratreshold Threshold
% 
% STalpha - if parametric T, critical val for STalpha ised
% STprop  - if pseudo T, 100*(1-STprop)%ile of correct perm's values used


%-Setup
%=======================================================================
global defaults
if isempty(defaults), spm_defaults; end


fprintf('\nSnPM: snpm_cp\n'),fprintf('%c','='*ones(1,72)),fprintf('\n')
disp('Initialising...');
bWin = ~isempty(spm_figure('FindWin','Interactive'));
s_SnPM_save = 's_SnPM_save CfgFile s_SnPMcfg_save S V df1 df MaxT ST_Ut STAT';

%-Check arguments & parameters from CfgFile
%-----------------------------------------------------------------------
if nargin == 0
  tmp = spm_select(1,'SnPMcfg.mat','Select SnPMcfg.mat CfgFile...');
  drawnow
  CWD = spm_str_manip(tmp,'hd');
end
if strcmp(CWD, '.')
  CWD=pwd;
end  
if ~strcmp(pwd,CWD)
  cd(CWD)
  CWD=pwd;
  fprintf('Changing directory to %s\n',CWD);
end
CfgFile = fullfile(CWD,'SnPMcfg.mat');

%-Load config file & catch all problem cases now
%-----------------------------------------------------------------------
load(CfgFile);

if strcmp(sDesFile, 'snpm_pi_OneSampT') || ...
        strcmp(sDesFile, 'snpm_pi_ANOVAwithinS')
    % Sign flipping
    nidm.ErrorModel_hasErrorDistribution = {'obo_nonparametricdistribution', 'obo_symmetricdistribution'};
    nidm.ErrorModel_errorVarianceHomogeneous = false;
    nidm.ErrorModel_varianceMapWiseDependence = 'nidm_IndependentParameter';
    nidm.ErrorModel_hasErrorDependence = 'nidm_IndependentError';
else
    % Permutation
    nidm.ErrorModel_hasErrorDistribution = 'obo_nonparametricdistribution';
    nidm.ErrorModel_errorVarianceHomogeneous = true;
    nidm.ErrorModel_varianceMapWiseDependence = 'nidm_IndependentParameter';
    % TODO: the 'obo_exchangeable' term is not yet in STATO
    nidm.ErrorModel_hasErrorDependence = 'obo_exchangeable';
    nidm.ErrorModel_dependenceMapWiseDependence = 'nidm_IndependentParameter';
end

% TODO: check this is correct
nidm.ModelParameterEstimation_withEstimationMethod = 'obo_ordinaryleastsquaresestimation';

if isempty([H C])
  error('SnPM:NoModel', 'No model specified; [H C] empty'); 
end
if ~isempty(H) & ~isempty(C)
    error('SnPM:HierarchicalAndCov', 'Cannot have both heirachical and covariate effects'); 
end
if size(CONT,2) ~= size([H C B G],2)
    error('SnPM:InvalidContrast','Contrast problem; wrong number of columns'); 
end

if size(CONT,1) > 1 % F-contrast
  nidm.Contrasts(1).contrastweightmatrix_value = CONT;
  warning('SnPM:FContrast', ...
          'F contrast!  F statistic images are being created.'); 
  STAT = 'F';
  nidm.Contrasts(1).StatisticMap_statisticType = 'obo_Fstatistic';
  
  if (CONT(1,:) == -CONT(2,:))
    CONT = CONT(1,:);
  end
  con_name = 'Positive';
  nidm.Contrasts(1).StatisticMap_contrastName = con_name;  
  con_neg_name = '';
else
  con_name = 'Positive';  
  nidm.Contrasts(1).StatisticMap_contrastName = con_name;  

  con_neg_name = 'Negative';
  nidm.Contrasts(2).StatisticMap_contrastName = con_neg_name;  

  STAT = 'T';
  if bVarSm
    nidm.Contrasts(1).StatisticMap_statisticType = 'nidm_smoothedtstatistic';
    nidm.Contrasts(2).StatisticMap_statisticType = 'nidm_smoothedtstatistic';
  else
    nidm.Contrasts(1).StatisticMap_statisticType = 'obo_tstatistic';
    nidm.Contrasts(2).StatisticMap_statisticType = 'obo_tstatistic';
  end
  nidm.Contrasts(1).ContrastMap_contrastName = ['Positive T-Contrast: [' mat2str(CONT) ']'];
  nidm.Contrasts(1).contrastweightmatrix_value = CONT;
  nidm.Contrasts(2).ContrastMap_contrastName = ['Negative T-Contrast: [' mat2str(-CONT) ']'];
  nidm.Contrasts(2).contrastweightmatrix_value = -CONT;
end

if rank(CONT)<size(CONT,1)
  [u s] = spm_svd(CONT'); % Kill zero-rank components
  CONT = full(u*sqrt(s))';
end
if ~bVolm & bVarSm & vFWHM(3)
  error('SnPM:ZSmoothVolume', 'Cannot z-smooth variance in non-volumetric mode'); 
end
if exist('bVarAlph')~=1
  bVarAlph=0; 
end
if bVarAlph & ~(~bVarSm & bVolm)
  error('SnPM:AlphaVolumePseudo', 'No pseudo t or nonvolumetric w/ variable alpha');
end
if ~bVolm & pU_ST_Ut>=0
  error('SnPM:STCSNotVolume', 'Must work volumetrically to computer STCS on-the-fly');
end
% Re-map files to avoid Endian headaches; note if NaN's available
NaNrep=0;
for i = 1:length(V)
    curr_pinfo = V(i).pinfo;% Added to keep scaling
    V(i) = spm_vol([V(i).fname ',' num2str(V(i).n)]);
    original_pinfo = V(i).pinfo;
    V(i).pinfo = curr_pinfo;
    NaNrep = NaNrep | spm_type(V(i).dt(1),'nanrep');
end

%-Delete files from previous analyses, if they exist
%-----------------------------------------------------------------------
files = {	'^ResMS\..{3}$','^beta_.{4}\..{3}','^GrandMean\..{3}', ...
        '^mask\..{3}', '^conse\..{3}',...
        '^con.{2}\..{3}','^con.{1}\..{3}',...
        '^lP_.{4}\..{3}',...
		'^lP.{1}\..{3}','^snpm.{2}\..{3}','^snpm.{1}\..{3}'};

for i=1:length(files)
  j = spm_select('List',pwd,files{i});
  for k=1:size(j,1)
    spm_unlink(deblank(j(k,:)));
  end
end

spm_unlink SnPM.mat SnPM_ST.mat SnPMt.mat SnPMucp.mat XYZ.mat SnPM_pp.mat ...
    SnPM_pp_Neg.mat STCS.mat


%-Parameters & Initialisation
%=======================================================================

%-Suprathreshold parameters
%-----------------------------------------------------------------------
STalpha = snpm_get_defaults('STalpha'); 
STprop  = snpm_get_defaults('STprop');

s_SnPM_save = [s_SnPM_save ' STalpha STprop'];  % Save for PP

%-Work out degrees of freedom
%-----------------------------------------------------------------------
q       = size([H C B G],1);		%-# observations
p       = size([H C B G],2);		%-# predictors
r       = rank([H C B G]);		%-Model degrees of freedom
df      = q - r;			%-Residual degrees of freedom
nPerm   = size(PiCond,1);		%-# permutations

nidm.NonParametricNullDistribution_numberOfPermutations = nPerm;
nidm.NonParametricNullDistribution_hasResamplingScheme = 'nidm_Permutation';
nidm.NonParametricNullDistribution_hasApproximationMethod = 'nidm_MonteCarlo';
nidm.NonParametricNullDistribution_maximumNumberOfPermutations = nPerm_max;

for i = 1:numel(nidm.Contrasts)
    nidm.Contrasts(i).StatisticMap_errorDegreesOfFreedom = df;
end

%-Get ORIGIN, etc
%-----------------------------------------------------------------------
DIM    = [V(1).dim(1)   V(1).dim(2)   V(1).dim(3)];
M=V(1).mat(1:3, 1:3);
VOX=sqrt(diag(M'*M))';
MAT    = V(1).mat;
IMAT   = inv(MAT);
ORIGIN = IMAT(1:3,4);


%-Var-alpha stuff
%-----------------------------------------------------------------------
bMask = 0;
if bVarAlph
  Vwt    = spm_vol(Pwt);
  MinwP  = repmat(Inf,nPerm,2);
  s_SnPM_save = [s_SnPM_save ' MinwP Pwt Vwt'];
  bMask = 1;
elseif exist('Pwt')==1
  Vwt    = spm_vol(Pwt);
  s_SnPM_save = [s_SnPM_save ' Pwt Vwt'];
  bMask = 1;
elseif ~isempty(MASK)
  Vwt    = spm_vol(MASK);
  s_SnPM_save = [s_SnPM_save ' Vwt'];
  bMask = 1;
end
% Add updated nidm to the saved variables
s_SnPM_save = [s_SnPM_save ' con_name con_neg_name nidm '];

%-Useful quantities - handy for later
%-----------------------------------------------------------------------
xdim     = DIM(1);			%-X dimension
ydim     = DIM(2);			%-Y dimension
zdim     = DIM(3);			%-Z dimension
PlDim    = xdim*ydim;			%-Plane size in voxels
VolDim   = xdim*ydim*zdim;		%-Volume size in voxels
if bVolm,
  WorkDim = VolDim;	%-Working dimension (if volumetric)
else
  WorkDim = PlDim;	%-Working dimension (if plane by plane)
end

%-Location vectors --> In units of mm <--
%-----------------------------------------------------------------------
[y x] = meshgrid([1:ydim],[1:xdim]');
x     = x(:)';
y     = y(:)';
z     = (1:zdim);

xyPl = [x;y];  % All x & y's in one plane

%-Initialize variables
%-----------------------------------------------------------------------
TH    = TH*ones(1,WorkDim);	%-Global activities
S     = 0;			%-Volume analyzed
MaxT  = repmat(-Inf,nPerm,2);	%-Max t
nP    = zeros(1,WorkDim);	%-Nonparam P's
XYZ_total=[];                   %-the variable for keeping all XYZ
%-If working plane by plane, preallocate Q & XYZ for speed/mem. efficiency
if ~bVolm, 
  Q    = zeros(1,PlDim);
  XYZ  = zeros(3,PlDim);
end

SmTime = 0;					%-Smoothing time
perm   = 0;

% Initialize structure template
%----------------------------------------
Vt=V(1);

%
%-Initialize image structures.
%
for ii=1:p
  fname{ii}= sprintf('beta_%04d.img',ii);
  descrip{ii}=sprintf('beta_%04d hats',ii);
  Vbeta(ii)=snpm_clone_vol(Vt,fname{ii},descrip{ii}); 
end  
nidm.DesignMatrix_regressorNames = descrip;
nidm.ParameterEstimateMaps = fname;

Vbeta = spm_create_vol(Vbeta);

VResMS=snpm_clone_vol(Vt,'ResMS.img','Residual sum-of-squares');
VResMS=spm_create_vol(VResMS);
nidm.ResidualMeanSquaresMap_atLocation = 'ResMS.img';

if bVarSm==0
  str = sprintf('%c_{%d} statistic',STAT,df);
else
  if STAT=='T'
    str = sprintf('SmVar T_{%d} statistic, %fx%fx%f VarSm',df,vFWHM);
  elseif STAT=='F'
    str = sprintf('SmVar F_{%d,%d} statistic, %fx%fx%f VarSm',...
		  [rank(CONT) df],vFWHM);
  end
end

Vmask=snpm_clone_vol(Vt,'mask.img',str);
Vmask=spm_create_vol(Vmask);
nidm.MaskMap_atLocation = 'mask.img';

Vgmean=snpm_clone_vol(Vt,'GrandMean.img','GrandMean');
Vgmean=spm_create_vol(Vgmean);
nidm.GrandMeanMap_atLocation = 'GrandMean.img';

if STAT=='T'
  VT_pos=snpm_clone_vol(Vt,'snpmT+.img',[str,' (+ve)']);
  VT_pos=spm_create_vol(VT_pos);
  
  VCON_pos=snpm_clone_vol(Vt,'con+.img',[str,' (+ve)']);
  VCON_pos=spm_create_vol(VCON_pos);
  
  VCONSE=snpm_clone_vol(Vt,'conse.img', str);
  VCONSE=spm_create_vol(VCONSE);
  
  nidm.Contrasts(1).ContrastStandardErrorMap_atLocation = 'conse.img';
  nidm.Contrasts(1).ContrastMap_atLocation = 'con+.img';
  nidm.Contrasts(1).StatisticMap_atLocation = 'snpmT+.img';
  
  VT_neg=snpm_clone_vol(Vt,'snpmT-.img',[str,' (-ve)']);
  VT_neg=spm_create_vol(VT_neg);
  
  VCON_neg=snpm_clone_vol(Vt,'con-.img',[str,' (+ve)']);
  VCON_neg=spm_create_vol(VCON_neg);
  
  % TODO: check it's fine to use same conse image for pos and neg contrasts  
  nidm.Contrasts(2).ContrastStandardErrorMap_atLocation = 'conse.img';
  nidm.Contrasts(2).ContrastMap_atLocation = 'con-.img';
  nidm.Contrasts(2).StatisticMap_atLocation = 'snpmT-.img';
elseif STAT=='F'
  VF=snpm_clone_vol(Vt,'snpmF.img',str);
  VF=spm_create_vol(VF);
  
  nidm.Contrasts(1).StatisticMap_atLocation = 'snpmF.img';

  VFnum=snpm_clone_vol(Vt,'snpmFnum.img',str);
  VFnum=spm_create_vol(VFnum);
  
  nidm.Contrasts(1).ContrastExplainedMeanSquareMap_atLocation = 'snpmFnum.img';
end

VlP_pos=snpm_clone_vol(Vt, 'lP+.img', '-log10(uncor. non-para. P, +ve)');
VlP_pos=spm_create_vol(VlP_pos);
VlP_FWE_pos=snpm_clone_vol(Vt, 'lP_FWE+.img','-log10(FWE-corr. P, +ve)');
VlP_FWE_pos=spm_create_vol(VlP_FWE_pos);
VlP_FDR_pos=snpm_clone_vol(Vt, 'lP_FDR+.img','-log10(FDR-corr. P, +ve)');
VlP_FDR_pos=spm_create_vol(VlP_FDR_pos);

if STAT=='T'
  VlP_neg=snpm_clone_vol(Vt, 'lP-.img', '-log10(uncor. non-para. P, -ve)');
  VlP_neg=spm_create_vol(VlP_neg);
  VlP_FWE_neg=snpm_clone_vol(Vt, 'lP_FWE-.img','-log10(FWE-corr. P, -ve)');
  VlP_FWE_neg=spm_create_vol(VlP_FWE_neg);
  VlP_FDR_neg=snpm_clone_vol(Vt, 'lP_FDR-.img','-log10(FDR-corr. P, -ve)');
  VlP_FDR_neg=spm_create_vol(VlP_FDR_neg);
end

if bVarAlph
  VlwP=snpm_clone_vol(Vt, 'lwP.img','-log10(weighted p-value)');
  VlwP=spm_create_vol(VlwP);
end  

%	
%-Initialize image data.
%
lP_pos_image=repmat(NaN,1,VolDim);
lP_FWE_pos_image=repmat(NaN,1, VolDim);
lP_FDR_pos_image=repmat(NaN,1, VolDim);
if STAT=='T'
  lP_neg_image=repmat(NaN,1,VolDim);
  lP_FWE_neg_image=repmat(NaN,1, VolDim);	
  lP_FDR_neg_image=repmat(NaN,1, VolDim);	
end

%=======================================================================
% - C O R R E C T   P E R M U T A T I O N
%=======================================================================
% Work out correct permuation completely. Separating the first
% permutation simplifies the permutation loop (fewer conditionals) and
% allows determination of pseudo-t threshold when saving supratheshold
% statistics.
disp('Working on correct permutation...');

SnPMt=[]; %Initialzie SnPMt,which will store the t's from correct permutation.

for i = 1:zdim
  
  %-Initialize the image data for this slice/volume
  %---------------------------------------------------------------------
  gmean_image=repmat(NaN,1,WorkDim);
  BETA_image=repmat(NaN,p,WorkDim);
  ResSS_image=repmat(NaN,1,WorkDim);
  mask_image=repmat(NaN,1,WorkDim);
  CONSE_image=repmat(NaN,1,WorkDim);
  if STAT=='T'
    T_pos_image=repmat(NaN,1,WorkDim);
    T_neg_image=repmat(NaN,1,WorkDim);
    CON_pos_image=repmat(NaN,1,WorkDim);
    CON_neg_image=repmat(NaN,1,WorkDim);
  elseif STAT=='F'
    F_image=repmat(NaN,1,WorkDim);
  end
  
  if bVarAlph
    lwP_image=repmat(NaN, 1, WorkDim); 
  end
  
  %-Form data matrix for this slice/volume
  %---------------------------------------------------------------------
  X     = zeros(q,WorkDim);
  if bMask, 
    Wt = zeros(1,WorkDim); 
  else
    Wt = 1;
  end
  if bVolm
    for j = 1:q
      for k = 1:zdim
        tmp    = spm_slice_vol(V(j),spm_matrix([0 0 k]), ...
			       [xdim ydim],0);
        X(j,(k-1)*PlDim+1:k*PlDim) = tmp(:)';
      end
    end
  else
    for j = 1:q
      tmp    = spm_slice_vol(V(j),spm_matrix([0 0 i]), ...
			     [xdim ydim],0);
      X(j,:) = tmp(:)';
    end
  end
  if bMask
    if bVolm
      for k = 1:zdim
        j = Vwt.mat\MAT*[xyPl;repmat(k,1,PlDim);ones(1,PlDim)];
        tmp    = spm_get_data(Vwt,j,false);
        tmp(~isfinite(tmp) | tmp<0) = 0;
        Wt(1,(k-1)*PlDim+1:k*PlDim) = tmp(:)';
      end
    else
      j = Vwt.mat\MAT*[xyPl;repmat(i,1,PlDim);ones(1,PlDim)];
      tmp = spm_get_data(Vwt,j,false);
      tmp(~isfinite(tmp) | tmp<0) = 0;
      Wt  = tmp(:)';
    end
  end
  
  
  %-Eliminate background voxels (based on threshold TH), and
  % eliminate voxels where there are no differences across scans.
  %---------------------------------------------------------------------
  if ImMASK & NaNrep==0
    Q = find(all(X>TH) & any(diff(X)) & Wt & all(X~=0));
  else
    Q = find(all(X>TH) & any(diff(X)) & Wt);
  end
  

  if length(Q)
   
    X     = X(:,Q);
    S     = S + length(Q); 			%-Volume 
    if bVolm
      XYZ   = [ x(rem(Q-1,PlDim)+1);          ...
		y(rem(Q-1,PlDim)+1);            ...
		z(ceil(Q/PlDim))      ]; %-Locations
    else
      XYZ   = [ x(rem(Q-1,PlDim)+1);         ...
		y(rem(Q-1,PlDim)+1);         ...
		z(i)*ones(1,length(Q))];	%-Locations
    end 
    
    % Convert Voxels to mm's
    XYZ = MAT*[XYZ;ones(1,length(Q))]; XYZ(4,:) = [];
    
    if (bMask)
      Wt    = Wt(1,Q);
    end
    
    perm = 1;

    %-Estimate parameters and sum of squares due to error.
    % Use pseudo inverse rather than BETA=inv(D'*D)*D'*X for 
    % D = DesMtx, to allow for non-unique designs. See matlab help.
    %-----------------------------------------------------------------
    gmean = mean(X);
    BETA  = pinv([H C B G])*X;
    ResSS = sum((X - [H C B G]*BETA).^2);
    
    %-Variance smoothing.
    % Blurred mask is used to truncate kernal to brain; if not
    % used variance at edges would be underestimated due to
    % convolution with zero activity out side the brain.
    %-----------------------------------------------------------------
    if bVarSm
      if bVolm
        SmResSS   = zeros(xdim, ydim, zdim);
        SmMask    = zeros(xdim, ydim, zdim);
        TmpVol    = zeros(xdim, ydim, zdim);
        TmpVol(Q) = ones(size(Q));
        % FWHM in voxels (and not in mm) as TmpVol is not a struct 
        spm_smooth(TmpVol,SmMask,vFWHM./VOX);
        TmpVol(Q) = ResSS;
        % FWHM in voxels (and not in mm) as TmpVol is not a struct 
        spm_smooth(TmpVol,SmResSS,vFWHM./VOX);

        ResSS     = SmResSS(Q)./SmMask(Q);
      else
        TmpPl     = zeros(xdim,ydim);
        TmpPl(Q)  = ones(size(Q));
        SmMask    = spm_conv(TmpPl, vFWHM(1)/VOX(1),vFWHM(2)/VOX(2));
        TmpPl(Q)  = ResSS;
        SmResSS   = spm_conv(TmpPl, vFWHM(1)/VOX(1),vFWHM(2)/VOX(2));

        ResSS     = SmResSS(Q)./SmMask(Q);
      end
    end
    
    %-Compute t-statistics for specified compounds of parameters
    %-----------------------------------------------------------
    T      = zeros(1,size(BETA,2));
    Fnum   = zeros(1,size(BETA,2));
    CON    = zeros(1,size(BETA,2));
    Co     = CONT;
    if STAT=='T'
      CON(1,:) = Co*BETA;
      CONSE(1,:) = sqrt((ResSS*(Co*pinv([H C B G]'*[H C B G])*Co'))/df);
      % t, as usual
      T(1,:) = Co*BETA./sqrt((ResSS*(Co*pinv([H C B G]'*[H C B G])*Co'))/df);
    else
      % F!
      pX   = pinv([H C B G]);
      T(1,:) = (sum(((Co*BETA)'*inv(Co*pinv([H C B G]'*[H C B G])*Co'))' .* ...
		    (Co*BETA),1)/rank(Co)) ./ (ResSS/df);
      Fnum(1,:) = (sum(((Co*BETA)'*inv(Co*pinv([H C B G]'*[H C B G])*Co'))' .* ...
		    (Co*BETA),1)/rank(Co));
    end	
    
    %-Save Max T statistic
    %-----------------------------------------------------------
    MaxT(perm,:) = max([ max(T(1,:)), -min(T(1,:));   ...
		    MaxT(perm,1),  MaxT(perm,2) ]);
    
    %-Save min weighted p-value
    %-----------------------------------------------------------
    if bVarAlph,
      MinwP(perm,:) = min([ min(Wt.*(1-spm_Tcdf(T(1,:),df))),     ...
		    min(Wt.*(1-spm_Tcdf(-T(1,:),df)));    ...
		    MinwP(perm,1), MinwP(perm,2) ]);
    end
    
    %-Save weighted p-value (later converted into corr'd wt'd p-val)
    %-----------------------------------------------------------
    if bVarAlph,
      wP = [Wt.*(1-spm_Tcdf( T(1,:),df));  ...
	    Wt.*(1-spm_Tcdf(-T(1,:),df))];
    end

    %-Adjustment (remove effects of no interest) & save
    %-----------------------------------------------------------
    XA = X - [zeros(size([H C])) B G]*BETA;
    
    %
    %- New! Write out data images.
    %- Input image data.
    gmean_image(:,Q)=gmean;
    BETA_image(:,Q)=BETA;
    ResSS_image(:,Q)=ResSS;
    mask_image(:,Q)=1;
    if STAT=='T'
      T_pos_image(:,Q)=T;
      T_neg_image(:,Q)=-T;
      
      CON_pos_image(:,Q)=CON;
      CON_neg_image(:,Q)=-CON;
      
      CONSE_image(:,Q)=CONSE;
    elseif STAT=='F'
      F_image(:,Q)=T;
      Fnum_image(:,Q)=Fnum;
    end
	
    if bVarAlph
      lwP_image(:,Q)=-log10(wP);
    end  
	
    if bVolm
      SnPMt=T; % save T's
    else
      SnPMt=[SnPMt,T]; % save T's.
    end

    XYZ_total=[XYZ_total, XYZ];
	
    %if bVarAlph,
    %  spm_append_96('SnPMwP',wP);      % wt'd p-val of Stat du jour
    %end
	    
     
  end %(if length(Q))
    
  %-The image of the volume or the slice should be written out no matter length(Q)=1
  %or 0. 
  if bVolm
    for ii=1:p
      BETA_vol=reshape(BETA_image(ii,:),DIM(1),DIM(2),DIM(3));
      spm_write_vol(Vbeta(ii),BETA_vol);
    end
    
    ResSS_vol=reshape(ResSS_image,DIM(1),DIM(2),DIM(3));
    spm_write_vol(VResMS, ResSS_vol);
    
    % Analysis mask
    mask_vol=reshape(mask_image,DIM(1),DIM(2),DIM(3));
    spm_write_vol(Vmask,mask_vol);

    % Grand mean
    gmean_vol = reshape(gmean_image,DIM(1),DIM(2),DIM(3));
    spm_write_vol(Vgmean,gmean_vol);
    
    if STAT=='T'
      T_pos_vol=reshape(T_pos_image,DIM(1),DIM(2),DIM(3));
      spm_write_vol(VT_pos,T_pos_vol);
      
      T_neg_vol=reshape(T_neg_image,DIM(1),DIM(2),DIM(3));
      spm_write_vol(VT_neg,T_neg_vol);

      CON_pos_vol=reshape(CON_pos_image,DIM(1),DIM(2),DIM(3));
      spm_write_vol(VCON_pos,CON_pos_vol);
      
      CON_neg_vol=reshape(CON_neg_image,DIM(1),DIM(2),DIM(3));
      spm_write_vol(VCON_neg,CON_neg_vol);
      
      CONSE_vol=reshape(CONSE_image,DIM(1),DIM(2),DIM(3));
      spm_write_vol(VCONSE,CONSE_vol);
      
    elseif STAT=='F'
      F_vol=reshape(F_image,DIM(1),DIM(2),DIM(3));
      spm_write_vol(VF,F_vol);
      
      Fnum_vol=reshape(Fnum_image,DIM(1),DIM(2),DIM(3));
      spm_write_vol(VFnum,Fnum_vol);
    end
	  
    if bVarAlph
      lwP_vol=reshape(lwP_image, DIM(1), DIM(2), DIM(3));
      spm_write_vol(VlwP, lwP_vol);
    end
	  
  else
    % Analysis mask
    mask_plate=reshape(mask_image,DIM(1),DIM(2));
    spm_write_plane(Vmask,mask_plate,i);  
      
    % Grand mean
    gmean_plate=reshape(gmean_image, DIM(1), DIM(2));
    spm_write_plane(Vgmean,gmean_plate,i);

    for ii=1:p
      BETA_plate=reshape(BETA_image(ii,:), DIM(1), DIM(2));
      spm_write_plane(Vbeta(ii),BETA_plate,i);
    end
    
    ResSS_plate=reshape(ResSS_image, DIM(1), DIM(2));
    spm_write_plane(VResMS,ResSS_plate,i);
	    
    if STAT=='T'
      T_pos_plate=reshape(T_pos_image, DIM(1), DIM(2));
      spm_write_plane(VT_pos,T_pos_plate,i);
      
      T_neg_plate=reshape(T_neg_image, DIM(1), DIM(2));
      spm_write_plane(VT_neg,T_neg_plate,i);
      
      CON_pos_plate=reshape(CON_pos_image,DIM(1),DIM(2));
      spm_write_plane(VCON_pos,CON_pos_plate,i);
      
      CON_neg_plate=reshape(CON_neg_image,DIM(1),DIM(2));
      spm_write_plane(VCON_neg,CON_neg_plate,i);
      
      CONSE_plate=reshape(CONSE_image,DIM(1),DIM(2));
      spm_write_plane(VCONSE,CONSE_plate,i);
	  
    elseif STAT=='F'
      F_plate=reshape(F_image, DIM(1), DIM(2));
      spm_write_plane(VF,F_plate,i);
      
      Fnum_plate=reshape(Fnum_image, DIM(1), DIM(2));
      spm_write_plane(VFnum,Fnum_plate,i);
    end  
	    
    if bVarAlph
      lwP_plate=reshape(lwP_image, DIM(1), DIM(2));
      spm_write_plane(VlwP, lwP_plate,i);
    end
  end	  
    
  %-Whole volume complete in one pass if volumetric
  if bVolm
    break
  end
end

% Make an error if actually 'no voxels in brain'. 
if perm==0, error('SnPM:NoVoxelsInBrain', 'No voxels in brain'); end

save SnPMt SnPMt

%save XYZ in a XYZ.mat file.
%===============
XYZ=XYZ_total;
save XYZ XYZ


%-Set SupraThreshold t-threshold
%=======================================================================
if bST 
  if pU_ST_Ut==-1 % No threshold has been set yet.
    if bVarSm
      load SnPMt
      SnPMt = sort(SnPMt')';
      ST_Ut = SnPMt(round((1-STprop)*length(SnPMt)));
      % clear SnPMt
    else
      ST_Ut = spm_invTcdf(1 - STalpha, df);
    end
  else            % A threshold has been set.    
    if bVarSm
        if (pU_ST_Ut < 1)
            ST_Ut = spm_invNcdf(1-pU_ST_Ut);
            warning('snpm_cp:pseudoTFormingThresholdP',...
                ['Pseudo-T cluster-forming threshold defined by '...
                'P-value using Gaussian approximation P=' num2str(pU_ST_Ut)...
                ' -> Z=' num2str(ST_Ut) '; actual Pseudo-T threshold '...
                'unknown but may be higher than ' num2str(ST_Ut) '.']);
        else
            ST_Ut=pU_ST_Ut;
        end
    else
      if (pU_ST_Ut>1)
        ST_Ut=pU_ST_Ut;
      else
        if STAT == 'T'
                ST_Ut = spm_invTcdf(1-pU_ST_Ut, df);
        else
                ST_Ut = spm_invFcdf(1-pU_ST_Ut, df1, df);
        end
      end   
    end 
  end  
  StartPerm = 1;   % redo 1st perm for ST stats
else
  ST_Ut = Inf;
  StartPerm = 2;
end

%-Save correctly labeled T's
if bVolm & (StartPerm==2)
  T0 = T;
  nPtmp = ones(size(T));
  if bhPerms
    nPtmp = nPtmp + (T0<=0);
  end
else
  StartPerm = 1;
  nPtmp=[];
end


%=======================================================================
% - C O M P U T E   F O R   P E R M U T A T I O N S
%=======================================================================
%-Cycle over planes (or just once for volumetric mode)

%-If working plane by plane, preallocate Q & XYZ for speed/mem. efficiency
if ~bVolm, 
  Q    = zeros(1,PlDim);
  XYZ  = zeros(3,PlDim);
end

%-Setup progress bar
if bWin && ~bVolm
  spm_progress_bar('Init',zdim,'Looping over (perms within) planes...','Plane')
elseif bWin
  spm_progress_bar('Init',nPerm,'Volumetric mode...','Permutation')
  spm_progress_bar('Set',StartPerm-1)
end
tic %-Start the clock: Timing code is commented with "clock" symbol: (>)

%-Loop over planes (breaks out after first loop if bVolm)
%-----------------------------------------------------------------------
nP = [];
for i = 1:zdim
    
  PlStart=toc;SmTime=0; %-Timestamp (>) 
	
  if bVolm
    disp('Working on the whole volume');
  else
    fprintf('\tPlane %3d: ',i); 
  end
    
  %-Form data matrix for this slice (done in correctPerm code above if bVolm)
  %---------------------------------------------------------------------
  if ~bVolm
    X     = zeros(q,PlDim);
    for j = 1:q
      tmp = spm_slice_vol(V(j),spm_matrix([0 0 i]),[xdim ydim],0);
      X(j,:) = tmp(:)';
    end
    if bMask
      j = Vwt.mat\MAT*[xyPl;repmat(i,1,PlDim);ones(1,PlDim)];
      tmp = spm_get_data(Vwt,j,false);
      tmp(~isfinite(tmp) | tmp<0) = 0;
      Wt  = tmp(:)';
    end
    
    %-Eliminate background voxels (based on global threshold TH),
    % and eliminate voxels where there are no differences across scans.
    %-----------------------------------------------------------------
    Q = find(all(X > TH) & any(diff(X)) & Wt);
  end % (if ~bVolm)
    
  if length(Q)
    if ~bVolm, 
      X = X(:,Q); 
    end	%-Already done if bVolm
    
    if bST & ~bVolm			%-XYZ already done if bVolm
      XYZ   = [ x(rem(Q-1,PlDim)+1);          ...
        y(rem(Q-1,PlDim)+1);          ...
        z(i)*ones(length(Q),1)];	%-Locations
      XYZ = MAT*[XYZ;ones(1,length(Q))]; 
      XYZ(4,:) = [];
    end 
		
    if bVarSm & ~bVolm			%-Smoothing & plane-by-plane
      SmStart = toc;			%-Timestamp (>)
      TmpPl     = zeros(xdim,ydim);
      TmpPl(Q)  = ones(size(Q));
      SmMask     = spm_conv(TmpPl, vFWHM(1)/VOX(1),vFWHM(2)/VOX(2));
      SmTime = SmTime + toc-SmStart;	%-Timestamp (>)
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initialize structure STCS 
    if bST & pU_ST_Ut>=0
      STCS = snpm_STcalc('init',nPerm); 
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	
    %-Loop over permutations
    %-----------------------------------------------------------------
    for perm = StartPerm:nPerm
      PmStart = toc;			%-Timestamp (>)

      if bVolm 
        SmTime=0;			%-Timestamp (>)
      else
        clear T BETA ResSS; 	%-Clean up
      end

      %-Rebuild H C for current permuation
      %-----------------------------------------------------------
      HC = eval(sHCform);
      
      %-Estimate parameters and sum of squares due to error
      %-Use pseudo inverse rather than BETA=inv(D'*D)*D'*X
      % for D = DesMtx, to allow for non-unique designs. 
      % See matlab help.
      %-----------------------------------------------------------
      BETA  = pinv([HC B G])*X;
      ResSS = sum((X - [HC B G]*BETA).^2);
      
      if bVarSm
        SmStart=toc;			%-Timestamp (>)
        if bVolm
          TmpVol(Q) = ResSS;
          % FWHM in voxels (and not in mm) as TmpVol is not a struct 
          spm_smooth(TmpVol,SmResSS,vFWHM./VOX);
          ResSS     = SmResSS(Q)./SmMask(Q);
        else
          TmpPl(Q)  = ResSS;
          SmResSS   = spm_conv(TmpPl,vFWHM(1)/VOX(1),vFWHM(2)/VOX(2));
          ResSS     = SmResSS(Q)./SmMask(Q);
        end
        SmTime = SmTime + toc-SmStart;	%-Timestamp (>)
      end
	    
      %-Compute t-statistics for specified contrast of parameters
      %-----------------------------------------------------------
      T      = zeros(1,size(BETA,2));
      Co     = CONT;
      if STAT=='T'
        % t, as usual
        T(1,:) = Co*BETA./sqrt((ResSS*(Co*pinv([HC B G]'*[HC B G])*Co'))/df);
      else
        % F!
        pX   = pinv([HC B G]);
        T(1,:) = (sum(((Co*BETA)'*inv(Co*pinv([HC B G]'*[HC B G])*Co'))' .* ...
                (Co*BETA),1)/size(Co,1)) ./ (ResSS/df);
      end	
      
      
      %-Save Max T statistic
      %-----------------------------------------------------------
      MaxT(perm,:) = max([ max(T(1,:)), -min(T(1,:));      ...
		    MaxT(perm,1), MaxT(perm,2) ]);
	    
      %-Update nonparametric P-value
      %-----------------------------------------------------------
      if (perm==1)
        T0 = T;
        nPtmp = ones(size(T));
        if bhPerms
          nPtmp = nPtmp + (T0<=0);
        end
      else
        if bhPerms
          nPtmp = nPtmp + (T>=T0) + (-T>=T0);   % NB: Worry if T0=T=0
                  % if STAT=='T', then T, 
                  % T0 >=0, so (-T>=T0) 
                  % will be empty.
        else
          nPtmp = nPtmp + (T>=T0);
        end
      end
      
      %-Save min weighted p-value
      %-----------------------------------------------------------
      if bVarAlph,
        MinwP(perm,:) = min([ min(Wt.*(1-spm_Tcdf( T(1,:),df))),     ...
              min(Wt.*(1-spm_Tcdf(-T(1,:),df)));    ...
              MinwP(perm,1), MinwP(perm,2) ]);
      end
      
      %-Save T,XYZ,perm for suprathreshold analysis
      %-----------------------------------------------------------
      if bST 

        if pU_ST_Ut==-1  % No threshold set - save mountain tops
          clear d1 d2
          d1 = find(T(1,:) >  ST_Ut);
          d2 = find(T(1,:) < -ST_Ut);
          spm_append_96('SnPM_ST',[                            ...
          XYZ(:,d1),               XYZ(:,d2);              ...
          T(1,d1),                 -T(1,d2);               ...
          perm*ones(1,length(d1)), -perm*ones(1,length(d2)) ...
          ],'Consider using ''set cluster-forming threshold now (fast)'' option');

        else  % pU_ST_Ut>=0 - threshold set

          clear d1 d2 SnPM_ST_Pos SnPM_ST_Neg
          d1 = find(T(1,:) >  ST_Ut);
          d2 = find(T(1,:) < -ST_Ut);
          SnPM_ST_Pos=[               ...
          XYZ(:,d1);              ...
          T(1,d1)];  

          SnPM_ST_Neg=[               ...
          XYZ(:,d2);              ...
          -T(1,d2)];

          if STAT== 'F'
            loop = 1;
          else
            loop = 1:2;
          end

          for isPos = loop %1 for positive; 2 for negative
            if isPos==1
              SnPM_ST = SnPM_ST_Pos;
            else
              SnPM_ST = SnPM_ST_Neg;
            end

            % consider Permuation NO. perm
            if ~isempty(SnPM_ST)
              Locs_mm=SnPM_ST(1:3,:);
              Locs_mm (4,:) = 1;
              Locs_vox = IMAT * Locs_mm;

              % Sometimes Locs_vox are not exactly integers and this raises an
              % error later in the code. Here check that the values are
              % integers with respect to a level of absolute tolerance (~10^-14)
              % and enforce Locs_vox to be integers.
              diffWithRounded = max(abs(Locs_vox(:)-round(Locs_vox(:))));
              tolerance = 10^-10;
              if diffWithRounded > tolerance
                 Locs_vox_alter = MAT\Locs_mm;
                 diffWithRounded_alter = max(abs(Locs_vox_alter(:)-round(Locs_vox(:))));
                 error('SnPM:NonIntegerLocs', ['''Locs_vox'' must be integers (difference is ' num2str(diffWithRounded) ...
                     ' or ' num2str(diffWithRounded_alter) ')']);
              else
                 Locs_vox = round(Locs_vox); 
              end

              STCS = snpm_STcalc('update',STCS, SnPM_ST(4,:),...
              Locs_vox(1:3,:),isPos,perm,pU_ST_Ut,df);

              %save perm 1 stats for use later -[X;Y;Z;T;perm;STCno]
              if (perm==1)
                tmp = spm_clusters(Locs_vox(1:3,:));
                STCstats=[SnPM_ST;perm*ones(1,size(SnPM_ST,2));tmp];
                if isPos==1
                  save SnPM_pp STCstats
                else
                  STCstats_Neg = STCstats;
                  save SnPM_pp_Neg STCstats_Neg
                end
              end			
            end  % if ~isempty(SnPM_ST) 
          end  % for isPos=loop    
        end % pU_ST_Ut==-1 

      end % bST
      
      %-Print status at each perm if bVolm (& stop maybe)
      %-----------------------------------------------------------
      if bVolm
	if bWin, spm_progress_bar('Set',perm), end
	%-Printout timing information (>)
	fprintf('\tPerm %4d: %3d''%2d" (%3d%% Sm)\n', perm, ...
		floor((toc-PmStart)/60),round(rem(toc-PmStart,60)), ...
		round(100*SmTime/(toc-PmStart)))
	
      end % (if bVolm)

    end 	% (for perm = StartPerm:nPerm) - Perm loop
    
    %- save STCS
    if bST & pU_ST_Ut>=0
      if bhPerms %Double the STCS variables.
	STCS = snpm_STcalc('double',STCS);
      end
	     	    
      save STCS STCS
    end
    
  end 	% (length(Q)) - Conditional on non-zero voxels
  nP = [nP, nPtmp];
  nPtmp = [];
  
  if bVolm 
    break
  end
    
  %-Print status at each plane if ~bVolm
  %-----------------------------------------------------------
  if bWin, spm_progress_bar('Set',i), end
  %-Printout timing information (>)
  fprintf('%3d''%2d" (%3d%% Sm)\n', ...
	  floor((toc-PlStart)/60),round(rem(toc-PlStart,60)), ...
	  round(100*SmTime/(toc-PlStart)))
  
end		% (for i = 1:zdim) - loop over planes

fprintf('\n\nPermutations are done. Writing out images.\n')

if bhPerms
  nP = nP/(2*nPerm);
else
  nP = nP/nPerm;
end
SnPMucp=nP;
save SnPMucp SnPMucp

%
% - write out lP+ and lP- images;
%
nP_pos=nP;
lP_pos=-log10(nP_pos);
lP_pos_image(spm_xyz2e(XYZ_total, Vt))=lP_pos;
lP_pos_vol=reshape(lP_pos_image,DIM(1),DIM(2),DIM(3));
spm_write_vol(VlP_pos, lP_pos_vol);

if STAT == 'T'
  if bhPerms
    nP_neg=1+1/(2*nPerm)-nP;
  else 
    nP_neg=1+1/nPerm-nP;
  end

  lP_neg=-log10(nP_neg);
  lP_neg_image(spm_xyz2e(XYZ_total, Vt))=lP_neg;
  lP_neg_vol=reshape(lP_neg_image,DIM(1),DIM(2),DIM(3));
  spm_write_vol(VlP_neg, lP_neg_vol);
end

%
% - write out lP_FWE+ and lP_FWE- images;
%
tol = 1e-4;	% Tolerance for comparing real numbers

cP_pos=zeros(size(nP));

if bhPerms
	MaxT   = [ MaxT; flipud(fliplr(MaxT)) ];
end
MaxT_pos=MaxT(:,1);

for t = MaxT_pos'
	%-FEW-corrected p is proportion of randomisation greater or
	% equal to statistic.
	%-Use a > b -tol rather than a >= b to avoid comparing
	% two reals for equality.
	cP_pos = cP_pos + (t > SnPMt -tol);
end

if STAT =='T'
  cP_neg=zeros(size(nP));
  MaxT_neg=MaxT(:,2);

  for t = MaxT_neg'
        cP_neg = cP_neg + (t > -SnPMt -tol);
  end
end

if bhPerms
  cP_pos = cP_pos / (2* nPerm);
  if STAT=='T'	
    cP_neg = cP_neg / (2* nPerm);
  end
else 
  cP_pos = cP_pos / nPerm;  
  if STAT=='T'
    cP_neg = cP_neg / nPerm;
  end
end

lP_FWE_pos=-log10(cP_pos);
lP_FWE_pos_image(spm_xyz2e(XYZ_total, Vt))=lP_FWE_pos;
lP_FWE_pos_vol=reshape(lP_FWE_pos_image,DIM(1),DIM(2),DIM(3));
spm_write_vol(VlP_FWE_pos, lP_FWE_pos_vol);

if STAT=='T'
  lP_FWE_neg=-log10(cP_neg);
  lP_FWE_neg_image(spm_xyz2e(XYZ_total, Vt))=lP_FWE_neg;
  lP_FWE_neg_vol=reshape(lP_FWE_neg_image,DIM(1),DIM(2),DIM(3));
  spm_write_vol(VlP_FWE_neg, lP_FWE_neg_vol);
end

%
% - write out lP_FDR+ and lP_FDR- images;
%
[snP_pos,I_pos]=sort(nP_pos);

Pfdr_pos=snpm_P_FDR([],[],'P',[],snP_pos');

Pfdr_pos(I_pos) = Pfdr_pos;

lP_FDR_pos=-log10(Pfdr_pos);

lP_FDR_pos_image(spm_xyz2e(XYZ_total, Vt))=lP_FDR_pos;

lP_FDR_pos_vol=reshape(lP_FDR_pos_image,DIM(1),DIM(2),DIM(3));
spm_write_vol(VlP_FDR_pos, lP_FDR_pos_vol);

if STAT =='T'
  [snP_neg,I_neg]=sort(nP_neg);
  Pfdr_neg=snpm_P_FDR([],[],'P',[],snP_neg');
  Pfdr_neg(I_neg) = Pfdr_neg;
  lP_FDR_neg=-log10(Pfdr_neg);
  lP_FDR_neg_image(spm_xyz2e(XYZ_total, Vt))=lP_FDR_neg;

  lP_FDR_neg_vol=reshape(lP_FDR_neg_image,DIM(1),DIM(2),DIM(3));
  spm_write_vol(VlP_FDR_neg, lP_FDR_neg_vol);
end


if bWin, spm_progress_bar('Clear'), end
%-Printout final timing information (>)
fprintf('\n\nThe run took %0.02f minutes\n', toc/60);

%-Cleanup
clear X

%-Save key variables
%=======================================================================
eval(['save SnPM ',s_SnPM_save])

%-Print quick summary info (allowing for STOPing)
%=======================================================================
if bhPerms
    Rank = sum([MaxT(1:perm,1);MaxT(1:perm,2)] >= MaxT(1,1));
else
    Rank = sum([MaxT(1:perm,1)] >= MaxT(1,1));
end
fprintf(['\nCorrect Perm has max t %g & rank %d out of %d ', ...
	'completed permutations\n'],MaxT(1,1),Rank,perm*(bhPerms+1));
fprintf('\n\tRun snpm_pp for full results\n\n');



function Vs = sf_close_vol(Vs)
% Don't need to close images in SPM5
return


