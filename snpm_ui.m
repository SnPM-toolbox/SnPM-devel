function snpm_ui(varargin)
% Set up for general linear model and permutation/rando analysis
% FORMAT snpm_ui
%        snpm_ui(job)
%_______________________________________________________________________
%
% snpm_ui sets up the parameters for a non-parametric
% permutation/randomisation analysis. The approach taken with SnPM
% analyses differs from that of SPM. Instead of intiating the analysis
% with one command, an analysis consists of 3 steps:
%
%   1. Configuration of Analysis           --- interactive
%   2. Calculation of Raw Statistics       --- noninteractive
%   3. Post Processing of Raw Statistics   --- interactive
%
% ( In SPM, all of these steps are done with a click to the "Statistics"  )
% ( button, though the 3rd step is often redone with "Results" or "SPM{Z}")
%
% The first step is embodied in this function, snpm_ui. snpm_ui
% configures the design matrix and calls "plug in" modules that specify
% how the relabeling is to be done for particular designs.
%
% The result of this function is a mat file, "SnPMcfg.mat", written
% to the present working directory (when using the batch system, the
% present directory is changed to the analysis directory specified). 
% This file contains all the
% parameters needed to perform the second step, which is in embodied in
% snpm_cp. Design parameters are displayed in the SPM graphics window,
% and are printed.
%
% When called with a job file, as configured by the SPM batch system, no
% user interaction is required.
%-----------------------------------------------------------------------
%
%-The Prompts Explained
%=======================================================================
%
% 'Select design type...': Choose from the available designs. Use the
% 'User Specifed PlugIn' option to supply your own PlugIn function. Use
% the 'keyboard' option to manually set the required PlugIn variables
% (as defined below under "PlugIn Must Supply the following"). The "User
% Specified" option allows you to select a custom PlugIn file. Note that
% due to a limitation in MatLab, only the name of the PlugIn is used,
% not the directory, so make sure the PlugIn is the first with that name
% on the MATLABPATH. type which <plugin-file-name> at the MatLab prompt
% before starting snpm_ui to check!
%
% - At this point you will be prompted by the PlugIn file;
% - see help for the PlugIn file you selected.
%
% 'FWHM(mm) for Variance smooth': Variance smoothing gives the
% nonparmetric approach more power over the parametric approach for
% low-df analyses.  If your design has fewer than 20 degrees of freedom
% variance smoothing is advised.  10 mm FWHM is a good starting point
% for the size of the smoothing kernal. For non-isotropic smoothing,
% enter three numbers: FWHM(x) FWHM(y) FWHM(z), all in millimeters.
%
% If there are enough scans and there is no variance smoothing in the z
% direction, you will be asked...
%
% '# scans: Work volumetrically?':  Volumetric means that the entire
% data set is loaded into memory; while this is more efficient than
% iterating over planes it is very memory intensive.
%
% ( Note: If you specify variance smoothing in the z-direction, SnPM     )
% ( (in snpm_cp.m) has to work volumetrically. Thus, for moderate to    )
% ( large numbers of scans there might not be enough memory to complete  )
% ( the calculations. This shouldn't be too much of a problem because    )
% ( variance smoothing is only necesary at low df, which usually         )
% ( corresponds to a small number of scans. Alternatively, specify only  )
% ( in-plane smoothing.                                                  )
%
% 'Collect Supra-Threshold stats?': In order to use the permutation
% test on supra-threshold cluster size you have to collect a
% substantial amount of additional data for each permutation.  If you
% want to look at cluster size answer yes, but have lots of free space
% (the more permutations, the more space needed). You can however
% delete the SnPM_ST.mat file containing the supra-threshold cluster
% data at a later date without affecting the ability to analyze at the
% voxel level.
%
% The remaining questions all parallel the standard parametric analysis
% of SPM; namely
% 'Select global normalisation' - account for global flow confound
% 'Select global calculation'   - how global flow is computed
% 'Threshold masking'           - threshold to determine voxels to analyze
% 'Explicit masking'            - mask image to set analysis
% 'Grand mean scaling'          - Whether to scale overall grand mean
%
%
% Programmers / Hackers help...
%=======================================================================
%
% Variables saved in SnPMcfg.mat
%-----------------------------------------------------------------------
% H             condition partition of DesMtx for correctly labeled data
% C             covariate partition of DesMtx for correctly labeled data
% B             block     partition of DesMtx for correctly labeled data
% G             confound  partition of DesMtx for correctly labeled data
% HCBGnames     string matrix of column names of [H C B G]
% P             string matrix of Filenames corresponding to observations
% PiCond        Permuted conditions matrix, one labelling per row, actual
%               labelling on first row
% sPiCond       String describing permutations in PiCond
% bhPerms       Flag indicating use of "half permutations" trick
% sHCform       String for computation of HC design matrix partitions
%               permutations indexed by perm in snpm_cp
% iGloNorm      Global normalisation code
% sGloNorm      String describing Global Normalisation option
% GM            Grand mean
% GMscale       Scaling coefficients for grand mean scaling
% GX            Global means for each scan (after scaling of Grand Mean)
% CONT          Contrast (only one)
% THRESH        Grey matter threshold, as %age of global mean
% TH            Grey matter thresholds for each image
% bVarSm        Flag for variance smoothing (Pseudo t-statistics)
% vFWHM         FWHM for variance smoothing ([0,0,0] if bVarSm=0)
% sVarSm        Sring describing variance Smoothing (empty if bVarSm=0)
% bVolm         Flag for volumetric computation (whole volume at once)
% bST           Flag for collection of superthreshold info 
% sDesFile      Name of PlugIn design file
% sDesign       Description of PlugIn design
% V             Memory mapping handles
% MASK          Filename of explicit mask image
% ImMASK        Implicit masking; 0=none; 1=zeros are equivalent to NaN
% nidm_json     json structure storing minimal information required for a
%               NIDM-Results export
% 
% df            degrees of freedom due to error
% sDesSave      String of PlugIn variables to save to cfg file
%               String itself and variables listed within it are saved.
% s_SnPMcfg_save string matrix of all variables saved in SnPMcfg
%
% PlugIn Must Supply the following:
%-----------------------------------------------------------------------
% P             - string matrix of Filenames corresponding to observations
% iGloNorm      - Global normalisation code, or allowable codes
%               - Names of columns of design matrix subpartitions
% PiCond        - Permuted conditions matrix, one labelling per row, actual
%                 labelling on first row
% sPiCond       - String describing permutations in PiCond
% sHCform       - String for computation of HC design matrix partitions
%                 permutations indexed by perm in snpm_cp
% CONT          - single contrast for examination, a row vector
% sDesign       - String defining the design [Defaults to PlugIn description]
% sDesSave      - String of PlugIn variables to save to cfg file  [Default '']
%
% In addition, non-null portions of the design matrix (i.e. non-null H
% C B G) must be specified along with the respective effect names
% (Hnames, Cnames &c.)
%
%_______________________________________________________________________
% Copyright (C) 2013 The University of Warwick
% Id: snpm_ui.m  SnPM13 2013/10/12
% Thomas Nichols, Camille Maumet

%------------------------------functions-called------------------------
% spm_DesMtx
% spm_clf
% spm_figure
% spm_select
% spm_global
% spm_hread
% spm_input
% snpm_MSA2x
% snpm_SSA2x
% snpm_SSC
% spm_str_manip
% spm_vol
%----------------------------- functions-called------------------------

%-Initialise workspace
%-----------------------------------------------------------------------
global defaults
if isempty(defaults), spm_defaults; end
global SnPMdefs
if isempty(SnPMdefs), snpm_defaults; end

if nargin==0
  % Traditional mode
  Finter = spm_figure('FindWin','Interactive');
  Fgraph = spm_figure('FindWin','Graphics');
  if isempty(Fgraph), Fgraph=spm_figure('Create','Graphics'); end
  spm_clf(Finter), spm_clf(Fgraph)
  set(Finter,'Name','SnPM Setup');
  BATCH=false;
else
  job=varargin{1};
  BATCH=true;
  if exist(job.dir{1})~=7
    mkdir(job.dir{1})
  end
  cd(job.dir{1})
end

nidm_json = containers.Map();

%-Definitions & Design parameters
%=======================================================================
sDesigns=str2mat(...
	'SingleSub: Two Sample T test; 2 conditions',...
	'SingleSub: Simple Regression; 1 covariate of interest',...
	'MultiSub: One Sample T test on differences; 1 condition',...
	'MultiSub: Simple Regression; 1 covariate of interest',...
	'MultiSub: Paired T test; 2 conditions, 1 scan/condition',...
	'MultiSub: Randomized Design Paired T test; 2 conditions',... 
        'MultiSub: Within Subject ANOVA; multiple scans/subject',... 
	'2 Groups: Test diff of response; 2 conditions, 1 scan/condition',...
	'2 Groups: Two Sample T test; 1 scan/subject',...
	'>2 Groups: Between Group ANOVA; 1 scan/subject',...
        'User Specified PlugIn',...
	'keyboard');

sDesFile=str2mat(...
        'snpm_pi_TwoSampTss', ...
	'snpm_pi_Corr1S', ...
	'snpm_pi_OneSampT', ...     
	'snpm_pi_Corr', ...         
	'snpm_pi_PairT', ...        
	'snpm_pi_PairTrand', ...
        'snpm_pi_ANOVAwithinS',... 
	'snpm_pi_TwoSampPairT', ... 
	'snpm_pi_TwoSampT', ...     
	'snpm_pi_ANOVAbtwnS',...   
    	'',...
	['clc, fprintf([''SnPM: Setup PlugIn variables by hand...\n',...
		'\tsee snpm_ui for PlugIn interface definitions\n',...
		'\tType return when done\n\n'']), keyboard']);

%%%
%%% These descriptions and codes must match those in config/snpm_cfg_ui.m
%%%

%-Global normalization                                    (GloNorm)
sGloNorm=str2mat(... 
	'<no global normalisation>',...				%-1
	'proportional scaling',...				%-2
	'AnCova',...						%-3
	'AnCova {subject-specific}',...				%-4
	'AnCova {study-specific}');				%-5

%-Global calculation options                               (GXcalc)
sGXcalc  = str2mat(...
    'omit',...							%-1
    'user specified',...					%-2
    'mean voxel value (within per image fullmean/8 mask)');	%-3

%-Grand mean scaling options                                (GMsca)
sGMsca = str2mat(...
    '<no grand Mean scaling>',...				%-1
    'scaling of overall grand mean');				%-2


%-Select design type
%=======================================================================
if BATCH
  sDesign=job.DesignName;
  sDesFile=strrep(job.DesignFile,'bch_ui','pi');
else
  DesType  = spm_input('Select design type...',1,'m',sDesigns);
  sDesign  = deblank(sDesigns(DesType,:));
  sDesFile = deblank(sDesFile(DesType,:));
  if isempty(sDesFile)
    sDesFile = spm_select(1,'^.*\.m$','Select SnPM design PlugIn Mfile...');
    sDesFile = spm_str_manip(sDesFile,'rtd');
  end
end

%-Variable initialisation prior to running PlugIn
%-----------------------------------------------------------------------
iStud=[];	% Study indicator vector
iSubj=[];	% Subject indicator vector
iCond=[];	% Condition indicator vector
iRepl=[];	% Replication indicator vector
iXblk=[];	% Exchangability block indicator vector
H=[]; Hnames='';% Condition partition & effect names
C=[]; Cnames='';% Covariates (of interest)
Cc=[];		% Covariates (of interest)	| Required only for covariate
Ccnames=[];	% Names of covariates		| by factor interactions
B=[]; Bnames=[];% Block partition & effect names
G=[];Gnames='';	% Covariates (no interest)
Gc=[];		% Covariates (no interest)	| Required only for covariate
Gcnames=[];	% Names of covariates		| by factor interactions
bST=0;		% Flag for collection of superthreshold info
bVarSm=0;	% Flag for variance smoothing
vFWHM=[0,0,0];	% FWHM for variance smoothing
sVarSm='';	% String describing Variance Smoothing
bVolm=0;	% Flag for volumetric computation
nMax4DefVol=snpm_get_defaults('nMax4DefVol');
                % Default to volumetric if less than nMax4DefVol scans
sPiCond='';	% String describing permutations in PiCond
bhPerms=0;	% Flag for half permutations. Rest are then their inverses
sDesSave='';	% String of PlugIn variables to save to cfg file
iGXcalc='23';   % Global calculation
iGMsca = '12';  % Grand mean scaling of globals
df1=1;          % For F stat, it will be the numerator df; For T stat, it
                % will be 1.    

%-Run PlugIn design specification module
%=======================================================================
eval(sDesFile);

%-Total #observations
nScan = size(P,1);

%-Get general analysis & data parameters
%=======================================================================

%-Ask about variance smoothing & volumetric computation
%-----------------------------------------------------------------------
vFWHM = job.vFWHM;
if length(vFWHM)==1
	vFWHM = vFWHM * ones(1,3);
elseif length(vFWHM)==2
	vFWHM = [vFWHM, 0];
else
	vFWHM = reshape(vFWHM(1:3),1,3);
end
if ~all(vFWHM==0)
  bVarSm=1; 
end
if bVarSm
  sVarSm = sprintf('Pseudo-t: Variance smoothed with FWHM [%dx%dx%d]  mm',vFWHM);
end

%-Decide upon volumetric operation
bVolm = job.bVolm;
if ~bVolm & (vFWHM(3)~=0)
warning('SnPM:MayRunOutOfMemory', ...
    sprintf(['Working volumetrically because of smoothing in z (%g).\n'... 
         'May run out of memory.'],vFWHM(3)));
bVolm=1;
end

%-Ask about collecting Supra-Threshold cluster statistics
%-----------------------------------------------------------------------
bST = ~isfield(job.ST,'ST_none');

% Add: get primary threshold for STC analysis if requested
if bST
  if ~bVolm
    warning('SnPM:CannotDefineThreshVolumetrically', ...
        sprintf('Note:  Cannot define threshold now, because not working volumetrically\n'));
    pU_ST_Ut=-1; % Define the threshold later
  else
    pU_ST_Ut = ~isfield(job.ST,'ST_later');
    if pU_ST_Ut    % Define the threshold right now
        pU_ST_Ut = job.ST.ST_U;
    else
      pU_ST_Ut=-1; % Define the threshold later
    end
  end
else
  pU_ST_Ut=NaN; 
end  


%-Global normalization options
%-----------------------------------------------------------------------
iGloNorm = job.globalm.glonorm;
sGloNorm=deblank(sGloNorm(iGloNorm,:));


%-Get value to be assigned to grand mean:
%-----------------------------------------------------------------------
if iGloNorm==2 % Proportional scaling
  iGMsca = 2;%CHANGED from 1 to 2	%-grand mean scaling implicit in PropSca GloNorm
  % TODO: must be set otherwise prop scaling alone in Batch mode lead to error  
  % job.globalm.gmsca.gmscv = ; % By default use mean
else % No global normalisation or ANCOVA
    iGMsca = 1 + isfield(job.globalm.gmsca,'gmsca_yes');
end

if (iGMsca==2) % CHANGED from 1 to 2 as should not ask for a value if grand mean scaling is not required.
  nidm_json('nidm_Data/nidm_grandMeanScaling') = true;
  if (iGloNorm==2) % Proportional scaling
    str = 'PropSca global mean to';
  else
    str = [strrep(sGMsca(iGMsca,:),'scaling of','scale'),' to'];
  end
  % As in SPM:
  switch char(fieldnames(job.globalm.gmsca))
      case 'gmsca_yes',
          % Proportionally scale to this value
          GM = job.globalm.gmsca.gmsca_yes.gmscv;
      case 'gmsca_no',
          GM = 50;
  end
  nidm_json('nidm_Data/nidm_targetIntensity') = GM;
elseif (iGMsca==1) % No grand mean scaling
  GM = 0;
  nidm_json('nidm_Data/nidm_grandMeanScaling') = false;
end


%-Get globals
%-----------------------------------------------------------------------
if (iGloNorm==1) && (iGMsca==1)
    % No need for globals, omit
    iGXcalc = 1;
elseif isfield(job.globalc,'g_omit')
    %iGXcalc=1;
    % Global needed and none requested, force iGXcalc to grand mean (3)
    iGXcalc=3;
elseif isfield(job.globalc,'g_user')
    iGXcalc=2;
    GX = job.globalc.g_user.global_uval;
    rg = GX;
    if length(GX) ~= nScan
        error('SnPM:InvalidGlobals', ['User-specified globals length [%d] doesn''t match number of' ...
            ' scans [%d]'],length(GX),nScan);
    end
elseif isfield(job.globalc,'g_mean')
    iGXcalc=3;
end
sGXcalc=deblank(sGXcalc(iGXcalc,:));


%-Get threshold defining voxels to analyse
%-----------------------------------------------------------------------
if isfield(job.masking.tm,'tm_none')
    iTHRESH=1;
    THRESH = -Inf;
    sThresh = 'None';
elseif isfield(job.masking.tm,'tmr')
    iTHRESH=2;
    THRESH=job.masking.tm.tmr.rthresh;
    sThresh=sprintf('Relative (%g)',THRESH);
    
    % As in SPM if relative thresholding then calculation of globals is
    % implied
    if iGXcalc~=2
        iGXcalc = 3;
    end
elseif isfield(job.masking.tm,'tma')
    iTHRESH=3;
    THRESH=job.masking.tm.tma.athresh;
    sThresh = sprintf('Absolute (%g)',THRESH);
end

%-Implicit Masking - Batch only!
%-----------------------------------------------------------------------
ImMASK=job.masking.im;

%-Get analysis mask
%-----------------------------------------------------------------------
MASK = job.masking.em{1};

%=======================================================================
%-Computation
%=======================================================================

%-Condition Cc & Gc "Shadow" partitions if no FxC interactions
% These store the covariate values for printing only
%-----------------------------------------------------------------------
if (isempty(Cc) & ~isempty(C)), Cc=C; Ccnames=Cnames; end
if (isempty(Gc) & ~isempty(G)), Gc=G; Gcnames=Gnames; end


%-Examine images
%=======================================================================

%-MMap image files
V = spm_vol(P);

%-Check compatability of images
spm_check_orientations(V);

%-Get ORIGIN, etc
DIM    = [V(1).dim(1)   V(1).dim(2)   V(1).dim(3)]';
M=V(1).mat(1:3, 1:3);
VOX=sqrt(diag(M'*M))';
MAT    = V(1).mat;
IMAT   = inv(MAT);
ORIGIN = IMAT(1:3,4);

%-Global calculation
%-----------------------------------------------------------------------
if iGXcalc==2
  %-User specified globals
elseif iGXcalc==3
  %-Compute global values
  rg     = zeros(nScan,1);
  for i  = 1:nScan, rg(i)=spm_global(V(i)); end
  if any(~isfinite(rg))
    disp(rg)
    error('SnPM:NaNGlobal', 'Global computation returned NaN! Cannot continue')
  end
  GX     = rg;
elseif iGXcalc==1
  rg     = [];
  GX     = [];
end



%-Scale scaling coefficients so that Grand mean, mean(GX), is = GM (if GM~=0)
% Since images are unmapped, this must be replicated in snpm_cp
% Done here to provide check on V's in snpm_cp
if GM ~= 0
	GMscale = GM/mean(GX);
	for i = 1:nScan, 
		V(i).pinfo(1:2,:)  = V(i).pinfo(1:2,:) * GMscale;
	end
	GX      = GX     * GMscale;
else
	GMscale = 1;
end

%-Compute Grey matter threshold for each image
if (iTHRESH==1) || (iTHRESH==3)
    % No threshold or Absolute threshold    
    TH    = repmat(THRESH,nScan,1);
elseif (iTHRESH==2)
    % Relative threshold
    TH    = THRESH * GX;
else 
    error('SnPM:InvalidiTHRESH', ['Wrong value for iTHRESH: ', num2str(iTHRESH)])
end


%-Construct Global part of covariates of no interest partition.
%-Centre global means if included in AnCova models, by mean correction.
%=======================================================================
Gc    = [Gc,GX];
if isempty(Gcnames), Gcnames = 'Global';
    else Gcnames = str2mat(Gcnames,'Global'); end

if iGloNorm == 1				%-No global adjustment
%-----------------------------------------------------------------------
elseif iGloNorm == 2				%-Proportional scaling
%-----------------------------------------------------------------------
% Since images are unmapped, this must be replicated in snpm_cp
% Done here to provide check on V's in snpm_cp
   for i = 1:nScan,
      V(i).pinfo(1:2,:) = GM*V(i).pinfo(1:2,:)/GX(i);
   end

elseif iGloNorm == 3				%-AnCova
%-----------------------------------------------------------------------
   G = [G,(GX - mean(GX))];
   if isempty(Gnames), Gnames = 'Global'; 
       else Gnames = str2mat(Gnames,'Global'); end

elseif iGloNorm == 4				%-AnCova by subject
%-----------------------------------------------------------------------
    [GL,GLnames] = spm_DesMtx([iSUBJ',GX-mean(GX)],'FxC',['SUBJ  ';'Global']);
    G = [G,GL];
    if isempty(Gnames), Gnames = GLnames;
        else Gnames = str2mat(Gnames,GLnames); end

elseif iGloNorm == 5				%-AnCova by study
%-----------------------------------------------------------------------
    [GL,GLnames] = spm_DesMtx([iStud',GX-mean(GX)],'FxC',['Stud  ';'Global']);
    G = [G,GL];
    if isempty(Gnames), Gnames = GLnames; 
        else Gnames = str2mat(Gnames,GLnames); end
else
%-----------------------------------------------------------------------
    error('SnPM:InvalidiGloNorm', sprintf('%cError: invalid iGloNorm option\n',7))
end % (if)


%=======================================================================


%-Ensure validity of contrast of condition effects, zero pad
%-----------------------------------------------------------------------
%-Only a single contrast
if size(CONT,1)==1
    if size(H,2)>1
        CONT(1:size(H,2)) = CONT(1:size(H,2)) - mean(CONT(1:size(H,2)));
    end
end
%-Zero pad for B & G partitions.
% (Note that we trust PlugIns to create valid contrasts for [H C])
CONT  = [CONT, zeros(size(CONT,1),size([B G],2))];

%-Construct full design matrix and name matrices for display
%-----------------------------------------------------------------------
[nHCBG,HCBGnames] = spm_DesMtx('Sca',H,Hnames,C,Cnames,B,Bnames,G,Gnames);
nidm_json('nidm_DesignMatrix/prov:value') = nHCBG;
nidm_json('nidm_DesignMatrix/nidm_regressorNames') = HCBGnames;

%-Setup is complete - save SnPMcfg Mat file
%-----------------------------------------------------------------------
s_SnPMcfg_save = ['s_SnPMcfg_save H C B G HCBGnames P PiCond ',...
	'sPiCond bhPerms sHCform iGloNorm sGloNorm GM rg GX GMscale CONT ',...
	'THRESH MASK ImMASK TH bVarSm vFWHM sVarSm bVolm bST sDesFile sDesign ',...
        'V pU_ST_Ut df1 nidm_json nPiCond_mx nidm_json ', ...
	'sDesSave ',sDesSave];
eval(['save SnPMcfg ',s_SnPMcfg_save])


if ~spm_get_defaults('cmdline')

    %=======================================================================
    %-Display parameters
    %=======================================================================

    if BATCH
      % OK, only _now_ activate the graphics... so if we crash here, there's
      % nothing lost
      Finter = spm_figure('FindWin','Interactive');
      Fgraph = spm_figure('FindWin','Graphics');
      if isempty(Fgraph), Fgraph=spm_figure('Create','Graphics'); end
      spm_clf(Finter), spm_clf(Fgraph)
    end


    %-Muck about a bit to set flags for various indicators - handy for later
    bMStud=~isempty(iStud);
    bMSubj=~isempty(iSubj);
    bMCond=~isempty(iCond);
    bMRepl=~isempty(iRepl);
    bMXblk=~isempty(iXblk);

    %-Compute common path components - all paths will begin with file separator
    %-----------------------------------------------------------------------
    d     = max(find(P(1,1:find(~all(P == ones(nScan,1)*P(1,:)), 1 )-1)==filesep)) - 1;
    CPath = P(1,1:d);
    Q     = P(:,d+1:size(P,2));

    %-Display data parameters
    %=======================================================================
    figure(Fgraph); spm_clf; axis off
    text(0.30,1.02,'Statistical analysis','Fontsize',16,'Fontweight','Bold');
    text(-0.10,0.85,'Scan Index','Rotation',90)
    if bMStud, text(-0.05,0.85,'Study',      'Rotation',90); end
    if bMSubj, text(+0.00,0.85,'Subject',    'Rotation',90); end
    if bMCond, text(+0.05,0.85,'Condition',  'Rotation',90); end
    if bMRepl, text(+0.10,0.85,'Replication','Rotation',90); end
    if bMXblk, text(+0.15,0.85,'Exchange Blk','Rotation',90); end
    x0    = 0.20; y0 = 0.83;
    dx    = 0.15; dy = 0.02;
    x     = x0;
    for i = 1:size(Cc,2)
        text(x + 0.02,0.85,Ccnames(i,:),'Rotation',90);
        x = x + dx; end
    for i = 1:size(Gc,2)
        text(x + 0.02,0.85,Gcnames(i,:),'Rotation',90);
    x = x + dx; end
    text(x,0.92,'Base directory:','FontSize',10,'Fontweight','Bold');
    text(x,0.90,CPath,'FontSize',10,'interpreter','none');
    text(x,0.87,'Filename Tails');
    y     = y0;

    for i = 1:nScan
        text(-0.12,y,sprintf('%02d :',i));
       if bMStud, text(-0.06,y,sprintf('%2d',iStud(i))); end
       if bMSubj, text(-0.01,y,sprintf('%2d',iSubj(i))); end
       if bMCond, text(+0.04,y,sprintf('%2d',iCond(i))); end
       if bMRepl, text(+0.09,y,sprintf('%2d',iRepl(i))); end
       if bMXblk, text(+0.14,y,sprintf('%2d',iXblk(i))); end
       x     = x0;
       for j = 1:size(Cc,2)
        text(x,y,sprintf('%-8.6g',Cc(i,j)),'FontSize',10)
        x = x + dx; end
       for j = 1:size(Gc,2)
        text(x,y,sprintf('%-8.6g',Gc(i,j)),'FontSize',10)
        x = x + dx; end
       text(x,y,Q(i,:),'FontSize',10,'interpreter','none');
       y     = y - dy;
       if y < 0;
        spm_print
        spm_clf; axis off
        y = y0;
        text(0.16,1.02,['Statistical analysis (continued)'],...
            'Fontsize',16,'Fontweight','Bold');
       end
    end

    %-Print miscellaneous data parameters
    %-----------------------------------------------------------------------
    y      = y - dy;
    dy     = dy*1.2;
    if (GM~=0)
        text(0,y,sprintf(['Images scaled to a grand mean of %g'],GM))
        y = y - dy;
    end
    text(0,y,sprintf(...
        'Analysis threshold is %3.0f%% of the whole brain mean',THRESH*100))
    spm_print


    %-Display design parameters
    %=======================================================================
    figure(Fgraph); spm_clf(Fgraph); axis off
    text(0.30,1.02,'Design Matrix','Fontsize',16,'Fontweight','Bold');

    %-Label the effects
    %-----------------------------------------------------------------------
    hDesMtx = axes('Position',[0.2 0.3 0.6 0.5]);
    image((nHCBG + 1)*32);
    ylabel('Observations')
    set(hDesMtx,'XTick',[],'XTickLabel','')
    hEfLabs = axes('Position',[0.2 0.82 0.6 0.1],'Visible','off');
    y     = 0.1;
    dx    = 1/size(nHCBG,2);
    for i = 1:size(nHCBG,2)
        text((i - 0.5)*dx,y,deblank(HCBGnames(i,:)),...
        'Fontsize',8,'Rotation',90)
    end


    %-Display non-parametric analysis summary
    %-----------------------------------------------------------------------
    hPramAxes=axes('Position',[0.05 0.08 0.8 0.20],'Visible','off');
    text(0,1.00,sDesign,'Fontsize',10);
    text(0,0.90,['SnPM design flie: ',sDesFile],'Fontsize',10);
    text(0,0.80,sPiCond,'Fontsize',10);
    text(0,0.70,['Global normalisation: ',deblank(sGloNorm)],'Fontsize',10);
    text(0,0.60,['Threshold masking: ',deblank(sThresh)],'Fontsize',10);

    %-Display parameter summary
    %-----------------------------------------------------------------------
    text(0,.5,'Parameters:','Fontsize',10,'Fontweight','Bold');
    text(0,.4,sprintf(['%d Condition + %d Covariate ',...
        '+ %d Block + %d Confound'],...
        size(H,2),size(C,2),size(B,2),size(G,2)),...
        'Fontsize',10);
    text(0,.3,sprintf(['= %d parameters, having %d degrees of freedom, ',...
        'giving %d residual df (%d scans).'],...
        size([H C B G],2),rank([H C B G]),nScan-rank([H C B G]),nScan),...
        'Fontsize',10);
    if (bVarSm) text(0,0.2,sVarSm,'Fontsize',10);
    end

    spm_print

    %-Clear interactive window
    %-----------------------------------------------------------------------
    spm_clf(Finter)
end