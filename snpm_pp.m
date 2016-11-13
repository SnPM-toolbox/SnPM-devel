function snpm_pp(CWD,varargin)
% SnPM post processing and results display
% FORMAT snpm_pp(CWD)
%
% CWD -	Directory containing SnPM results files
%
% If CWD is not specified then user is prompted to locate results file SnPM.mat
%_______________________________________________________________________
%
% snpm_pp is the PostProcessing function for the SnPM nonParametric
% statistical analysis. SnPM statistical analyses are split into three
% stages; Setup, Compute & Assess. This is the third stage.
% Nonparametric randomisation distributions are read in from MatLab
% *.mat files, with which the observed statistic image is assessed
% according to user defined parameters. It is the SnPM equivalent of
% the "Results" section of SPM, albeit with reduced features.
%
% Voxel level corrected p-values are computed from the permutation
% distribution of the maximal statistic. If suprathreshold cluster
% statistics were collected in the computation stage (and the large
% SnPM_ST.mat file hasn't been deleted!), then assessment by
% suprathreshold cluster size is also available, using a user-specified
% primary threshold.
%
% Instructions:
%=======================================================================
%
% You are prompted for the following:
%
% (1) ResultsDir: If the results directory wasn't specified on the command
%     line, you are prompted to locate the SnPM results file SnPM.mat.
%     The directory in which this file resides is taken to be the
%     results directory, which must contain *all* the files listed
%     below ("SnPM files required").
%
%     Results (spm.ps & any requested image files) are written in the
%     present working directory, *not* the directory containing the
%     results of the SnPM computations.
%
%                           ----------------
%
% (2) +/-: Having located and loaded the results files, you are asked to
%     chose between "Positive or negative effects?". SnPM, like SPM,
%     only implements single tailed tests. Choose "+ve" if you wish to
%     assess the statistic image for large values, indicating evidence
%     against the null hypothesis in favour of a positive alternative
%     (activation, or positive slope in a covariate analysis).
%
%     Choose "-ve" to assess the negative contrast, i.e. to look for
%     evidence against the null hypothesis in favour of a negative
%     alternative (de-activation, or a negative slope in a covariate
%     analysis). The "-ve" option negates the statistic image and
%     contrast, acting as if the negative of the actual contrast was
%     entered.
%
%     A two-sided test may be constructed by doing two separate
%     analyses, one for each tail, at half the chosen significance
%     level, doubling the resulting p-values.
%     ( Strictly speaking, this is not equivalent to a rigorous two-sided )
%     ( non-parametric test using the permutation distribution of the     )
%     ( absolute maximum statistic, but it'll do!                         )
%
%                           ----------------
%
% (3) WriteFiles: All image files have been written by snpm_cp, and there
% is only one file that may be written here:  Users are asked 'write
% filtered statistic image?'. 
%                           ----------------
%
% Next come parameters for the assessment of the statistic image...
%
% (4) alpha: (p-value for filtering)
%     First, you are asked 'Use corrected threshold?' You can choose to set
%     the threshold as FWE-corrected('FWE'), FDR-corrected('FDR') or
%     uncorrected('None').  A FWE threshold controls the chance of one or
%     more false positives; a FDR threshold controls the expected
%     fraction of false positives among the voxels detected. 
%  
%     Then you specify the threshold, the statistical significance
%     level at which you wish to assess the evidence against the null
%     hypothesis. In SPM this is called "filtering by corrected
%     p-value".  SnPM will only show you voxels (& suprathreshold
%     regions if you choose) that are significant based on the method you
%     choose at level \alpha.  I.e. only voxels (& regions) with
%     corrected (or uncorrected) p-value less than \alpha are shown to 
%     you.
%
%     Setting \alpha to 1 will show you all voxels with a positive statistic.
%
%
% (5) SpatEx: If you collected supra-threshold cluster statistics during
%     the SnPM computation phase, you are offered the option to assess
%     the statistic image by supra-threshold cluster size (spatial
%     extent).
%
% 5a) ST_Ut: If you chose to asses spatial extent, you are now prompted
%     for the primary threshold. This is the threshold applied to the
%     statistic image for the identification of supra-threshold
%     clusters.
%
%     The acceptable range is limited.  SnPM has to collect
%     suprathreshold information for every relabelling. Rather that
%     pre-specify the primary threshold, information is recorded for
%     each voxel exceeding a low threshold (set in snpm_cp) for every
%     permutation. From this, suprathreshold cluster statistics can be
%     generated for any threshold higher than the low recording
%     threshold. This presents a lower limit on the possible primary
%     threshold.
%
%     The upper limit (if specified) corresponds to the statistic value
%     at which voxels become individually significant at the chosen
%     level (\alpha).  There is little point perusing a suprathreshold
%     cluster analysis at a threshold at which the voxels are
%     individually significant.
%
%     If the statistics are t-statistics, then you can also specify the
%     threshold via the upper tail probability of the t-distribution.
%
%     (NB: For the moment, \alpha=1 precludes suprathreshold analysis, )
%     (    since all voxels are significant at \alpha=1.               )
%
%
% That's it. SnPM will now compute the appropriate significances,
% reporting its progress in the MatLab command window. Note that
% computing suprathreshold cluster size probabilities can take a long
% time, particularly for low thresholds or large numbers of
% relabellings. Eventually, the Graphics window will come up and the
% results displayed.
%     
% - Results
%=======================================================================
%
% The format of the results page is similar to that of SPM:
%
% A Maximum Intensity Projection (MIP) of the statistic image is shown
% top left: Only voxels significant (corrected) at the chosen level
% \alpha are shown. (If suprathreshold cluster size is being assessed,
% then clusters are shown if they have significant size *or* if they
% contain voxels themselves significant at the voxel level.) The MIP is
% labelled SnPM{t} or SnPM{Pseudo-t}, the latter indicating that
% variance smoothing was carried out.
%
% On the top right a graphical representation of the Design matrix is
% shown, with the contrast illustrated above.
%
% The lower half of the output contains the table of p-values and
% statistics, and the footnote of analysis parameters. As with SPM, the
% MIP is tabulated by clusters of voxels, showing the maximum voxel
% within each cluster, along with at most three other local maxima
% within the cluster. The table has the following columns:
% 
%
% At Cluster Level
%
% * Pcorrected: If "spatial extent" has been assessed, FWE-corrected
%   p-values for region size will be shown.  The corrected P-value for
%   the suprathreshold cluster size is the probability (conditional on
%   the data) of the experiment giving a suprathreshold cluster of size
%   as or more extreme anywhere in the statistic image. This is the
%   proportion of the permutation distribution of the maximal
%   suprathreshold cluster size exceeding (or equalling) the observed
%   size of the current cluster. 
%
% * k: The size (in voxels) of the cluster.
%
%
% At Voxel Level
%
% * PFWE-corr: FWE-corrected non-parametric P-values.  This is the
%   probability of the experiment giving a voxel statistic this extreme
%   anywhere in the statistic image. This is the proportion of the
%   permutation distribution of the maximal statistic exceeding (or
%   equalling) the observed statistic value.
%
% * PFDR-corr: FDR-corrected non-parametric P-values.  This is the
%   smallest FDR level for which the voxel would be significant.  A level
%   q FDR threshold ensures that, on average, no more than q*100% of the
%   suprathreshold voxels will be false positives.
%
% * t / Pseudo-t: The statistic value.
%
% * Puncorrected: uncorrected non-parametric P-values.  This is the
%   result of a non-parametric permutation test completed at this
%   individual voxel.  This P-value is is the proportion of the
%   permutation distribution (at *this* voxel) exceeding (or
%   equalling) the observed statistic value.
%
%
% * {x,y,z} mm: Locations of local maxima.
%
% The SnPM parameters footnote contains the following information:
%
% * Primary threshold: If assessing "spatial extent", the primary
%   threshold used for identification of suprathreshold clusters is
%   printed. If using t-statistics (as opposed to Pseudo-t's), the
%   corresponding upper tail probability is also given.
%
% * Critical STCS: The critical suprathreshold cluster size. This is
%   size above which suprathreshold clusters have significant size at
%   level \alpha It is computed as the 100(1-alpha)%-ile of the
%   permutation distribution of the maximal suprathreshold cluster
%   size. Only shown when assessing "spatial extent".
%
% * alpha: The test level specified.
%
% * Critical threshold: The critical statistic level. This is the value 
%   above which voxels are significant (corrected) at level \alpha.  It
%   is computed as the 100(1-alpha)%-ile of the permutation
%   distribution of the maximal statistic.
%
% * df: The degrees of freedom of the t-statistic. This is printed even
%   if
%   variance smoothing is used, as a guide.
%
% * Volume & voxel dimensions:
%
% * Design: Description of the design
%
% * Perms: Description of the exchangability and permutations used.
%
%
%
% - SnPM files required:
%=======================================================================
% snpm_pp loads parameters and results from the following files, which
% must all be in the same directory:
%       SnPMcfg.mat    - SnPM design configuration
%       SnPM.mat        - SnPM analysis & permutation distribution
%       SnPMt.mat       - Pointlist of (Pseudo) t-statisic for actual labelling
%       XYZ.mat         - Co-ordinates of pointlist
%       SnPM_ST.mat (*) - Suprathreshold cluster statistics (if required)
%
% Further details of the actual variables required from these files are
% given in the main body of snpm_pp
%
% (*) The SnPM_ST.mat file containing the suprathreshold cluster
% information for each of the relabellings can be very large, and is
% only needed if a suprathreshold cluster size test is required. If
% such an analysis is not required, but suprathreshold cluster stats
% were collected, then this file may be deleted, without compromising
% further voxel-level analyses.
%
%_______________________________________________________________________
% Copyright (C) 2013 The University of Warwick
% Id: snpm_pp.m  SnPM13 2013/10/12
% Andrew Holmes, Thomas Nichols, Jun Ding, Darren Gitelman

%-----------------------------functions-called------------------------
% spm_DesMtx
% spm_Tcdf
% spm_clf
% spm_clusters
% spm_figure
% spm_select
% spm_hwrite
% spm_input
% spm_invTcdf
% spm_max
% snpm_mip
% snpm_pp
% spm_str_manip
% spm_t2z
% spm_type
% spm_xyz2e
%-----------------------------functions-called------------------------

%-Variable "decoder" - Following files/variables are required:
%=======================================================================
% NB: Mat files contain additional variables beyond those required here.
%     See function that wrote each file for full definitions.
% SnPM design configuration file: 			SnPMcfg.mat
%-----------------------------------------------------------------------
% - saved by snpm_ui
% H             condition partition of DesMtx for correctly labeled data
% C             covariate partition of DesMtx for correctly labeled data
% B             block     partition of DesMtx for correctly labeled data
% G             confound  partition of DesMtx for correctly labeled data
% HCBGnames     string matrix of column names of [H C B G]
% PiCond        Permuted conditions matrix, one labelling per row, actual
%               labelling on first row
% sPiCond       String describing permutations in PiCond
% bhPerms       Flag indicating use of "half permutations" trick
% CONT          Contrast (only one)
% bVarSm        Flag for variance smoothing (Pseudo t-statistics)
% sVarSm        Sring describing variance Smoothing (empty if bVarSm=0)
% bST           Flag for collection of superthreshold info 
% sDesign       Description of PlugIn design
% 
% SnPM analysis & permutation distribution file:	SnPM.mat
%-----------------------------------------------------------------------
% - saved by snpm_cp
% S		- Volume analyzed (in voxels)
% V             Image file handles (see spm_vol)
% df            Residual degrees of freedom of raw t-statistic
% MaxT          2xnPerm matrix of [max;min] t-statistics per perm
% ST_Ut         Threshold above which suprathreshold info was collected.
%               Voxel locations, t and perm are saved in SnPM_ST.mat for
%               t's greater than ST_Ut. ST_Ut=Inf if not saving STCdata
%
% Pointlist of (Pseudo) t-statisic for actual labelling:SnPMt.mat
%-----------------------------------------------------------------------
% - saved by snpm_cp
% SnPMt		- 1xS matrix of voxel statistics (t or pseudo-t)
%
% Co-ordinates of pointlist:				XYZ.mat
%-----------------------------------------------------------------------
% - saved by snpm_cp
% XYZ		- 3xS matrix of co-ordinates [x;y;z] of voxels on SnPMt
%
% Suprathreshold cluster statistics:			SnPM_ST.mat
%-----------------------------------------------------------------------
% - saved by snpm_cp
% SnPM_ST	- Suprathreshold cluster statistics, see snpm_cp.m
% NB: This file is only required for suprathreshold cluster size analysis



%-Setup
%=======================================================================
global defaults
if isempty(defaults), spm_defaults; end
global SnPMdefs
if isempty(SnPMdefs), snpm_defaults; end
MLver=version;MLver=MLver(1);

fprintf('\nSnPM: snpm_pp\n'),fprintf('%c','='*ones(1,72)),fprintf('\n')

%-Initialise variables & constants
%-----------------------------------------------------------------------
tol = 1e-4;	% Tolerance for comparing real numbers
		% Two reals with abs(a-b)<tol are considered equal
		% ( Reals have to be compared for equality when        )
		% ( computing FWE-corrected p-values                   )

Dis = SnPMdefs.Results_distmin;  % mm distance between sub regions
Num = SnPMdefs.Results_nbmax;    % Maximum number of sub-regions
units = {'mm' 'mm' 'mm'};

%-SetUp figure window
%-----------------------------------------------------------------------
Finter = spm_figure('FindWin','Interactive');
Fgraph = spm_figure('FindWin','Graphics');
if isempty(Fgraph), Fgraph=spm_figure('Create','Graphics'); end
spm_clf(Finter), spm_clf(Fgraph)
set(Finter,'Name','SnPM PostProcess');


%-Get Data
%=======================================================================
% Get analysis directory
if nargin==0
  tmp  = spm_select(1,'SnPM.mat','Select SnPM.mat for analysis...');
  CWD  = spm_str_manip(tmp,'hd');
end
if nargin>1
  job=varargin{1};
  BATCH=true;
else
  BATCH=false;
end

%-Skip reports in BATCH; will create them later
if BATCH
  Report = {job.Report};
else
  % Only can specify multiple reports outside of BATCH mode; otherwise the results will fly by,
  % one after another without pausing
  Report = {'FWEreport','FDRreport','MIPtable'};
end


%-Load Config file & SnPM permutation data
load(fullfile(CWD,'SnPMcfg'))
load(fullfile(CWD,'SnPM'))
load(fullfile(CWD,'SnPMucp'))


%-Ask whether positive or negative effects be analysed
%-----------------------------------------------------------------------
if BATCH
  if STAT == 'T'
    bNeg = job.Tsign==-1;
  else
    bNeg = 0;
  end
else
  if STAT == 'T'
    bNeg = spm_input('Positive or negative effects?',1,'b','+ve|-ve',[0,1],1);
  else
    bNeg = 0;  
    str = 'Positive effects';
    spm_input('F-statistic, so only:','+1','b',str,1);
  end
end

%-Take MaxT for increases or decreases according to bNeg
MaxT = MaxT(:,bNeg+1);
nPerm = size(PiCond,1);  %nPerm is consistent with the one in snpm_cp
nPermReal = size(MaxT,1); %different with nPerm when bhPerms==1
[StMaxT, iStMaxT] = sort(MaxT);

%-Load statistic image
%-----------------------------------------------------------------------
load(fullfile(CWD,'SnPMt'))
load(fullfile(CWD,'XYZ'))
XYZ0=XYZ;

%-Negate if looking at negative contrast
%-----------------------------------------------------------------------
if bNeg
	SnPMt    = -SnPMt;
	CONT     = -CONT;
end

%-Get ORIGIN, etc
DIM    = [V(1).dim(1) V(1).dim(2) V(1).dim(3)];
M=V(1).mat(1:3, 1:3);
VOX=sqrt(diag(M'*M))';
MAT    = V(1).mat;
IMAT   = inv(MAT);
ORIGIN = IMAT(1:3,4);

% Template vol structure
Vs0 = V(1);

% Vs0 = struct('fname',	'',...
% 	     'dim',	[DIM,spm_type('float')],...
% 	     'mat',	MAT,...
% 	     'pinfo',	[1 0 0]',...
% 	     'descrip',	'');

% Process Nonparmaetric P-values
if bNeg 
  % Here, nPermReal has already been doubled if bhPerms=1
  SnPMucp = 1+1/nPermReal-SnPMucp;
end  
sSnPMucp = sort(SnPMucp);



%-Write out filtered statistic image?  (Get's done later)
%-----------------------------------------------------------------------
if BATCH
  if isfield(job.WriteFiltImg,'WF_no')
    WrtFlt=0;
  else
    WrtFlt=1;
    WrtFltFn=job.WriteFiltImg.name;
    
    if isempty(spm_str_manip(WrtFltFn, 'e'))
        WrtFltFn = [WrtFltFn '.nii'];
    end
  end
else
  WrtFlt = spm_input('Write filtered statistic img?','+1','y/n',[1,0],2);
  if WrtFlt
    WrtFltFn = 'SnPMt_filtered';
    WrtFltFn=spm_input('Filename ?','+1','s',WrtFltFn);
    WrtFltFn = [WrtFltFn, '.img'];
  end
end


%-Get inference parameters
%=======================================================================

% Map of options from Batch System
% [Root]     > .Thr
% Voxel-Level    > .Vox
% . Significance  > .VoxSig
% . . Uncorrected Nonparametric P | Uncorrected T or F | FDR Corrected | FWE Corrected 
%    > .Pth                            .TFth                .FDRth          .FWEth
% Cluster-Level  > .Clus
% . Cluster size statistic > .ClusSize
% . . Cluster-Forming Threshold > .CFth
% . . Significance Level > .ClusSig
% . . . Uncorrected k | FWE Corrected 
%   >   .Cth           .FWEthC

%-Get corrected threshold
%-----------------------------------------------------------------------
u         = NaN;   % Statistic Image threshold
alpha_ucp = NaN;   % Uncorrected P-value image threshold
C_MaxT    = NaN;   % Statistic image threshold set by alph_FWE
C_STCS    = NaN;   % Cluster size threshold (set directly by uncorrected
		   % threshold or by alph_FWE)
		   
alph_FWE  = NaN;   % FWE rate of a specified u threshold
alph_FDR  = NaN;   % FDR rate of a specified alpha_ucp
  
if BATCH
    bSpatEx = isfield(job.Thr,'Clus');
    if ~bSpatEx
        % Voxel-wise inference
        tmp = fieldnames(job.Thr.Vox.VoxSig);
        switch tmp{1}
        case 'Pth'
            alpha_ucp = BoundCheck(job.Thr.Vox.VoxSig.Pth,[0 1],'Invalid Uncorrected P');
            alph_FDR  = snpm_P_FDR(alpha_ucp,[],'P',[],sSnPMucp');
        case 'TFth'
            u         = BoundCheck(job.Thr.Vox.VoxSig.TFth,[0 Inf],'Negative Threshold!');
            alph_FWE  = sum(MaxT > u -tol) / nPermReal;
        case 'FDRth'
            alph_FDR  = BoundCheck(job.Thr.Vox.VoxSig.FDRth,[0 1],'Invalid FDR level');
            alpha_ucp = snpm_uc_FDR(alph_FDR,[],'P',[],sSnPMucp');
        case 'FWEth'
            alph_FWE  = BoundCheck(job.Thr.Vox.VoxSig.FWEth,[0 1],'Invalid FWE level');
            iFWE      = ceil((1-alph_FWE)*nPermReal);
            if alph_FWE<1
                C_MaxT=StMaxT(iFWE);
            else
                C_MaxT = 0;
            end
            u = C_MaxT;
        end
    else
        % Cluster-wise inference
        if exist(fullfile(CWD,'SnPM_ST.mat'))~=2 & exist(fullfile(CWD,'STCS.mat'))~=2
            error(['SnPM:NoClusterInfo', 'ERROR: Cluster-wise inference requested, but no cluster information saved.\n',...
            'Re-configure analysis changing "Cluster inference" to "Yes" and re-run.\n'])
        end
        %%% Sort out the cluster-forming threshold
        if pU_ST_Ut==-1  % No threshold was set in snpm_ui.
            if isnan(job.Thr.Clus.ClusSize.CFth)
                error('SnPM:NoClusterFormingThresh', 'ERROR: Cluster-forming threshold set to NaN in results with "slow" cluster inference method used in compoutation.  \nRe-run results specifying a cluster-forming threshold.\n')
            end
            % Save original ST_Ut
            ST_Ut_0 = ST_Ut;
            CFth=job.Thr.Clus.ClusSize.CFth;
            if (CFth<=0)
                error('SnPM:InvalidClusterFormingThresh', 'ERROR: Cluster-forming threshold must be strictly positive.\nRe-run results with a cluster-forming threshold greater than 0.\n')
            end
            if bVarSm
                %-If using pseudo-statistics then can't use (uncorrected) 
                % upper tail p-values to specify primary threshold
                pCFth = NaN;
                if (CFth<1)
                    warning('snpm_pp:pseudoTFormingThresholdP',...
                        ['Pseudo-T cluster-forming threshold defined by '...
                        'P-value using Gaussian approximation P=' num2str(pU_ST_Ut)...
                        ' -> Z=' num2str(ST_Ut) '; actual Pseudo-T threshold '...
                        'unknown but may be higher than ' num2str(ST_Ut) '.']);
                    pCFth = CFth;
                    CFth = spm_invNcdf(1-CFth);
                    
%                     error(sprintf('ERROR: Cluster-forming threshold specified as a P-value (%g), but uncorrected P-values are unavailable for the pseudo t (smoothed variance t-test).  \nRe-run results with a cluster-forming threshold greater than 1.\n',ST_Ut))
                end
                
                if (CFth < ST_Ut)%(CFth>=ST_Ut-tol)
                    if isnan(pCFth)
                        error('SnPM:InvalidClusterFormingThresh', sprintf('ERROR: Cluster-forming threshold of %0.2f specified, but statistic image information only saved for %0.2f and greater. \nRe-run results with a cluster-forming threshold of %0.2f or higher.  (Alternatively, increase SnPMdefs.STprop in snpm_defaults.m, re-start SnPM, and re-compute analysis.)\n',CFth,ST_Ut,ST_Ut))
                    else
                        error('SnPM:InvalidClusterFormingThresh', sprintf('ERROR: Cluster-forming threshold of P=%0.4f (T=%0.2f) specified, but statistic image information only saved for %0.2f and greater. \nRe-run results with a cluster-forming P-value threshold of %0.2f or lower.  (Alternatively, increase SnPMdefs.STalpha in snpm_defaults.m, re-start SnPM, and re-compute analysis.)\n',pCFth,CFth,ST_Ut,p_ST_Ut))
                    end
                end
            else
                %-Statistic image is t with df degrees of freedom
                p_ST_Ut  = STalpha;
                if (CFth < 1)
                    pCFth = CFth;
                    CFth = spm_invTcdf(1-CFth,df);
                else
                    pCFth = NaN;
                    if (abs(CFth-ST_Ut)<=tol)
                    CFth=ST_Ut; % If tmp is very close to ST_Ut, set tmp equal to ST_Ut.
                    end
                end

                if (CFth < ST_Ut) %(CFth>=ST_Ut-tol)
                    if isnan(pCFth) % statistic-value cluster-forming threshold
                        error('SnPM:InvalidClusterFormingThresh', sprintf('ERROR: Cluster-forming threshold of %0.2f specified, but statistic image information only saved for %0.2f and greater. \nRe-run results with a cluster-forming threshold of %0.2f or higher.  (Alternatively, increase SnPMdefs.STalpha in snpm_defaults.m, re-start SnPM, and re-compute analysis.)\n',CFth,ST_Ut,ST_Ut))
                    else
                        error('SnPM:InvalidClusterFormingThresh', sprintf('ERROR: Cluster-forming threshold of P=%0.4f (T=%0.2f) specified, but statistic image information only saved for %0.2f and greater. \nRe-run results with a cluster-forming P-value threshold of %0.2f or lower.  (Alternatively, increase SnPMdefs.STalpha in snpm_defaults.m, re-start SnPM, and re-compute analysis.)\n',pCFth,CFth,ST_Ut,p_ST_Ut))
                    end
                end
            end
            if (abs(CFth-ST_Ut)<=tol)
                CFth = ST_Ut; % If tmp is very close to ST_Ut, set tmp equal to ST_Ut.
            end
                ST_Ut = CFth;
        else % Threshold *was* set in snpm_ui.
            if ~isnan(job.Thr.Clus.ClusSize.CFth)
                error('SnPM:InvalidClusterFormingThresh', sprintf('ERROR: Cluster-forming threshold of T=%0.2f was already set during analysis configuration; hence, in results, cluster-forming threshold must be left as "NaN".\nRe-run results with cluster-forming threshold set to NaN.\n',ST_Ut))
            end
        end
        u=ST_Ut; % Flag use of a statistic-value threshold
        % Inference details...
        tmp = fieldnames(job.Thr.Clus.ClusSize.ClusSig);
        switch tmp{1}
            case 'Cth'
                C_STCS = job.Thr.Clus.ClusSize.ClusSig.Cth;
            case 'PthC'
                alpha_ucp = BoundCheck(job.Thr.Clus.ClusSize.ClusSig.PthC,[0 1],'Invalid uncorrected P(k)');
            case 'FWEthC'
                alph_FWE  = BoundCheck(job.Thr.Clus.ClusSize.ClusSig.FWEthC,[0 1],'Invalid FWE level (cluster-level inference)');
                iFWE      = ceil((1-alph_FWE)*nPermReal);
        end
    end % END: Cluster-wise inference

else  % GUI, interative inference specification

  str_img =[STAT,'|P'];
  switch spm_input('Results for which img?','+1','b',str_img,[],1)
    
   case 'P' % Use the P-image
    bSpatEx = 0; % Cluster-wise inference won't be performed anyway. 
    str = 'FDR|None';
    switch spm_input('Use corrected threshold?','+1','b',str,[],1)
      
     case 'FDR' % False discovery rate
      %---------------------------------------------------------------	
      alph_FDR = spm_input('FDR-Corrected p value threshold','+0','r',0.05,1,[0,1]);
      alpha_ucp = snpm_uc_FDR(alph_FDR,[],'P',[],sSnPMucp');
		
     otherwise  %-Uncorrected: no adjustment
      %%%% Now ask: Threshold statistic image or Uncorr P-value Image ?
      %%%% If stats image, maybe do ST, no FDR-level of thresh;
      %%%% If P-value image, no ST, can find FDR-level of thresh

      % p for conjunctions is p of the conjunction SPM
      %---------------------------------------------------------------
      alpha_ucp = spm_input('Uncorrected p value threshold','+0','r',0.01,1,[0,1]);
      alph_FDR = snpm_P_FDR(alpha_ucp,[],'P',[],sSnPMucp');
    end 
  
   case STAT % Use the T-image
    %-Ask whether SupraThreshold cluster size test required
    %----------------------------------------------------------------------- 
    %-To have cluster size inference, need
    %  1. Spatial extent information was collected (bST=1),
    %  2. SnPM_ST.mat or STCS.mat file exists
    bSpatEx = bST & (exist(fullfile(CWD,'SnPM_ST.mat'))==2|exist(fullfile(CWD,'STCS.mat'))==2);
    
    if bSpatEx
      str = 'Voxelwise|Clusterwise';
      bSpatEx = spm_input('Inference method?','+1','b',str,[1 2],1)==2;
    else
      str = 'Voxelwise';
      spm_input('Inference method:','+1','b',str,1);
    end
    
    if ~bSpatEx % Voxel-wise inference
      str = 'FWE|None';
      switch spm_input('Voxelwise: Use corrected thresh?','+1','b',str,[],1)
	
       case 'FWE' % family-wise false positive rate
	%---------------------------------------------------------------
	alph_FWE  = spm_input('FWE-Corrected p value threshold','+0','r',0.05,1,[0,1]);
	iFWE=ceil((1-alph_FWE)*nPermReal);
	if alph_FWE<1
	  C_MaxT=StMaxT(iFWE);
	else
	  C_MaxT = 0;
	end
	u = C_MaxT;
	
       otherwise  %-NB: no adjustment
	%%%% Now ask: Threshold statistic image or Uncorr P-value Image ?
	%%%% If stats image, maybe do ST, no FDR-level of thresh;
	%%%% If P-value image, no ST, can find FDR-level of thresh

	% p for conjunctions is p of the conjunction SPM
   	%---------------------------------------------------------------
	if bVarSm, str = 'pseudo t'; else, str = sprintf('t_{%d}',df); end
	u  = spm_input(['threshold (',str,')'],'+0','r',0.01,1,[0,Inf]);
	alph_FWE = sum(MaxT > u -tol) / nPermReal;
      end
    
    else % Cluster-wise inference
    
      if pU_ST_Ut==-1  % No threshold was set in snpm_ui.
        %-Get primary threshold for STC analysis if requested
	%-----------------------------------------------------------------------
	% Save original ST_Ut
	ST_Ut_0 = ST_Ut;
	%-Threshold must be greater or equal to that (ST_Ut) used to collect
	% suprathreshold data in snpm_cp
	%-If a test level alpha has been set, then it there's no sense in having
	% the threshold greater than C_MaxT, above which voxels are individually 
	% significant
	tmp = 0;
	if bVarSm
	  %-If using pseudo-statistics then can't use (uncorrected) 
	  % upper tail p-values to specify primary threshold
	  while ~(tmp>=ST_Ut-tol)
	    tmp = spm_input(sprintf(...
		'Clus-def thresh(pseudo t>%4.2f)',ST_Ut),'+0');
	  end
	  if (abs(tmp-ST_Ut)<=tol)
	  tmp=ST_Ut; % If tmp is very close to ST_Ut, set tmp equal to ST_Ut.
	  end
	else
	  %-Statistic image is t with df degrees of freedom
	  p_ST_Ut  = STalpha;
	  while ~( tmp>=ST_Ut-tol | (tmp>0 & tmp<=p_ST_Ut))
	    tmp = spm_input(sprintf(...
		'Clus-def thresh(p<=%4.2fIt>=%4.2f)',p_ST_Ut,ST_Ut),'+0','r',ST_Ut,1); 
	  end
	  clear p_ST_Ut
	  if (tmp < 1)
	    tmp = spm_invTcdf(1-tmp,df);
	  else 
	    if (abs(tmp-ST_Ut)<=tol)
	      tmp=ST_Ut; % If tmp is very close to ST_Ut, set tmp equal to ST_Ut.
	    end
	  end
	end
	ST_Ut = tmp;
      end
      u=ST_Ut; % Flag use of a statistic-value threshold

      str = 'FWE|Uncorr';
      switch spm_input('Clusterwise: Use corrected thresh?','+1','b',str,[],1)
       case 'FWE' % family-wise false positive rate
        %---------------------------------------------------------------
	alph_FWE  = spm_input('FWE-Corrected p value threshold','+0','r',0.05,1,[0,1]);
	iFWE=ceil((1-alph_FWE)*nPermReal);
	
       case 'Uncorr' %Uncorrected cluster size threshold
	
	str = 'ClusterSize|P-value';
	switch spm_input('Define uncorrected','+1','b',str,[],1)
	  
	 case 'ClusterSize'
	  C_STCS = spm_input('Uncorr cluster size threshold','+0','w',0,1);
	  
	 case 'P-value'
	  alpha_ucp = spm_input('Uncorrected p value threshold','+0','r',0.01,1,[0,1]);
	  
	end
      end
    end	
  end
end



% Workflow for statistical Inference:
%
% Results for T or P
%
% 1) P
%
%   case 'FDR': alph_FDR is set by spm_input; alpha_ucp is calculated from
%               alph_FDR and is used as the threshold.
%   case 'None': alpha_ucp is set by spm_input; alph_FDR is calculated from
%                alpha_ucp. alpha_ucp is used as the threshold.
% 
% 2) T
%
%   a) Voxel-wise
%      case 'FWE': alph_FWE is set by spm_input; u=C_MaxT; iFWE is
%                  calculated. u is used as the threshold.
%      case 'None': u is set by spm_input; alph_FWE is calculated from u.
%                   u is used as the threshold. 
%
%   b) Cluster-wise
%      If cluster defining threshold not set, ask for pU_ST_Ut
%      case 'FWE': alph_FWE is set by spm_input; iFWE is calculated from
%                  alph_FWE. C_STCS is calculated from iFWE and is used
%                  as the threshold. 
%      case 'Uncorr': 
%               i) 'ClusterSize', C_STCS is defined directly and used as  
%                   the threshold; 
%               ii) 'P-value', alpha_ucp is set by spm_input. C_STCS is
%                   calculated from alpha_ucp and used as the threshold. 
%



%=======================================================================
%- C O M P U T A T I O N
%=======================================================================
set(Finter,'Pointer','Watch')

%-Calculate distribution of Maximum Suprathreshold Cluster size
%-Calculate critical Suprathreshold Cluster Size
%=======================================================================
if bSpatEx
	fprintf('Working on spatial extent...\n');
	%-Compute suprathreshold voxels - check there are some
	%---------------------------------------------------------------
	fprintf('\tComputing suprathreshold voxels...');
	Q     = find(SnPMt > ST_Ut);
	SnPMt = SnPMt(Q);
	XYZ   = XYZ(:,Q);
	
	if isempty(Q)
	  set(Finter,'Pointer','Arrow')
	  figure(Fgraph)
	  axis off
	  text(0,0.97,CWD,'Fontsize',16,'FontWeight','Bold');
	  str=sprintf('No voxels above threshold %4.2f\n',ST_Ut);
	  text(0,0.93,str);
	  fprintf(['WARNING: ' str])
	  if length(strmatch('FWEreport',Report))>0
	    ShowDist(MaxT,C_MaxT,alph_FWE,[],[],[],'max');
	    if ~BATCH
	      if spm_input('Review permutation distributions.',1,'bd',...
			   'Print & Continue|Continue',[1,0],1)
		spm_print
	      end
	      spm_clf(Fgraph)
	    end
	  end
	  if length(strmatch('FDRreport',Report))>0
	    axis off
	    text(0,0.97,'Uncorrected P Permutation Distributions','Fontsize',16,'FontWeight',...
		 'Bold');
	    ShowDist(SnPMucp,alpha_ucp,alph_FDR,[],[],[],'uncor');
	    if ~BATCH
	      if spm_input('Review permutation distributions.',1,'bd',...
			   'Print|Done',[1,0],1)
		spm_print
	      end
	    end
	  end
	  return
	end
	fprintf('done\n')

	%-Load & condition statistics
	%---------------------------------------------------------------
        if pU_ST_Ut==-1 % No threshold was set in snpm_ui.
	  fprintf('\tLoading & conditioning SupraThreshold statistics...');
      try
        load(fullfile(CWD,'SnPM_ST'))
      catch exception
          if strcmp(exception.identifier, 'MATLAB:load:unableToReadMatFile')
              warning('SnPM:SnPMSTFileNotLOaded', ...
                  ['SnPM_ST file can not be loaded. Consider using' ...
                  ' ''set cluster-forming threshold now (fast)'' option' ...
                  ' in SnPM ''Specify''.']);
          end
          % Rethrow exception           
          throw(exception)
      end
	  %-SnPM_ST stores columns of [x;y;z;abs(t);perm] with perm negative
	  % where the exceedence was t < -ST_Ut_0
	  %-Trim statistics according to threshold ST_Ut, if ST_Ut > ST_Ut_0
	  tmp = find(SnPM_ST(4,:)>ST_Ut);
	  SnPM_ST = SnPM_ST(:,tmp);
	  clear tmp;	    
	  
	  SnPM_ST_Pos = SnPM_ST(:,SnPM_ST(5,:)>0);
	  SnPM_ST_Neg = SnPM_ST(:,SnPM_ST(5,:)<0);
	  SnPM_ST_Neg(5,:) = -SnPM_ST_Neg(5,:);
	  
	  fprintf('done\n')
	  
	  %-Calculate distribution of Maximum SupraThreshold Cluster size
	  %---------------------------------------------------------------
	  fprintf('\tComputing dist. of max SupraThreshold cluster size: ');
	  STCS = snpm_STcalc('init',nPerm); 
	  fprintf('\nPerms left:     ');
	  for i = nPerm:-1:1
	    if (rem(i,10)==0)
	      fprintf('\b\b\b\b%-4u',i)
	      drawnow
	    end
	    
	    if STAT== 'F'
	      loop = 1;
	    else
	      loop = 1:2;
	    end
	    
        for isPos= loop  %1 for positive; 2 for negative
            if isPos==1
                SnPM_ST = SnPM_ST_Pos;
            else
                SnPM_ST = SnPM_ST_Neg;
            end
            tQ = (SnPM_ST(5,:)==i);
            if any(tQ)
                %-Compute cluster labellings for this perm
                Locs_mm = SnPM_ST(1:3,tQ);
                Locs_mm (4,:) = 1;
                Locs_vox = IMAT * Locs_mm;
                
                % Sometimes Locs_vox are not exactly integers and this raises an
                % error later in the code. Here check that the values are
                % integers with respect to a level of absolute tolerance (~10^-14)
                % and enforce Locs_vox to be integers.
                % (As in snpm_cp)                
                diffWithRounded = max(abs(Locs_vox(:)-round(Locs_vox(:))));
                tolerance = 10^-10;
                if diffWithRounded > tolerance
                 Locs_vox_alter = MAT\Locs_mm;
                 diffWithRounded_alter = max(abs(Locs_vox_alter(:)-round(Locs_vox(:))));
                 error('SnPM:NonIntegerLocsvox', ['''Locs_vox'' must be integers (difference is ' num2str(diffWithRounded) ...
                     ' or ' num2str(diffWithRounded_alter) ')']);
                else
                 Locs_vox = round(Locs_vox); 
                end

                STCS = snpm_STcalc('update',STCS, SnPM_ST(4,tQ),...
                   Locs_vox(1:3,:),isPos,i,ST_Ut,df);
            end
            if i==1
                %-Save perm 1 stats for use later - [X;Y;Z;T;perm;STCno]
                tmp = spm_clusters(Locs_vox(1:3,:));
                if isPos==1
                    STCstats_Pos = [ SnPM_ST(:,tQ); tmp];
                    if bNeg==0
                        STCstats=STCstats_Pos;
                    end
                else
                    STCstats_Neg = [ SnPM_ST(:,tQ); tmp];
                    if bNeg==1
                        STCstats=STCstats_Neg;
                    end
                end
            end
        end
    end
	  fprintf('\b\b\b\bdone\n');
	  
	  if bhPerms   %Double the STCS variables.
	    STCS = snpm_STcalc('double',STCS);
	  end
	  
	  %-Get the stats from STCS structure that will be used later
	  STCS_MxK = STCS.MxK(:,bNeg+1);
	  STCS_K=cat(1,STCS.K{:,bNeg+1});	   
	  
       else % A threshold was set in snpm_ui.
	  %-Load & condition statistics
	  %---------------------------------------------------------------
	  fprintf('\tLoading SupraThreshold statistics...');
	  load(fullfile(CWD,'STCS'))
	  STCS_MxK = STCS.MxK(:,bNeg+1);
	  STCS_K=cat(1,STCS.K{:,bNeg+1});
	  if (bNeg==0)
	    load(fullfile(CWD,'SnPM_pp'))
	  else
	    load(fullfile(CWD,'SnPM_pp_Neg'))
	    STCstats = STCstats_Neg;
	  end       
	end	
        
	%-Compute critical SupraThreshold Cluster size
	if isnan(C_STCS)
	  
	  if ~isnan(alph_FWE)
	    
	    % STCS_MxK: a vector of maximum cluster sizes of all the permutations.
	    % STCS_MxK = STCS.MxK(:,bNeg+1);  
	    [StMaxSTCS, iStMaxSTCS] = sort(STCS_MxK);
	    if alph_FWE < 1
		C_STCS = StMaxSTCS(iFWE);
	    else
		C_STCS = 0;
	    end
	  
	  elseif ~isnan(alpha_ucp)
	    
	    % STCS_K: a vector of all the cluster sizes of all the
            % permutations.
	    
	    [StKSTCS, iStKSTCS] = sort(STCS_K);
	    if alpha_ucp < 1
	        iucp=ceil((1-alpha_ucp)*length(STCS_K));
		C_STCS = StKSTCS(iucp);
	    else
	        C_STCS = 0;
	    end
	  
	  end
	
	end
	
	%-Check XYZ for points > ST_Ut in perm 1 matches
	% XYZ computed above for SnPMt > ST_Ut
	%if pU_ST_Ut==-1 
	% if ~all(all( SnPM_ST(1:3,SnPM_ST(5,:)==1) == XYZ ))
	%	error('SnPM:InvalidSTXYZ', 'ST XYZ don''t match between STCS & thresh')
	% end
	%else
	 if ~all(all( STCstats(1:3,:) == XYZ ))
		error('SnPM:InvalidSTXYZ', 'ST XYZ don''t match between STCS & thresh')
	 end
	%end
end

%-Save some time consuming results
%-----------------------------------------------------------------------
if bSpatEx & pU_ST_Ut==-1
  save SnPM_pp STCstats_Pos
  if STAT == 'T'
     save SnPM_pp_Neg STCstats_Neg
  end
  save STCS STCS
end


%-Filter data at specified corrected p-value alpha
%=======================================================================
if bSpatEx
    %-Analysing spatial extent
    % if alph_FWE is set then only do it when alph_FWE<1
    % otherwise(alph_FWE=NaN), do the analysis.
    if (alph_FWE<1|isnan(alph_FWE))
	%-Filter on significance of cluster size
	%---------------------------------------------------------------
	fprintf('Filtering statistic image, clusterwise...');
	nSTC     = max(STCstats(6,:));
	STCS_K1  = diff(find([diff([0,sort(STCstats(6,:))]),1]));
	Q        = [];
	for i = 1:nSTC
		tQ = find(STCstats(6,:)==i);
		% Note, we discard a cluster even if it has a significant
                % voxel in it
		if ( STCS_K1(i) > C_STCS )
			Q        = [Q tQ];
		end
	end
	if ~isempty(Q)
		SnPMt    = SnPMt(Q);
		XYZ      = XYZ(:,Q);
		STCstats = STCstats(:,Q);
	end
	fprintf('done\n')
    end
else
	if ~isnan(u) 
	  % Thresholding based on t-values
	  Q =  find(SnPMt > u);
	elseif ~isnan(alpha_ucp)
	  % Thresholding based on uncorrected P-values
	  Q = find(SnPMucp < alpha_ucp);
	else
	  error('SnPM:CodingError', 'Coding error')
	end
	if ~isnan(u)
	  fprintf('Filtering statistic image, voxelwise...');
	else
	  fprintf('Filtering uncorr. p image, voxelwise...');
	end
	if length(Q)
		SnPMt = SnPMt(Q);
		XYZ   = XYZ(:,Q);
	end
	fprintf('done\n')
end



%-Return if there are no voxels
%-----------------------------------------------------------------------
if isempty(Q)
	set(Finter,'Pointer','Arrow')
	figure(Fgraph)
	axis off
	text(0,0.97,CWD,'Fontsize',16,'FontWeight','Bold');
	tmp='voxels'; if bSpatEx, tmp='suprathreshold clusters'; end
	str='';
	if ~isnan(u)
	  if bSpatEx
	    str=sprintf(...
		'No %s significant at k>=%d (alpha=%6.4f FWE-corrected)\n',...
		tmp,C_STCS,alph_FWE);
	  else
	    str=sprintf(...
		'No %s significant at u>=%0.2f (alpha=%6.4f FWE-corrected)\n',...
		tmp,u,alph_FWE);
	  end
	else
	  str=sprintf(...
	      'No %s significant at alpha=%6.4f (alpha=%6.4f FDR-corrected)\n',...
	      tmp,alpha_ucp,alph_FDR);
	end
	if ~isempty(str)
	  text(0,0.93,str);
	  fprintf(['WARNING: ' str])
	end
	if length(strmatch('FWEreport',Report))>0
	  if bSpatEx,
	    ShowDist(MaxT,C_MaxT,alph_FWE,STCS_MxK,C_STCS,alph_FWE,'max');
	  else	   
	    ShowDist(MaxT,C_MaxT,alph_FWE,[],[],[],'max');
	  end
	  if ~BATCH
	    if spm_input('Review permutation distributions.',1,'bd',...
			 'Print & Continue|Continue',[1,0],1)
	      spm_print
	    end
	  end
	end
	if length(strmatch('FDRreport',Report))>0
	  spm_clf(Fgraph)
	  axis off
	  text(0,0.97,'Uncorrected P Permutation Distributions','Fontsize',16,'FontWeight',...
	       'Bold');
	  ShowDist(SnPMucp,alpha_ucp,alph_FDR,[],[],[],'uncor');
	  if ~BATCH
	    if spm_input('Review permutation distributions.',1,'bd',...
			 'Print|Done',[1,0],1)
	      spm_print
	    end  
	  end
	end
	return
end

%-Characterize local excursions in terms of maxima:
% #voxels STC_N; MaxTs STC_SnPMt; locations STC_XYZ, & region# STC_r
%-----------------------------------------------------------------------
%===== SnPM99 change =============================================
TempXYZmm = XYZ;
TempXYZmm(4,:) = 1;
TempXYZvoxel = IMAT*TempXYZmm;
TempXYZvoxel= TempXYZvoxel(1:3,:);

[STC_N, STC_SnPMt, STC_XYZ, STC_r] = spm_max(SnPMt,TempXYZvoxel);

TempXYZvoxel = STC_XYZ;
TempXYZvoxel(4,:) = 1;
TempXYZmm = MAT * TempXYZvoxel;
STC_XYZ = TempXYZmm(1:3,:);
%===== SnPM99 change =============================================

%-Compute corrected significances for local maxima, & regions (if required)
%-----------------------------------------------------------------------
Pt = ones(size(STC_r));
for i = 1:length(STC_r)
	%-Use a > b -tol rather than a >= b to avoid comparing reals
	Pt(i) = sum(MaxT > STC_SnPMt(i) -tol) / nPermReal;
end

%Since we have already considered both bNeg =0 and 1 situations in
%ucp, the Pu has the same formula for bNeg=0 and 1.

        ucp = zeros(1,prod(DIM));
	ucp(spm_xyz2e(XYZ0,V)) = SnPMucp;
        Pu  = ucp(spm_xyz2e(STC_XYZ,V)) ;
        Pfdr = snpm_P_FDR(Pu,[],'P',[],sSnPMucp');
if bSpatEx
	%-Compute FWE-corrected p-values and uncorrected p-values for region size: pSTSC_SS
	Pn    = ones(size(STC_r));     % FWE-corrected p-values
        Pun   = ones(size(STC_r));     % Uncorrected p-values
	for i = 1:length(STC_r)
	  Pn(i) = sum(STCS_MxK>=STC_N(i)) / nPermReal;
          Pun(i) = sum(STCS_K>=STC_N(i)) / length(STCS_K);	  
	end
end

% Display only if *not* in command line mode
if ~spm_get_defaults('cmdline')
    
%=======================================================================
%-D I S P L A Y :   Max report
%=======================================================================

if length(strmatch('FWEreport',Report))>0

  spm_clf(Fgraph)
  figure(Fgraph)
  axis off
  if (bSpatEx)
    text(0,0.97,'Permutation Distribution','Fontsize',16,'FontWeight', ...
	 'Bold');
    
    ShowDist(MaxT,C_MaxT,alph_FWE,STCS_MxK,C_STCS,alph_FWE,'max');
  else	     
    text(0,0.97,'Permutation Distributions','Fontsize',16,'FontWeight','Bold');
    ShowDist(MaxT,C_MaxT,alph_FWE,[],[],[],'max');
  end
  
  if ~BATCH
    if spm_input('Review permutation distributions.',1,'bd',...
		 'Print & Continue|Continue',[1,0],1)
      spm_print
    end
  end
end


%=======================================================================
%-D I S P L A Y :   FDR (uncorrected P-values) report
%=======================================================================

if length(strmatch('FDRreport',Report))>0

  spm_clf(Fgraph)
  figure(Fgraph)
  axis off
  text(0,0.97,'Uncorrected P Permutation Distributions','Fontsize',16,'FontWeight',...
       'Bold');
  ShowDist(SnPMucp,alpha_ucp,alph_FDR,[],[],[],'uncor');

  if ~BATCH
    if spm_input('Review permutation distributions.',1,'bd',...
                  'Print & Continue|Continue',[1,0],1)
      spm_print
    end
  end

end



%=======================================================================
%-D I S P L A Y :   Maximium intenisty projection of SPM{Z}
%=======================================================================

if length(strmatch('MIPtable',Report))>0

  spm_clf(Fgraph)
  figure(Fgraph)
  axis off
  
  hmip = axes('Position',[0.05 0.5 0.5 0.5]);
  snpm_mip(SnPMt,XYZ,MAT,DIM); axis image
  if bVarSm
    TitlStr='SnPM{Pseudo-t}';
  else
    TitlStr=sprintf('SnPM{%s}',STAT);
  end
  title(TitlStr,'FontSize',16,'Fontweight','Bold')

  
  %-Design matrix and contrast
  %=======================================================================
  hDesMtx = axes('Position',[0.65 0.6 0.2 0.2]);
  image((spm_DesMtx('Sca', [H,C,B,G],HCBGnames) + 1)*32)
  xlabel 'Design Matrix'
  set(hDesMtx,'XTick',[],'XTickLabel','')
  nPar   = size([H,C,B,G],2);
  hConAxes = axes('Position',[0.65 0.81 0.2 0.1]);
  if STAT == 'F'
    imagesc(CONT,[-1 1]);
    set(gca,'Tag','ConGrphAx',...
	    'Box','on','TickDir','out',...
	    'XTick',spm_DesRep('ScanTick',nPar,10),'XTickLabel','',...
	    'XLim',	[0,nPar]+0.5,...
	    'YTick',[1:size(CONT,1)],....
	    'YTickLabel','',...
	    'YLim',	[0,size(CONT,1)]+0.5	)
  else
    h = bar(CONT(1,:), 'FaceColor',[1 1 1]*.8, 'BarWidth', 1);  
    tX = get(h,'XData'); tY = get(h,'YData');
    bar_width = get(h, 'BarWidth');
    set(gca,'Xlim',[min(tX(:))-bar_width/2 max(tX(:))+bar_width/2]) 
    axis off
  end
  title 'contrast'; 
  

  %-Table of regional effects
  %=======================================================================
  %-Table headings
  %-----------------------------------------------------------------------
  hTable = axes('Position',[0.1 0.1 0.8 0.46],...
		'YLim',[0,27],'YLimMode','manual',...
		'DefaultTextInterpreter','Tex',...
		'DefaultTextVerticalAlignment','Baseline',...
		'Visible','off');
  %	      'DefaultHorizontalAlignment','right',...
  y = 26;
  dy=1;
  text(0,y,['P values & statistics:   ',spm_str_manip(CWD,'a40')],...
       'FontSize',12,'FontWeight','Bold','Interpreter','none');
  y  = y -dy;
  line([0 1],[y y],'LineWidth',3,'Color','r')
  y = y -dy;
  
  tCol       = [  0.00      0.12    0.26             ...	%-Cluster
		  0.34      0.46    0.63      0.72   ...  %-Voxel
		  0.88     0.94    1.00];        	%-XYZ
  
  PF    = spm_platform('fonts');   %-Font names (for this platform)


  %-Construct table header
  %-----------------------------------------------------------------------
  set(gca,'DefaultTextFontName',PF.helvetica,'DefaultTextFontSize',10)
  
  Hp = [];
  h  = text(0.10,y,	'cluster-level','FontSize',10,'HorizontalAlignment','Center');		
  h  = line([tCol(1),0.30],[1,1]*(y-dy/4),'LineWidth',0.5,'Color','r');	
  h  = text(tCol(1),y-9*dy/8,	'\itp_{FWE-corr}');        Hp = [Hp,h];
  h  = text(tCol(2),y-9*dy/8,	'\itp_{uncorr}');      Hp = [Hp,h];
  h  = text(tCol(3),y-9*dy/8,	'\itk ');			
  
  text(0.50,y,		'voxel-level','FontSize',10,'HorizontalAlignment','Center');
  line([tCol(4),0.80],[1,1]*(y-dy/4),'LineWidth',0.5,'Color','r');
  h  = text(tCol(4),y-9*dy/8,	'\itp_{FWE-corr}');	
  h  = text(tCol(5),y-9*dy/8,	'\itp_{FDR-corr}');	
  
  if ~bVarSm
    h  = text(tCol(6),y-9*dy/8,	sprintf('\\it%s',STAT));
  else 
    h  = text(tCol(6)-0.02,y-9*dy/8,	'Pseudo-t');
  end 
  
  h  = text(tCol(7),y-9*dy/8,	'\itp_{uncorr}');
  
  text(tCol(8),y-dy/2,'{x,y,z} mm','FontSize',10);
  
  
  y     = y - 7*dy/4;
  line([0 1],[y y],'LineWidth',1,'Color','r')
  y     = y - 5*dy/4;
  
  Fmtst = {	'%0.4f', '%0.4f', '%0.0f', ...                  %-Cluster
		'%0.4f', '%0.4f', '%6.2f','%0.4f', ...		%-Voxel
		'%3.0f', '%3.0f', '%3.0f'};			%-XYZ
  
  %-Column Locations
  %-----------------------------------------------------------------------
  %tCol       = [  0.07      0.30  ...			%-Cluster
  %	           0.50      0.62      0.77           ...  %-Voxel
  %                0.86 0.93 1.00];			%-XYZ
  
  
  %-Table Headers
  %----------------------------------------------------------------------
  TabDat.tit = sprintf('%s: p-values adjusted for search volume',TitlStr);
  
  TabDat.hdr = {...
      'cluster',  'p(FWE-corr)',  '\itp\rm_{FWE-corr}';...
      'cluster',  'p(uncorr)',    '\itp\rm_{uncorr}';...
      'cluster',  'k',            '\itk';...
      'peak',     'p(FWE-corr)',  '\itp\rm_{FWE-corr}';...
      'peak',     'p(FDR-corr)',  '\itq\rm_{FDR-corr}';...
      'peak',      STAT,          sprintf('\\it%s',STAT);...
      'peak',     'p(uncorr)',    '\itp\rm_{uncorr}';...
      '',         'x,y,z {mm}',   [units{:}]}';...

  TabDat.fmt = Fmtst;
  %-Table filtering note
  %----------------------------------------------------------------------
  if isinf(Num)
    TabDat.str = sprintf('table shows all local maxima more than %.1fmm apart',Dis);
  else
    TabDat.str = sprintf('table shows %d local maxima more than %.1fmm apart',Num,Dis);
    end 
  TabDat.dat = cell(0,10);

  StrAttr = {'Fontsize',10,'ButtonDownFcn','get(gcbo, ''UserData'')',...
	     'HorizontalAlignment','right'};
  StrAttrB = {StrAttr{:},'FontWeight','Bold'};
  if ~bSpatEx
      set(Hp,'Visible','off')
  end

  %-List of maxima
  %-----------------------------------------------------------------------
  r = 0;
  bUsed = zeros(size(STC_SnPMt));
  while max(STC_SnPMt.*(~bUsed)) & (y > 3)
    
    [null, i] = max(STC_SnPMt.*(~bUsed));	% Largest t value
    j         = find(STC_r == STC_r(i));	% Maxima in same region
    r         = r + 1;				% Next row
    
    %-Print region and largest maximum
    %-------------------------------------------------------------------
    
    %	text(0.00,y,sprintf('%0.0f',r),'UserData',r,StrAttrB{:})
    if bSpatEx
      TabDat.dat(r,1:2) = {Pn(i),Pun(i)};
      if (y>3) 
	text(tCol(1)+0.09,y,sprintf(Fmtst{1},Pn(i)),'UserData',Pn(i), ...
	     StrAttrB{:})
	text(tCol(2)+0.09,y,sprintf(Fmtst{2},Pun(i)),'UserData',Pun(i), ...
	     StrAttrB{:})
      end
    end
    
    TabDat.dat(r,3:10)={STC_N(i),Pt(i),Pfdr(i),STC_SnPMt(i),Pu(i),STC_XYZ(1,i),STC_XYZ(2,i),STC_XYZ(3,i)};

    if (y>3)
      text(tCol(3)+0.04,y,sprintf(Fmtst{3},STC_N(i)),'UserData',STC_N(i),StrAttrB{:})
      text(tCol(4)+0.08,y,sprintf(Fmtst{4},Pt(i)),'UserData',Pt(i),StrAttrB{:})
      text(tCol(5)+0.09,y,sprintf(Fmtst{5},Pfdr(i)),'UserData',Pfdr(i),StrAttrB{:})
      text(tCol(6)+0.04,y,sprintf(Fmtst{6},STC_SnPMt(i)),'UserData',STC_SnPMt(i),StrAttrB{:})
      
      text(tCol(7)+0.09,y,sprintf(Fmtst{7},Pu(i)),'UserData',Pu(i),StrAttr{:})
      text(tCol(8),y,sprintf(Fmtst{8},STC_XYZ(1,i)),'UserData',STC_XYZ(:,i),StrAttrB{:})
      text(tCol(9),y,sprintf(Fmtst{9},STC_XYZ(2,i)),'UserData',STC_XYZ(:,i),StrAttrB{:})
      text(tCol(10),y,sprintf(Fmtst{10},STC_XYZ(3,i)),'UserData',STC_XYZ(:,i),StrAttrB{:})
      y = y -1;
    end
    
    %-Print up to Num secondary maxima (> Dis apart)
    %-------------------------------------------------------------------
    [null, k] = sort(-STC_SnPMt(j));	% Sort on t value
    D         = i;
    for i = 1:length(k)
      d     = j(k(i));
      if min( sqrt( sum((STC_XYZ(:,D) - ...
			 STC_XYZ(:,d)*ones(1,size(D,2))).^2) ) ) > Dis;
	if length(D) < Num
	  r = r + 1;				% Next row

	  TabDat.dat(r,4:10)={Pt(i),Pfdr(i),STC_SnPMt(i),Pu(i),STC_XYZ(1,i),STC_XYZ(2,i),STC_XYZ(3,i)};

	  if (y>3)
	    text(tCol(4)+0.08,y,sprintf(Fmtst{4}, Pt(d)),'UserData',Pt(d),StrAttr{:})
	    text(tCol(5)+0.09,y,sprintf(Fmtst{5}, Pfdr(d)),'UserData',Pfdr(d),StrAttr{:})
	    text(tCol(6)+0.04,y,sprintf(Fmtst{6}, STC_SnPMt(d)), 'UserData',STC_SnPMt(d),StrAttr{:})
	    
	    text(tCol(7)+0.09,y,sprintf(Fmtst{7}, Pu(d)),'UserData',Pu(d),StrAttr{:})
	    text(tCol(8),y,sprintf(Fmtst{8}, STC_XYZ(1,d)),'UserData',STC_XYZ(:,d),StrAttr{:})
	    text(tCol(9),y,sprintf(Fmtst{9}, STC_XYZ(2,d)),'UserData',STC_XYZ(:,d),StrAttr{:})
	    text(tCol(10),y,sprintf(Fmtst{10}, STC_XYZ(3,d)),'UserData',STC_XYZ(:,d),StrAttr{:})
	  
	    y = y -1;
	  end
	  D = [D d];

	end
      end
    end
    
    bUsed(j) = (bUsed(j) | 1 );		%-Mark maxima as "used"
  end
  clear i j k D d r
  
  
  %-Footnote with SnPM parameters
  %=======================================================================
  TabDat.ftr    = cell(0,2);
  line([0,1],[0.5,0.5],'LineWidth',1,'Color','r')
  r = 1;
  y = 0;
  if bSpatEx
    if ~bVarSm
      TabDat.ftr{r,1}='Cluster-defining thresh. = %7.4f (p = %6.4f)';  
      TabDat.ftr{r,2}=[ST_Ut,spm_Tcdf(-ST_Ut,df)];
    else
      TabDat.ftr{r,1}='Cluster-defining thresh. = %7.4f';              
      TabDat.ftr{r,2}=ST_Ut;
    end
    text(0,y,sprintf(TabDat.ftr{r,1},TabDat.ftr{r,2}),'FontSize',8)
    r = r + 1;

    TabDat.ftr{r,1}='Critical STCS = %d voxels';    TabDat.ftr{r,2}=C_STCS;
    text(0.7,y,sprintf(TabDat.ftr{r,1},TabDat.ftr{r,2}),'FontSize',8)
    y = y -0.8; r = r + 1;

    if ~isnan(alph_FWE)
      TabDat.ftr{r,1}='Cluster threshold: FWE-corr. P value= %6.4f';     
      TabDat.ftr{r,2}=alph_FWE;
    elseif ~isnan(alpha_ucp)
      TabDat.ftr{r,1}='Cluster threshold: Uncorr. P value= %6.4f';       
      TabDat.ftr{r,2}=alpha_ucp;
    else
      TabDat.ftr{r,1}='Cluster threshold: Uncorr. cluster size STCS= %d';
      TabDat.ftr{r,2}=C_STCS;
    end
    text(0,y,sprintf(TabDat.ftr{r,1},TabDat.ftr{r,2}), 'FontSize',8)
    r = r + 1;
  else
    %make the format similar as spm.
    %text(0,y,sprintf('alpha = %6.4f, df = %d',alpha,df),'FontSize',8)
    if ~isnan(u)
      TabDat.ftr{r,1}='Height threshold: statistic u= %6.2f (%0.4f FWE)';  
      TabDat.ftr{r,2}=[u,alph_FWE];
    else
      TabDat.ftr{r,1}='Height threshold: Nonparam. P value alpha= %0.4f (%0.4f FDR)';  
      TabDat.ftr{r,2}=[alpha_ucp, alph_FDR];
     end
    text(0,y,sprintf(TabDat.ftr{r,1},TabDat.ftr{r,2}),'FontSize',8)
    r = r + 1;
  end
  
  if bVarSm
    TabDat.ftr{r,1}='Degrees of freedom = n/a (variance smoothing)';  
    TabDat.ftr{r,2}=[];
  else
    TabDat.ftr{r,1}='Degrees of freedom = [%d %d]';  
    TabDat.ftr{r,2}=[df1,df];
  end
  text(0.7,y,sprintf(TabDat.ftr{r,1},TabDat.ftr{r,2}),'FontSize',8)
  y=y-0.8;  r = r + 1;

  TabDat.ftr{r,1}='';  
  TabDat.ftr{r,2}=[];

  TabDat.ftr{r,1}='Design: %s';  
  TabDat.ftr{r,2}=sDesign;
  text(0,y,sprintf(TabDat.ftr{r,1},TabDat.ftr{r,2}),'FontSize',8);
  y=y-0.8;  r = r + 1;

  TabDat.ftr{r,1}='Search vol: %d cmm, %d voxels';  
  TabDat.ftr{r,2}=[S*abs(prod(VOX)),S];
  text(0,y,sprintf(TabDat.ftr{r,1},TabDat.ftr{r,2}), 'FontSize',8)
  y=y-0.8;  r = r + 1;

  TabDat.ftr{r,1}='Voxel size: [%5.2f, %5.2f, %5.2f] mm';  
  TabDat.ftr{r,2}=abs(VOX);
  text(0.7,y,sprintf(TabDat.ftr{r,1},TabDat.ftr{r,2}), 'FontSize', 8)
  r = r + 1;
  
  TabDat.ftr{r,1}='Perms: %s';  
  TabDat.ftr{r,2}=sPiCond;
  text(0,y,sprintf(TabDat.ftr{r,1},TabDat.ftr{r,2}),'FontSize',8);
  r = r + 1;

  if bVarSm
    TabDat.ftr{r,1}=sVarSm;  
    TabDat.ftr{r,2}=[];
    text(0,y,sprintf(TabDat.ftr{r,1},TabDat.ftr{r,2}),'FontSize',8)
    y = y -0.8;  r = r + 1;
  end
  
  save TabDat TabDat

  if ~BATCH
    %Set a button, so the user can decide whether to print the page of
    %results to spm2.ps.
    if spm_input('Review results.',1,'bd','Print|Done',[1,0],1)
      spm_print
    end 
  end
  
  set(Finter,'Pointer','Arrow')

end

end

%- Image output?
%=======================================================================
%-Write out filtered SnPMt?
if WrtFlt
	
  Fname = WrtFltFn;
  
  %-Set name string
  %---------------------------------------------------------------
  if ~bVarSm
    tmp = sprintf('SPMt - %d df',df);
  else
    tmp = 'SnPMt - pseudo t';
  end
  
  %-Reconstruct filtered image from XYZ & SnPMt
  %---------------------------------------------------------------
  t = zeros(1,prod(DIM))*NaN;
  t(spm_xyz2e(XYZ,V)) = SnPMt;
  if ~bSpatEx
    if ~isnan(u)
      tmp=sprintf('%s p<%10g FWE-corr @ voxel level',tmp,alph_FWE);
    else
      tmp=sprintf('%s p<%10g uncorr @ voxel level',tmp,alpha_ucp);
    end
  else
    if ~isnan(alph_FWE)
      tmp=sprintf('%s p<%10g corrected @ cluster level, u=%4.2f',...
		  tmp,alph_FWE,ST_Ut);
    elseif ~isnan(alpha_ucp)
      tmp=sprintf('%s p<%10g uncorrected @ cluster level, u=%4.2f',...
		  tmp,alpha_ucp,ST_Ut);
    else
      tmp=sprintf('%s cluster size >%d @ cluster level, u=%4.2f',...
		  tmp,C_STCS,ST_Ut);
    end    
  end
  
  %-Write out to image file
  %---------------------------------------------------------------
  Vs = snpm_clone_vol(Vs0, Fname, tmp); 
  Vs =  spm_create_vol(Vs);
  t = reshape(t,DIM);
  for p=1:Vs.dim(3)
    Vs = spm_write_plane(Vs,t(:,:,p),p);
  end
  Vs =  sf_close_vol(Vs);
  clear t
end

%-Reset Interactive Window
%-----------------------------------------------------------------------
spm_figure('Clear','Interactive')



function ShowDist(T,cT,aT,C,cC,aC,Typ)
% Shows permutation distribution
% FORMAT ShowDist(T,cT,aT,C,cC,aC,Typ)
%
% T    - Voxelwise statistic permutation dist
% cT   - Voxelwise statistic threshold
% aT   - (FWE or FDR) size of cT
% C    - Cluster size statistic permutation dist
% cC   - Cluster size statistic threshold
% aC   - (FWE) size of cC
% Typ  - 'uncor' for uncorrected values or 'max' for corrected
% 
% Display permutation distributions on current figure, using position
% pos.
%
% We assume that the observed (aka correctly labeled data) are the first
% element in the distribution.
%

if nargin<4, C   =[]; end
if nargin<5, cC  =[]; end
if nargin<6, aC  =[]; end
if nargin<7, Typ ='max'; end

figure(spm_figure('FindWin','Graphics'));

if strcmp(Typ,'max')
  TitStr = 'Permutation Distribution:  Maximum Statistic';
elseif  strcmp(Typ,'uncor')
  TitStr = 'Distribution: Observed Uncorrected P-values';
else
  error('SnPM:UnknownString', 'Unknown type string')
end  

pos1 = [0.125 0.50 0.75 0.3];
pos2 = [0.125 0.08 0.75 0.3];

%
% Display intensity perm dist
%
% Bin width rule from Scott, "Multivariate Density Estimation", 1992, pg 55.
%
axes('position',pos1)
BinWd  = 3.5*std(T)*length(T)^(-1/3);  
nBin   = floor((max(T)-min(T))/BinWd)+1;
nBin   = min(max(nBin,10),50);

hist(T,nBin);
set(get(gca,'Children'),'FaceColor',[.5 .5 .5]) 
title([ TitStr],'FontSize',14)
Ylim = get(gca,'ylim'); Xlim = get(gca,'xlim');
if strcmp(Typ,'max')
  str = sprintf('Observed Maximum Statistic = %g   ',T(1));
  if ~isnan(cT)
    str = {str, ...
	   sprintf('Voxelwise statistic threshold = %g (%0.4f FWE)',cT,aT)};
    text(cT+diff(Xlim)*0.01,Ylim(2)*0.85,'Threshold','FontSize',10);
  end
  xlabel(str)
  line(T(1)*[1 1],Ylim.*[1 0.95],'LineStyle','--','color',[1 0 0]);
  text(T(1)+diff(Xlim)*0.01,Ylim(2)*0.95,'Observed','FontSize',10)
  if (Xlim(1)<=cT) & (cT<=Xlim(2))
  line(cT*[1 1],Ylim.*[1 0.85],'LineStyle','-');
  end
else
  str = sprintf('Smallest Observed Uncorr. P-value = %g    ',min(T));
  if ~isnan(cT)
    str = {str, ...
	   sprintf('Voxelwise P-value threshold = %g (%0.4f FDR)',cT,aT)};
    text(cT+diff(Xlim)*0.01,Ylim(2)*0.85,...
	 sprintf('Threshold (%g)',cT),'FontSize',10);
  end
  xlabel(str);
  snpm_abline('h',length(T)/nBin);
end


if strcmp(Typ,'uncor')

  %
  % Show FDR plot
  %

  S = length(T);
  Ts = sort(T);
  figure(spm_figure('FindWin','Graphics'));
  axes('position',pos2)
  loglog((1:S)/S,Ts,'r-o')
  set(get(gca,'children'),'LineWidth',1,'MarkerSize',2)
  snpm_abline(0,1);
  snpm_abline(0,aT,'LineStyle','-','color','red')
  if ~isnan(cT) & (cT~=0)
    snpm_abline('h',cT,'LineStyle','-','Color','blue');
    text(10^(log10(1/S)*0.9),10^(log10(cT)*.85),...
	 sprintf('%3.3g FDR: p=%g',aT,cT),'color','red','FontSize',10)
  end
  title('FDR illustration: loglog pp-plot','FontSize',20)
  axis image
  xlabel('index/(number of voxels)'); 
  ylabel('Ordered uncorrected nonparametric p-value')


elseif ~isempty(C)

  %
  % Display cluster size perm dist
  %

  axes('position',pos2)
  rC = C.^(1/3);
  BinWd  = 3.5*std(rC)*length(rC)^(-1/3);
  nBin   = floor((max(rC)-min(rC))/BinWd)+1;
  nBin   = min(max(nBin,10),50);

  hist(rC,nBin);
  set(get(gca,'Children'),'FaceColor',[.5 .5 .5]) 
  if strcmp(Typ,'max')
    title('Permutation Distribution: Maximum Cluster Size (cuberoot plot)','FontSize',14)
    xlabel(sprintf('Max Cluster Size Observed = %g  Cluster Threshold = %g',C(1),cC))
  else
    warning('SnPM:UncorrectedClustSize', 'Uncorrected cluster size not supported yet!')
    title('Permutation Distribution: (Uncorrected) Cuberoot Cluster Size','FontSize',14)
    xlabel(sprintf('Max Cluster Size Observed = %g  Cluster Threshold = %g',C(1),cC))
  end
  Ylim = get(gca,'ylim'); Xlim = get(gca,'xlim');
  set(gca,'Xticklabel',num2str(str2double(cellstr(get(gca,'Xticklabel'))).^3))
  line(rC(1)*[1 1],Ylim.*[1 0.95],'LineStyle','--','color',[1 0 0]);
  text(rC(1)+diff(Xlim)*0.01,Ylim(2)*0.95,'Observed','FontSize',10)
  line(cC.^(1/3)*[1 1],Ylim.*[1 0.85],'LineStyle','-');
  text(cC.^(1/3)+diff(Xlim)*0.01,Ylim(2)*0.85,'Threshold','FontSize',10)

end

return



function Vs = sf_close_vol(Vs)
switch spm('ver')
 case 'SPM99'
  % n/a; File closed by spm_write_plane
 case 'SPM2'
  % Vs = spm_close_vol(Vs);
  % Don't need to close images in SPM5??
end         
return

function v = BoundCheck(val,Range,ErrMsg)
if val<Range(1) | val>Range(2)
  error('SnPM:RangeError', [ErrMsg ': ' num2str(val)])
else
  v=val;
end
return