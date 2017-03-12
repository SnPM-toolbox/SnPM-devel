function matlabbatch = snpm_auto(SPM,varargin)
% Converts existing parametric SPM analysis' SPM.mat into a job structure
% FORMAT job = snpm_auto(SPM,[dir,nPerm,vFWHM,bVol])
%_______________________________________________________________________
% Arguments
%     SPM    - SPM structure or filepath to SPM.mat, to be base of new SnPM run
%     dir    - Directory where new SnPM to be run
%     nPerm  - Number of permutations to use
%     vFWHM  - Variance FWHM smoothing (a 3 vector; defaults to [0 0 0])
%     bVol   - Use 3D (1) or 2D (0) processing
%

if nargin>=2; dir=varargin{1};   else dir=''; end
if nargin>=3; nPerm=varargin{2}; else nPerm=[]; end
if nargin>=4; vFWHM=varargin{3}; else vFWHM=[]; end
if nargin>=5; bVol=varargin{4};  else bVol=[]; end

snpm_defaults

if ~isstruct(SPM)
  load(SPM)
end
% Set up 'default' job structure, with as much from SPM.mat as possible
job = JobDefault(SPM,dir,nPerm,vFWHM,bVol);
  

%
% Decide what type of plug it is!
% 

X  = SPM.xX.X;
Design = SPM.xsDes.Design;
iid = IsIID(SPM);
ErrMsg = '';

if IsLevel1(SPM) 

  ErrMsg = 'First level fMRI time series model detected';

% One-sample t-test
% snpm_pi_OneSampT.m
elseif strcmp(Design,'One sample t-test')

  if ~iid
    ErrMsg = 'One-sample model detected, but has dependence or het. var.';
  else
    job.DesignName = 'MultiSub: One Sample T test on differences; 1 condition';
    job.DesignFile = 'snpm_pi_OneSampT';
  end

% Two-sample t-test
% snpm_pi_TwoSampT.m
elseif strcmp(Design,'Two-sample t-test')

  if ~iid
    ErrMsg='Two-sample model detected, but has dependence or het. var.';
  else
    job.DesignName = '2 Groups: Two Sample T test; 1 scan per subject';
    job.DesignFile = 'snpm_pi_TwoSampT';
  end
  
% Paired t-test
% snpm_pi_PairT.m
elseif strcmp(Design,'Paired t-test')

  % Should always be iid

  job.DesignName = 'MultiSub: Paired T test; 2 conditions, 1 scan per condition';
  job.DesignFile = 'snpm_pi_PairT';
  
% Is it a between-subject ANOVA?
% snpm_pi_ANOVAbtwnS.m
elseif strcmp(Design,'ANOVA')

  if ~iid
    ErrMsg='ANOVA (between subjects) model detected, but has dependence or het. var.';
  else
    job.DesignName = 'Between Subject ANOVA, k groups, 1 scan per subject';
    job.DesignFile = 'snpm_pi_ANOVAbtwnS';
  end

% Is it a within-subject ANOVA?
% snpm_pi_ANOVAwithinS.m
elseif strcmp(Design,'ANOVA - within subject')

  % Works even with dependence
  job.DesignName = 'Within Subject ANOVA, k diffs/contrasts per subject';
  job.DesignFile = 'snpm_pi_ANOVAwithinS';

% Is it a simple correlation
% snpm_pi_Corr.m
elseif (size(X,2)==2 && sum(all(diff(X)==0))==1) 

  if ~iid
    ErrMsg='Simple correlation model detected, but has dependence or het. var.';
  else
    job.DesignName = 'MultiSub: Simple Regression (correlation); single covariate of interest, 1 scan per subject';
    job.DesignFile = 'snpm_pi_Corr';
  end

% % Simple correlation extracted from a multiple regression
% % snpm_pi_Corr.m
% elseif strcmp(Design,'Multiple regression')
% 
%   % If you can find any singleton T contrasts that are *not* on the
%   % intercept... go for it
% 
%   if ~iid
%     ErrMsg='Simple correlation (simple regression) model detected, but has dependence or het. var.';
%   else
%     job.DesignName = '';
%     job.DesignFile = '';
%   end

end

if ~isempty(job.DesignName) 
  fprintf(['Design detected: ' job.DesignName]);
elseif ~isempty(ErrMsg)
  error(ErrMsg)
else
  error(['Design not compatible (' Design ')']);
end

% Don't worry about these:
% snpm_pi_TwoSampPairT.m  --  Within subject ANOVA (2 obs) + between subject factor (2 levels)
% snpm_pi_TwoSampTss.m    --  Single subject two-sample T, for PET (has EBs)
% snpm_pi_Corr1S.m        --  Single subject correlation, for PET (has EBs)
% snpm_pi_PairTrand.m     --  Paired T-test, randomisation instead of permutation

end

matlabbatch{1}.spm.tools.snpm.des.OneSampT = job;

return

function job = JobDefault(SPM,dir,nPerm,vFWHM,bVol)
global SnPMdefs

cov = struct('c',[],...
	     'cname','');
masking = struct('tm',struct('tm_none',1),...
		 'im',1,...
		 'em',{{''}});
globalm = struct('gmsca',struct('gmsca_no',1),...
		 'glonorm',1);
job.DesignName = '';
job.DesignFile = '';
job.dir = {};
job.P = {};
job.cov = repmat(cov,0,0);
job.nPerm = SnPMdefs.nPerm;
job.vFWHM = [0 0 0];
job.bVolm = SnPMdefs.bVolm;
job.ST = struct('ST_none',0);
job.masking = masking;
job.globalc = struct('g_omit',1);
job.globalm = globalm;

% Set stuff from SPM.mat -- or options set by user

if isempty(dir)
  job.dir = {fullfile(SPM.swd,'SnPM')};
else
  job.dir = {dir};
end
if ~isempty(nPerm)
  job.nPerm = nPerm;
end
if ~isempty(vFWHM);
  job.vFWHM = vFWHM;
end
if ~isempty(bVol)
  job.bVol = bVol;
end

job.P = SPM.xY.P;

return

function Test = IsIID(SPM)

if ~isfield(SPM,'xVi')
  Test = true
elseif ~isfield(SPM.xVi,'Vi')
  Test = true;
elseif length(SPM.xVi.Vi)==1 && all(all(full(SPM.xVi.Vi{1})==eye(size(SPM.xVi.Vi{1}))))
  Test = true;
else
  Test = false;
end

return

function Test = IsLevel1(SPM)

Test = isfield(SPM.xY,'RT');

return
