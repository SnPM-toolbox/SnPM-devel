function job = snpm_auto(SPM,varargin)
% FUNCTION job = snpm_auto(SPM,[dir,nPerm,vFWHM,bVol])
%
% Converts existing parametric SPM analysis' SPM.mat into a job structure

if nargin<2; dir=varargin{1};   else dir=''; end
if nargin<3; nPerm=varargin{2}; else nPerm=[]; end
if nargin<4; vFWHM=varargin{3}; else vFWHM=[]; end
if nargin<5; bVol=varargin{4};  else bVol=[]; end

snpm_defaults

load(SPM)
% Set up 'default' job structure, with as much from SPM.mat as possible
job = JobDefault(SPM,dir,nPerm,vFWHM,bVol);
  

%
% Decide what type of plug it is!
% 

% Is it one-sample t-test?
if size(SPM.xX.X,2)==1 && all(SPM.xX.X==SPM.xX.X(1)) % Yes!
  job.DesignName = 'MultiSub: One Sample T test on differences; 1 condition';
  job.DesignFile = 'snpm_pi_OneSampT';
  
% Is it a two-sample t-test
% elseif ????  % 

% Is it a simple correlation
% elseif ????  % 


return



function job = JobDefault(SPM,dir,nPerm,vFWHM,bVol)
global SnPMdefs


cov = struct('c',[],...
	     'cname','');
masking = struct('tm',struct('tm_none',1),...
		 'im',1,...
		 'em',{''});
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
  job.dir = dir;
end
if ~isempty(nPerm)
  job.nPerm = nPerm;
end
if ~isempthy(vFWHM);
  job.vFWHM = vFWHM;
end
if ~isempthy(bVol)
  job.bVol = bVol;
end

job.P = SPM.xY.P;

% set up threshold/masking/global 


%%%% STILL TO DO BY TOM
