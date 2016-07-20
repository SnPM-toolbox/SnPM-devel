% Mfile snpm_pi_TwoSampTss
% SnPM PlugIn design module - SingleSubj simple activation, 2cond
% SingleSub: Two Sample T test; 2 conditions, replications
% FORMAT snpm_pi_TwoSampTss 
%
% See body of snpm_ui for definition of PlugIn interface.
%_______________________________________________________________________
%
% snpm_pi_TwoSampTss is a PlugIn for the SnPM design set-up program,
% creating design and permutation matrix appropriate for single
% subject, two condition activation (with replication) experiments.
%
%-Number of permutations
%=======================================================================
% There are nScan-choose-nRepl possible permutations, where 
% nScan is the number of scans and nRepl is the number of replications
% (nScan = 2*nRepl).  Matlab doesn't have a choose function but
% you can use this expression
%
%	prod(1:nScan)/prod(1:nRepl)^2
% 
%
%-Prompts
%=======================================================================
% '# replications per condition':  Here you specify how many times
% each of the two conditions were repeated.
%
% 'Size of exchangability block':  This is the number of adjacent
% scans that you believe to be exchangeable. The most common cause of 
% non-exchangeability is a temporal confound. Exchangibility blocks
% are required to be the same size, a divisor of the number of scans.
% See snpm.man for more information and references.
%
% 'Select scans in time order':  Enter the scans to be analyzed. 
% It is important to input the scans in time order so that temporal
% effects can be accounted for.
%
% 'Enter conditions index: (B/A)':  Using A's to indicate activation
% scans and B's to indicate baseline, enter a sequence of 2*nRepl
% letters, nRepl A's and nRepl B's, where nRepl is the number of 
% replications. Spaces are permitted.
%
%
%-Variable "decoder" - This PlugIn supplies the following:
%=======================================================================
% - core -
% P             - string matrix of Filenames corresponding to observations
% iGloNorm      - Global normalisation code, or allowable codes
%               - Names of columns of design matrix subpartitions
% PiCond        - Permuted conditions matrix, one labelling per row, actual
%                 labelling on first row
% sPiCond       - String describing permutations in PiCond
% sHCform       - String for computation of HC design matrix partitions
%                 permutations indexed by perm in snpm_cp
% CONT          - single contrast for examination, a row vector
% sDesign       - String defining the design
% sDesSave      - String of PlugIn variables to save to cfg file
%
% - design -
% H,Hnames      - Condition partition of design matrix, & effect names
% B,Bnames      - Block partition (constant term), & effect names
%
% - extra -
% iCond         - Condition indicator vector
% iRepl         - Replication indicator vector
% Xblk          - Size of exchangability block
%
%_______________________________________________________________________
% Copyright (C) 2013 The University of Warwick
% Id: snpm_pi_TwoSampTss.m  SnPM13 2013/10/12
% Thomas Nichols & Andrew Holmes, Camille Maumet

%-----------------------------functions-called------------------------
% spm_DesMtx
% spm_select
% spm_input
%-----------------------------functions-called------------------------

%-Initialisation
%-----------------------------------------------------------------------
nCond    = 2;			% Number of conditions
iGloNorm = '123';		% Allowable Global norm. codes
sDesSave = 'iCond iRepl Xblk';	% PlugIn variables to save in cfg file

%-Get number of replications per condition - 2 x nRepl design
nRepl    = job.Tss_repc;%spm_input('# replications per condition','+1');

global SnPMdefs
if SnPMdefs.shuffle_seed
    % Shuffle seed of random number generator
    try
        rng('shuffle');
    catch
        % Old syntax        
        rand('seed',sum(100*clock));
    end
end

%-Work out exchangability blocks - Assumme Xblks of equal size
%-----------------------------------------------------------------------
%-Valid sizes of X-blk (Xblk) are integer divisors of nCond*nRepl that
% are multiples of nCond
tmp      = nRepl./[nRepl:-1:1];
Xblk     = nCond * tmp(floor(tmp)==ceil(tmp));
tmp      = int2str(Xblk(1));
for i=2:length(Xblk), tmp=str2mat(tmp,int2str(Xblk(i))); end
Xblk     = job.TwosampTss_Block; %spm_input('Size of exchangability block','+1','b',tmp,Xblk);
nXblk    = nCond*nRepl/Xblk;
iXblk    = meshgrid(1:nXblk,1:Xblk); iXblk = iXblk(:)';


%-Compute permutations of conditions
%=======================================================================
%-Compute permutations for a single exchangability block
%-----------------------------------------------------------------------
%-Generate all labellings of Xblk scans as +/- 1
% XblkPiCond=[];
% for i=0:Xblk-1
% 	XblkPiCond=[ones(2^i,1),XblkPiCond;-ones(2^i,1),XblkPiCond];
% end
% %-Trim to labellings with balance of conditions
% XblkPiCond=XblkPiCond(sum(XblkPiCond')==0,:);
nOfPerm = nchoosek(Xblk,Xblk/2);
XblkPiCond = -ones(nOfPerm, Xblk);
% Label affected to group label "-1"
alternativeGroup = nchoosek(1:Xblk,Xblk/2);
XblkPiCond(sub2ind(size(XblkPiCond), repmat(1:nOfPerm, Xblk/2,1)', alternativeGroup)) = 1;

%-Now build up the complete set of possibe labellings: PiConds
%-----------------------------------------------------------------------
nXblkPiCond=size(XblkPiCond,1);
PiCond=XblkPiCond;
for i=2:nXblk
	tmp=size(PiCond,1);
	[iSup iInf] = meshgrid(1:nXblkPiCond,1:tmp);
	iSup = iSup(:)';
	iInf = iInf(:)';
	PiCond=[XblkPiCond(iSup,:),PiCond(iInf,:)];
end
%-Check PiConds sum to zero within Xblks
if ~all(all(PiCond*spm_DesMtx(iXblk)==0))
	error('SnPM:InvalidPiCond', 'Invalid PiCond computed!'), end

%-Take half of PiCond. This is always possible.
%-----------------------------------------------------------------------
% Here, PiCond should *always* satisfy:
% all(all(PiCond(PiCond(:,1)==1,:)==flipud(-PiCond(PiCond(:,1)==-1,:))))
PiCond=PiCond(PiCond(:,1)==1,:);
bhPerms=1;

%-Get filenames and iCond, the condition labels
%=======================================================================
nScan = nCond*nRepl;
P = strvcat(job.P);% spm_select(nScan,'image','Select scans in time order');
perm=[]; while(isempty(perm))
    %tmp=['Enter conditions index: (B/A) [',int2str(nCond*nRepl),']'];
    %iCond = spm_input(tmp,'+1','s');
    %-Convert A/B notation to +/-1 vector - assume A-B is of interest
    %iCond = abs(upper(iCond(~isspace(iCond))));
    iCond = job.condidx;
    iCond = iCond-min(iCond); iCond = -iCond/max([1,iCond])*2+1;    
    %-Check validity of iCond
    % All valid iConds for this design are in the PiCond matrix
    % (or possibly in -PiCond if bhPerms)
    if length(iCond)==nCond*nRepl
	perm = find(all((meshgrid(iCond,1:size(PiCond,1))==PiCond)'));
	if (bhPerms), perm=[perm,...
	    -find(all((meshgrid(iCond,1:size(PiCond,1))==-PiCond)'))]; end
	if length(perm)>1, error('SnPM:InvalidPiCond', 'Multiple iCond in PiCond'), end
	if isempty(perm), error('SnPM:InvalidiCond', 'Invalid iCond for this design'), end
    else
	error('SnPM:InvalidIndices', ['Enter indicies for ',int2str(nCond*nRepl),' scans'])
    end
end

%-Build iRepl
%-----------------------------------------------------------------------
iRepl=cumsum(iCond==-1).*(iCond==-1) + cumsum(iCond==1).*(iCond==1);


%-Shuffle PiCond so actual iCond is first, negate if necc. (on bhPerms)
%-----------------------------------------------------------------------
if (perm<0), PiCond=-PiCond; perm=-perm; end
%-Lame last ditch check just to make sure I know what's going on! ****
if ~all(iCond==PiCond(perm,:)), error('SnPM:InvalidiCond', 'iCond~=PiCond(perm,:)'), end
%-Actual labelling must be at top of PiCond
if (perm~=1)
	PiCond(perm,:)=[];
	PiCond=[iCond;PiCond];
end
PiCond=[PiCond(1,:);PiCond(randperm(size(PiCond,1)-1)+1,:)];


%-Form non-null design matrix partitions (Globals handled later)
%=======================================================================
%-Form for HC computation at permutation perm
sHCform    = 'spm_DesMtx(PiCond(perm,:),''+0m'',''Cond'')';
%-Condition partition
[H,Hnames] = spm_DesMtx(iCond,'+0m','Cond');
%-Contrast of condition effects
% (spm_DesMtx puts condition effects in index order)
CONT       = [1,-1];
%-Include constant term in block partition 
B=ones(nScan,1); Bnames='Const';


%-Design description
%-----------------------------------------------------------------------
sDesign = sprintf('SingleSub: Two Sample T test; 2 conditions, replications: %d(cond)x%d(repl)',nCond,nRepl);
sPiCond = sprintf('%d permutations of conditions within exchangability blocks of size %d, bhPerms=%d',size(PiCond,1)*(bhPerms+1),Xblk,bhPerms);
