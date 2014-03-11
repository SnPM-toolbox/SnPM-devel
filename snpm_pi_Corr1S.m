% Mfile snpm_pi_Corr1S
% SnPM PlugIn design module - SingleSubj simple correlation
% SingleSub: Simple Regression (correlation); single covariate of interest
% FORMAT snpm_pi_Corr1S 
%
% See body of snpm_ui for definition of PlugIU interface.
%_______________________________________________________________________
%
% snpm_pi_Corr1S is a PlugIn for the SnPM design set-up program,
% creating design and permutation matrix appropriate for single-
% subject, correlation design.
%
%-Number of permutations
%=======================================================================
% There are nScan! (nScan factoral) possible permutations, where nScan
% is the number of scans of the subject. You can compute this using
% the gamma function in Matlab: nScan! is gamma(nScan+1); or by direct
% computation as prod(1:nScan)
%
%-Prompts
%=======================================================================
%
% 'Select scans in time order':  Enter the scans to be analyzed.
% It is important to input the scans in time order so that temporal
% effects can be accounted for.
%
% 'Enter covariate values': These are the values of the experimental
% covariate.
% 
% 'Size of exchangability block':  This is the number of adjacent
% scans that you believe to be exchangeable.  The most common cause 
% of nonexchangeability is a temporal confound (e.g. while 12 adjacent
% scans might not be free from a temporal effect, 4 adjacent scans
% could be regarded as having neglible temporal effect). See snpm.man
% for more information and references.
%
% '### Perms. Use approx. test?':  If there are a large number of
% permutations it may not be necessary (or possible!) to compute 
% all the permutations.  A common guideline is that 10,000 permutations
% are sufficient to characterize the permutation distribution well.
% More permutations will probably not change the results much.
% If you answer yes you will be prompted for the number...
% '# perms. to use?'
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
% C,Cnames      - Condition partition of design matrix, & effect names
% B,Bnames      - Block partition (constant term), & effect names
%
% - extra -
% CovInt        - Specified covariate of interest
% Xblk          - Size of exchangability block
%_______________________________________________________________________
% @(#)snpm_SSC.m	3.2 Andrew Holmes 04/06/08
%	$Id: snpm_pi_Corr1S.m,v 8.1 2009/01/29 15:02:57 nichols Exp $	

%-----------------------------functions-called------------------------
% spm_DesMtx
% spm_select
% spm_input
%-----------------------------functions-called------------------------


%-Initialisation
%-----------------------------------------------------------------------
iGloNorm = '123';		% Allowable Global norm. codes
sDesSave = 'CovInt Xblk';	% PlugIn variables to save in cfg file



%-Get filenames of scans
%-----------------------------------------------------------------------
P     = spm_select(Inf,'image','Select scans in time order');
nScan = size(P,1);

%-Get covariate
%-----------------------------------------------------------------------
CovInt = [];
while ~all(size(CovInt)==[nScan,1])
	CovInt = spm_input(sprintf('Enter covariate values [%d]',nScan),'+1');
	CovInt = CovInt(:);
end

%-Centre covariate if required
%-----------------------------------------------------------------------
%**** bCCovInt = (spm_input('Centre this covariate ?','+0','y/n')=='y');
bCCovInt = 1;	% Always centre covariate

%-Work out exchangability blocks
%-----------------------------------------------------------------------
%-Valid sizes of X-blk (Xblk) are integer divisors of nScan that are >1
tmp      = nScan./[nScan:-1:1];
Xblk     = tmp( floor(tmp)==ceil(tmp) & tmp>1 );
tmp      = int2str(Xblk(1));
for i=2:length(Xblk), tmp=str2mat(tmp,int2str(Xblk(i))); end
if length(Xblk)>1
  Xblk     = spm_input('Size of exchangability block','+0','b',tmp,Xblk);
end
nXblk    = (nScan/Xblk);
iXblk    = meshgrid(1:nXblk,1:Xblk); iXblk = iXblk(:)';


%-Work out how many perms, and ask about approximate tests
%-----------------------------------------------------------------------
%-NB: n! == gamma(n+1)
nPiCond = round(exp(nXblk*gammaln(Xblk+1)));
bAproxTst = spm_input(sprintf('%d Perms. Use approx. test?',nPiCond),...
							'+1','y/n')=='y';
if bAproxTst
    tmp = 0;
    while ((tmp>nPiCond) | (tmp==0) )
	tmp = spm_input(sprintf('# perms. to use? (Max %d)',nPiCond),'+0');
	tmp = floor(max([0,tmp]));
    end
    if (tmp==nPiCond), bAproxTst=0; else, nPiCond=tmp; end
end


%-Compute permutations of conditions
%=======================================================================

if bAproxTst
	%-Approximate test :
	% Build up random subset of all (within Xblk) permutations
	%===============================================================
	rand('seed',sum(100*clock))	%-Initialise random number generator
	PiCond      = zeros(nPiCond,nScan);
	PiCond(1,:) = 1+rem([0:Xblk*nXblk-1],Xblk);
	for i = 2:nPiCond
		%-Generate a new random permutation - see randperm
		[null,p] = sort(rand(Xblk,nXblk)); p = p(:)';
		%-Check it's not already in PiCond
		while any(all((meshgrid(p,1:i-1)==PiCond(1:i-1,:))'))
			[null,p] = sort(rand(Xblk,nXblk)); p = p(:)';
		end
		PiCond(i,:) = p;
	end
	clear p

else
	%-Full permutation test :
	% Build up exhaustive matrix of permutations
	%===============================================================
	%-Compute permutations for a single exchangability block
	%---------------------------------------------------------------
	%-Initialise XblkPiCond & remaining numbers
	XblkPiCond = [];
	lef = [1:Xblk]';
	%-Loop through numbers left to add to permutations, accumulating PiCond
	for i = Xblk:-1:1
		%-Expand XblkPiCond & lef
		tmp = round(exp(gammaln(Xblk+1)-gammaln(i+1)));
		Exp = meshgrid(1:tmp,1:i); Exp = Exp(:)';
		if ~isempty(XblkPiCond), XblkPiCond = XblkPiCond(:,Exp); end
		lef = lef(:,Exp);
		%-Work out sampling for lef
		tmp1 = round(exp(gammaln(Xblk+1)-gammaln(i+1)));
		tmp2 = round(exp(gammaln(Xblk+1)-gammaln(i)));
		sam = 1+rem(0:i*tmp1-1,i) + ([1:tmp2]-1)*i;
		%-Add samplings from lef to XblkPiCond
		XblkPiCond   = [XblkPiCond; lef(sam)];
		%-Delete sampled items from lef & condition size
		lef(sam) = [];
		tmp = round(exp(gammaln(Xblk+1)-gammaln((i-1)+1)));
		lef = reshape(lef,(i-1),tmp);
		%NB:gamma(Xblk+1)/gamma((i-1)+1) == size(XblkPiCond,2);
	end
	clear lef Exp sam i
	%-Reorientate so permutations are in rows
	XblkPiCond = XblkPiCond';

	%-Now build up the complete set of permutations: PiConds
	%---------------------------------------------------------------
	nXblkPiCond=size(XblkPiCond,1);
	PiCond=XblkPiCond;
	for i=2:nXblk
		tmp=size(PiCond,1);
		[iSup iInf] = meshgrid(1:nXblkPiCond,1:tmp);
		iSup = iSup(:)';
		iInf = iInf(:)';
		PiCond=[XblkPiCond(iSup,:),PiCond(iInf,:)];
	end
end

%-Check, condition and randomise PiCond
%-----------------------------------------------------------------------
%-Check PiConds sum within Xblks to sum of first nXblk natural numbers
if ~all(all(PiCond*spm_DesMtx(iXblk)== (Xblk+1)*Xblk/2 ))
	error('Invalid PiCond computed!'), end
%-Convert to full permutations from permutations within blocks
nPiCond = size(PiCond,1);
PiCond = PiCond + meshgrid((iXblk-1)*Xblk,1:nPiCond);
%-Randomise order of PiConds (except first) to allow interim analysis
rand('seed',sum(100*clock))	%-Initialise random number generator
PiCond=[PiCond(1,:);PiCond(randperm(nPiCond-1)+1,:)];
%-Check first permutation is null permutation
if ~all(PiCond(1,:)==[1:nScan])
	error('PiCond(1,:)~=[1:nScan]'); end


%-Form non-null design matrix partitions (Globals handled later)
%=======================================================================
CONT       = [1];		% Contrast of condition effects

%-Covariate partition & Form for HC computation at permutation perm
if bCCovInt
	C       = CovInt - mean(CovInt);
	Cnames  = 'CovInt(centred)';
	Cc      = CovInt;
	Ccnames = 'CovInt';
	sHCform = ['spm_DesMtx(CovInt(PiCond(perm,:))-mean(CovInt),',...
			'''C'',''CovInt(centred)'')'];
else
	C       = CovInt;
	Cnames  = 'CovInt';
	sHCform = 'spm_DesMtx(CovInt(PiCond(perm,:)),''C'',''CovInt'')';
end

%-Include constant term in block partition 
B=ones(nScan,1); Bnames='Const';


%-Design description
%-----------------------------------------------------------------------
sDesign = sprintf('SingleSub: Simple Regression (correlation); single covariate of interest: %d scans',nScan);
sPiCond = sprintf('Permutations of covariate within exchangability blocks of size %d, %d permutations',Xblk,size(PiCond,1));
