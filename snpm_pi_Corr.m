% Mfile snpm_pi_Corr
% SnPM PlugIn design module - MultiSubj simple correlation, 1 scan per subj.
% MultiSub: Simple Regression (correlation); single covariate of interest, 1 scan per subject
% FORMAT snpm_pi_Corr 
%
% See body of snpm_ui for definition of PlugIn interface.
%_______________________________________________________________________
%
% snpm_pi_Corr is a PlugIn for the SnPM design set-up program,
% creating design and permutation matrix appropriate for multi-
% subject, correlation design, where there is only one scan per subject.
%
%-Number of permutations
%=======================================================================
% There are nSubj! (nSubj factoral) possible permutations, where nSubj
% is the number of subjects. You can compute this using the gamma
% function in Matlab: nSubj! is gamma(nSubj+1); or by direct 
% computation as prod(1:nSubj)
%
%-Prompts
%=======================================================================
%
% 'Select scans':  Enter the scans to be analyzed.
%
% 'Enter covariate values': These are the values of the experimental
% covariate.
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
% nSubj          - Size of exchangability block
%_______________________________________________________________________
% Copyright (C) 2013 The University of Warwick
% Id: snpm_pi_Corr.m  SnPM13 2013/10/12
% Thomas Nichols & Andrew Holmes, Emma Thomas
% Based on snpm_SSC v3.2, 04/06/08

%-----------------------------functions-called------------------------
% spm_DesMtx
% spm_select
% spm_input
%-----------------------------functions-called------------------------


%-Initialisation
%-----------------------------------------------------------------------
iGloNorm = '123';		% Allowable Global norm. codes
sDesSave = 'CovInt';		% PlugIn variables to save in cfg file

global TEST;

%-Get filenames of scans
%-----------------------------------------------------------------------
if BATCH
  P = strvcat(job.P);  % Sill to use string arrays, but need it for compatibility
else
  P     = spm_select(Inf,'image','Select scans, 1 per subj.');
end
nSubj = size(P,1);

%-Get covariate
%-----------------------------------------------------------------------
if BATCH
  CovInt = job.CovInt(:);
  if ~all(size(CovInt)==[nSubj,1])
    error(sprintf('Covariate [%d,1] doesn''t match number of subjects [%d]',...
		  size(CovInt,1),nSubj))
  end
else
  CovInt = [];
  while ~all(size(CovInt)==[nSubj,1])
    CovInt = spm_input(sprintf('Enter covariate values [%d]',nSubj),'+1');
    CovInt = CovInt(:);
  end
end

%-Centre covariate if required
%-----------------------------------------------------------------------
%**** bCCovInt = (spm_input('Centre this covariate ?','+0','y/n')=='y');
bCCovInt = 1;	% Always centre covariate

%-Work out how many perms, and ask about approximate tests
%-----------------------------------------------------------------------
%-NB: n! == gamma(n+1)
nPiCond_mx = gamma(nSubj+1);
fprintf('\nNOTE: Number of possible permutations for this design is %d\n',nPiCond_mx)
if BATCH
  if job.nPerm >= nPiCond_mx
    bAproxTst=0;
    nPiCond = nPiCond_mx;
    fprintf('NOTE: Requested %d perms, but only using %d (maximum possible).\n',job.nPerm,nPiCond_mx);
  else
    bAproxTst=1;
    nPiCond = job.nPerm;
  end
else
  bAproxTst = spm_input(sprintf('%d Perms. Use approx. test?',nPiCond_mx),...
			'+1','y/n')=='y';
  if bAproxTst
    tmp = 0;
    while ((tmp>nPiCond_mx) | (tmp==0) )
      tmp = spm_input(sprintf('# perms. to use? (Max %d)',nPiCond_mx),'+0');
      tmp = floor(max([0,tmp]));
    end
    nPiCond=tmp; 
    if (tmp==nPiCond_mx), bAproxTst=0; end
  else
    nPiCond=nPiCond_mx;
  end
end
snpm_check_nperm(nPiCond,nPiCond_mx);


%-Compute permutations of conditions
%=======================================================================

if bAproxTst
	%-Approximate test :
	% Build up random subset of all (within nSubj) permutations
	%===============================================================
    if isempty(TEST) || ~TEST % When testing the code we need a fixed seed
        rand('seed',sum(100*clock))	%-Initialise random number generator
    end
	PiCond      = zeros(nPiCond,nSubj);
	PiCond(1,:) = 1+rem([0:nSubj-1],nSubj);
	for i = 2:nPiCond
		%-Generate a new random permutation - see randperm
		[null,p] = sort(rand(nSubj,1)); p = p(:)';
		%-Check it's not already in PiCond
		while any(all((meshgrid(p,1:i-1)==PiCond(1:i-1,:))'))
			[null,p] = sort(rand(nSubj,1)); p = p(:)';
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
	lef = [1:nSubj]';
	%-Loop through numbers left to add to permutations, accumulating PiCond
	for i = nSubj:-1:1
		%-Expand XblkPiCond & lef
		tmp = round(exp(gammaln(nSubj+1)-gammaln(i+1)));
		Exp = meshgrid(1:tmp,1:i); Exp = Exp(:)';
		if ~isempty(XblkPiCond), XblkPiCond = XblkPiCond(:,Exp); end
		lef = lef(:,Exp);
		%-Work out sampling for lef
		tmp1 = round(exp(gammaln(nSubj+1)-gammaln(i+1)));
		tmp2 = round(exp(gammaln(nSubj+1)-gammaln(i)));
		sam = 1+rem(0:i*tmp1-1,i) + ([1:tmp2]-1)*i;
		%-Add samplings from lef to XblkPiCond
		XblkPiCond   = [XblkPiCond; lef(sam)];
		%-Delete sampled items from lef & condition size
		lef(sam) = [];
		tmp = round(exp(gammaln(nSubj+1)-gammaln((i-1)+1)));
		lef = reshape(lef,(i-1),tmp);
		%NB:gamma(nSubj+1)/gamma((i-1)+1) == size(XblkPiCond,2);
	end
	clear lef Exp sam i
	%-Reorientate so permutations are in rows
	XblkPiCond = XblkPiCond';
	PiCond=XblkPiCond;
end

%-Check, condition and randomise PiCond
%-----------------------------------------------------------------------
%-Check PiConds sum within Xblks to sum to 1
if ~all(all(sum(PiCond,2)== (nSubj+1)*nSubj/2 ))
	error('Invalid PiCond computed!'), end
%-Convert to full permutations from permutations within blocks
nPiCond = size(PiCond,1);
%-Randomise order of PiConds (except first) to allow interim analysis
if isempty(TEST) || ~TEST % When testing the code we need a fixed seed
    rand('seed',sum(100*clock))	%-Initialise random number generator
end
PiCond=[PiCond(1,:);PiCond(randperm(nPiCond-1)+1,:)];
%-Check first permutation is null permutation
if ~all(PiCond(1,:)==[1:nSubj])
	error('PiCond(1,:)~=[1:nSubj]'); end


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
B=ones(nSubj,1); Bnames='Const';


%-Design description
%-----------------------------------------------------------------------
sDesign = sprintf('MultiSub: Simple Regression (correlation); single covariate of interest, 1 scan per subject: %d subj',nSubj);
sPiCond = sprintf('Permutations of covariate: %d',size(PiCond,1));
