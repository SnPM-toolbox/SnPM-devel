% Mfile snpm_pi_PairTrand
% PlugIn design module for snpm_ui
% - MultiSubj simple activation design: 2 (cond)
% MultiSub: 2 conditions, replications - permutation test
% FORMAT snpm_pi_PairTrand 
%
% See body of snpm_ui for definition of PlugIn interface.
%_______________________________________________________________________
%
% snpm_pi_PairTrand is a PlugIn for the SnPM design set-up program,
% creating design and permutation matrix appropriate for multi-
% subject, two condition with replication design, where where
% the condition labels are permuted over subject.
%
% This PlugIn is only designed for studies where there are two
% sets of labels, each the A/B complement of the other 
% (e.g. ABABAB and BABABA), where half of the subjects have one
% labeling and the other half have the other labeling.  Of course,
% there must be an even number of subjects.
%
%-Number of permutations
%=======================================================================
% There are (nSubj)-choose-(nSubj/2) possible permutations, where 
% nSubj is the number of subjects.  Matlab doesn't have a choose
% function but you can use this expression
%
%       prod(1:nSubj)/prod(1:nSubj/2)^2
% 
%
%-Prompts
%=======================================================================
% 
% '# subjects': Number of subjects to analyze
%
% '# replications per condition':  Here you specify how many times
% each of the two conditions were repeated.
% 
% For each subject you will be prompted:
%
% 'Subject #: Select scans in time order': Enter this subject's scans.
% It is important to input the scans in time order so that temporal
% effects can be accounted for.
%
% 'Enter conditions index: (B/A)':  Using A's to indicate activation
% scans and B's to indicate baseline, enter a sequence of 2*nRepl
% letters, nRepl A's and nRepl B's, where nRepl is the number of 
% replications. Spaces are permitted.
% There can only be two possible condition indicies:  That of the
% first subject and the A<->B flip of the first subject.
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
% H,Hnames      - Condition partition of design matrix, & effect names
% B,Bnames      - Block partition (exchangability blocks), & effect names
%
% - extra -
% iCond         - Condition indicator vector
% iRepl         - Replication indicator vector
% PiSubj        - +/-1 flip conditions indicator for subjects,
%                 relative to first subject.
%_______________________________________________________________________
% Based on snpm_SSA2x.m v1.2 by Andrew Holmes
% @(#)snpm_MSA2x.m	3.3  Thomas Nichols 04/06/08
%	$Id: snpm_pi_PairTrand.m,v 8.1 2009/01/29 15:02:57 nichols Exp $	

%-----------------------------functions-called------------------------
% spm_DesMtx
% spm_select
% spm_input
%-----------------------------functions-called------------------------

%-Initialisation
%-----------------------------------------------------------------------
nCond    = 2;			% Number of conditions
iGloNorm = '123';		% Allowable Global norm. codes
sDesSave = 'iCond iRepl PiSubj';
				% PlugIn variables to save in cfg file
                
if snpm_get_defaults('shuffle_seed')
    % Shuffle seed of random number generator
    try
        rng('shuffle');
    catch
        % Old syntax        
        rand('seed',sum(100*clock));
    end
end                

%-Get number of subjects
nSubj    = spm_input('# subjects','+1');
if (nSubj==1), error('SnPM:SingleSubject','Use single subject plug for single subjects'); end    

%-Get number of replications per condition - 2 x nRepl design
nRepl    = spm_input('# replications per condition','+1');
nScan    = nRepl*nCond;

%-Get filenames and iCond, the condition labels
%=======================================================================
P     = [];
iCond = [];
iSubj = [];
iSubjC= [];
iRepl = [];
for subj=1:nSubj
    tmp = ['Subject ',int2str(subj),': Select scans in time order'];
    P = str2mat(P,spm_select(nCond*nRepl,'image',tmp));
    Cond=[];
    while(isempty(Cond))
	tmp=['Enter conditions index: (A/B) [',int2str(nCond*nRepl),']'];
	tmpCond = spm_input(tmp,'+0','s');
	%-Convert a/b notation to +/- vector    
	tmpCond(isspace(tmpCond)) = [];
	tmpCond = abs(upper(tmpCond));
	tmp = tmpCond;
	tmpCond = tmpCond-min(tmpCond);
	tmpCond = tmpCond/max([1,tmpCond])*2-1;
	%-Check validity of tmpCond
	if length(tmpCond)==nScan
	    if length(find(diff(sort(tmpCond)))) ~= nCond-1
		error('SnPM:InvalidnCond', 'Exactly ',[int2str(nCond), ...
			' conditions must be supplied']);
	    elseif sum(tmpCond)~=0
		error('SnPM:InvalidnRepl',['Exactly ',int2str(nRepl),' As and ', ...
			    int2str(nRepl),' Bs must be supplied']);
	    elseif isempty(iCond)
		sCond=setstr(tmp([1,diff(sort(tmp))]~=0));
		Cond = tmpCond;		
	    elseif any(iCond(1:nScan)~=tmpCond) && ...
			any(iCond(1:nScan)~=(-tmpCond))
		error('SnPM:InvalidiCond',['Conditions index must be same as ', ...
			    'first subject, or flipped']);
	    else		
		Cond = tmpCond;		
	    end
	else
	    error('SnPM:InvalidnIndices',['Enter indicies for ',int2str(nCond*nRepl),' scans'])
	end
    end
    iCond = [iCond, Cond];
    iSubj = [iSubj, subj*ones(1,nRepl*nCond)];
    iSubjC= [iSubjC, Cond(1)];
    iRepl = [iRepl, cumsum(Cond==-1).*(Cond==-1) + ...
		cumsum(Cond==1).*(Cond==1)];
end
% Force Subj 1 to have +1 label
if (iSubjC(1)==-1)
    iSubjC = -iSubjC;
    iCond  = -iCond;
end    
P(1,:) = [];
iSUBJ = iSubj;

%-Check for randomized design
%-----------------------------------------------------------------------
% "Flip"ness is defined relative to subject 1's ordering
nFlip = sum(iSubjC==-1);
if nFlip==nSubj/2
    % Balanced design... everything is cool
else
    % not balanced
    warning('SnPM:Unbalanced', ...
            'Design not balanced; results may be confounded with time');
end

%-Compute permutations of subjects
%=======================================================================
%-Work out how many perms, and ask about approximate tests
%-----------------------------------------------------------------------
nPiSubj = prod(nSubj-nFlip+1:nSubj)/prod(1:nFlip); % nSubj-choose-nFlip
if (spm_input(sprintf('%d Perms. Use approx. test?',nPiSubj),'+1','y/n')=='y')
    bAproxTst = 1;
    tmp = 0;
    while ((tmp>nPiSubj) || (tmp==0))
	tmp = spm_input(sprintf('# perms. to use? (Max %d)',nPiSubj),'+0');
	tmp = floor(max([0,tmp]));
	if rem(tmp,2)
	    error('SnPM:OddPermutations',['Number of perms must be even']);
	    tmp=0;	    
	end	    
    end
    if (tmp==nPiSubj)
	bAproxTst = 0;
    else
	nPiSubj=tmp;
    end
else
    bAproxTst = 0;
end

%-Compute permutations of subjects
%=======================================================================
%-Generate all labelings of subjects as +/-1
% NB: we take the first subject to be unflipped (+1) and work out the
% perms for the other n-1 subjs; this saves throwing away half of our
% work later on.
PiSubj=[];
for i=0:nSubj-2
    PiSubj=[ones(2^i,1),PiSubj;-ones(2^i,1),PiSubj];
end

% Keep only those that specify flipping nFlip
PiSubj(sum((PiSubj==-1)')~=nFlip,:) = [];

% Move correct labeling to top
d = find(all((PiSubj==meshgrid(iSubjC(2:nSubj),1:size(PiSubj,1)))'));
if (length(d)~=1)
    error('SnPM:CorrectLabelMissing', 'Internal error: Correct labeling is not in the perms'); end
PiSubj(d,:) = PiSubj(1,:);
PiSubj(1,:) = iSubjC(2:nSubj);
% Pick a random sample of PiSubj if bAprxTst
if bAproxTst
    tmp = randperm(size(PiSubj,1)-1);
    tmp = [1,tmp(1:nPiSubj/2-1)+1];
    PiSubj = PiSubj(tmp,:);
end

%-Slap on "unflipped" column for first subject
%----------------------------------------------------------------------
PiSubj=[ones(size(PiSubj,1),1),PiSubj];
bhPerms=1;

%-Build conditions perumutions - PiCond
%=======================================================================
% is there a better way?
tmp1=zeros(nSubj*nScan,nSubj);
tmp2=meshgrid(1:nScan,1:nSubj)';
tmp3=meshgrid(1:nSubj,1:nScan); tmp3 = tmp3(:)';
tmp1((tmp3-1)*(nSubj+1)*nScan + tmp2(:)') = iCond(tmp2(:));
PiCond = PiSubj*tmp1';
clear tmp1 tmp2 tmp3

%-Build correct perm design matrix, partition B
%=======================================================================
%-Use implicit SumToZero constraints via relative block effects & pinv.
%-See spm_DesMtx for more information on this.
[B,Bnames] = spm_DesMtx(iSUBJ,'+0m','Subj');

%-Form non-null design matrix partitions (Globals handled later)
%=======================================================================
%-Form for HC computation at permutation perm
sHCform    = 'spm_DesMtx(PiCond(perm,:),''-'',''Cond'')';
%-Condition partition
[H] = spm_DesMtx(iCond,'-','Cond');
Hnames=[];
for i=1:nCond
    Hnames=str2mat(Hnames,['Cond_',sCond(i)]); end
Hnames(1,:) = [];    
%-Contrast of condition effects
% (spm_DesMtx puts condition effects in index order)
CONT       = [1,-1];

%-Design description
%-----------------------------------------------------------------------
sDesign = sprintf('MultiSub: 2 conditions, replications - permutation test: %d(subj) %d(cond)x%d(repl)%s', ...
	nSubj,nCond,nRepl);
sPiCond = sprintf('Permutations of conditions by subject, %d permutations, bhPerms=%d', ...
	size(PiCond,1)*(bhPerms+1),bhPerms);
