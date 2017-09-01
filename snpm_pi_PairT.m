% Mfile snpm_pi_PairT
% PlugIn design module for snpm_ui
% - MultiSubj nonrandomized act design: 2 (cond) no replications
% MultiSub: Paired T test; 2 conditions, 1 scan per condition
% FORMAT snpm_pi_PairT
%
% See body of snpm_ui for definition of PlugIn interface.
%_______________________________________________________________________
%
%
% snpm_pi_PairT is a PlugIn for the SnPM design set-up program,
% creating design and permutation matrix appropriate for multi-
% subject, two condition no replication design, where the condition
% labels have been applied to the subjects in a nonrandomized fasion.  
%
% With only two scans per subject, there are only two possible sets of
% labels: AB and BA.  No restriction is placed on how many subjects
% received A first.  That is, unlike snpm_MSA2x, this PlugIn handles
% designs where all subjects receive condition A first. 
%
% Since there is no randomization we must justify exchangibility with
% assumptions under the null hypothesis.  In particular, we must assume
% that at each voxel, the distribution of the data is the same for A and
% B scans and for all subjects.  This is equivalent to assuming that
% distribution of A-B is symmetric and the same for all subjects. Note
% that that an unmodeled temporal effect would violated this assumption
% (no such assumption is needed when a randomization design is used; if
% randomization was used snpm_MSA2x should be used). 
%
% (If there are replications, it is recommended that a "first level")
% (model is fit, reducing each condition to a single summary image. )
%
%-Number of permutations
%=======================================================================
% There are 2^nSubj possible labelings, where nSubj is the number of 
% subjects.  This is because each subject can have two possible states,
% flipped or unflipped.
%
% For example, 2^7 = 128 and 2^8 = 256.  Hence at least eight subjects are
% needed to characterize the permutation distribution well, and 10 or more
% are best.
%
%-Prompts
%=======================================================================
% 
% '# subjects': Number of subjects to analyze
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
% Copyright (C) 2013 The University of Warwick
% Id: snpm_pi_PairT.m  SnPM13 2013/10/12
% Thomas Nichols, Camille Maumet
% Based on snpm_MSA2x.m v1.5

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

if snpm_get_defaults('shuffle_seed')
    % Shuffle seed of random number generator
    try
        rng('shuffle');
    catch
        % Old syntax        
        rand('seed',sum(100*clock));
    end
end

				% PlugIn variables to save in cfg file

%-Get number of subjects
% nSubj    = spm_input('# subjects','+1');
% if (nSubj==1), error('SnPM:SingleSubject', 'Use single subject plug for single subjects'); end    
nSubj = numel(job.fsubject);

%-Only consider one replication -- basically a RFX machine.
nRepl  = 1;

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
    P = str2mat(P, str2mat(job.fsubject(subj).scans)); %str2mat(P,spm_select(nCond*nRepl,'image',tmp));
    Cond=[];
    while(isempty(Cond))
%         tmp=['Enter conditions index: (A/B) [',int2str(nCond*nRepl),']'];
%         tmpCond = spm_input(tmp,'+0','s');
%         %-Convert a/b notation to +/- vector    
%         tmpCond(isspace(tmpCond)) = [];
%         tmpCond = abs(upper(tmpCond));
%         tmp = tmpCond;
        tmpCond = job.fsubject(subj).scindex;
        tmpCond = tmpCond-min(tmpCond);
        tmpCond = tmpCond/max([1,tmpCond])*2-1;
        %-Check validity of tmpCond
        if length(tmpCond)==nScan
            if length(find(diff(sort(tmpCond)))) ~= nCond-1
            error('SnPM:InvalidnCond', 'Exactly ',[int2str(nCond), ...
                ' conditions must be supplied']);
            elseif sum(tmpCond)~=0
            error('SnPM:InvalidnRepl', ['Exactly ',int2str(nRepl),' As and ', ...
                    int2str(nRepl),' Bs must be supplied']);
            elseif isempty(iCond)
            sCond=setstr(tmp([1,diff(sort(tmp))]~=0));
            Cond = tmpCond;		
            elseif any(iCond(1:nScan)~=tmpCond) && ...
                any(iCond(1:nScan)~=(-tmpCond))
            error('SnPM:InvalidiCond', ['Conditions index must be same as', ...
                    'first subject, or flipped']);
            else		
            Cond = tmpCond;		
            end
        else
            error('SnPM:InvalidiIndices', ['Enter indicies for ',int2str(nCond*nRepl),' scans'])
        end
    end
    iCond = [iCond, Cond];
    iSubj = [iSubj, subj*ones(1,nScan)];
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

%-Compute permutations of subjects
%=======================================================================
%-Work out how many perms, and ask about approximate tests
%-----------------------------------------------------------------------
if nSubj <= 52
  nPiSubj_mx = 2^nSubj;
else
  nPiSubj_mx = Inf;
end
nPiSubj = job.nPerm;
if job.nPerm >= nPiSubj_mx
    bAproxTst=0;
    if job.nPerm > nPiSubj_mx
        fprintf('NOTE: %d permutations requested, only %d possible.\n',job.nPerm, nPiSubj_mx)
        nPiSubj = nPiSubj_mx;
    end
else
    bAproxTst=1;
end
if rem(nPiSubj,2)
    error('SnPM:OddPermutations', ['Number of perms must be even']);
    nPiSubj = 0;	    
end	  

% if (spm_input(sprintf('%d Perms. Use approx. test?',nPiSubj),'+1','y/n')=='y')
%     bAproxTst = 1;
%     tmp = 0;
%     while ((tmp>nPiSubj) | (tmp==0))
% 	tmp = spm_input(sprintf('# perms. to use? (Max %d)',nPiSubj),'+0');
% 	tmp = floor(max([0,tmp]));
% 	if rem(tmp,2)
% 	    error('SnPM:OddPermutations', ['Number of perms must be even']);
% 	    tmp=0;	    
% 	end	    
%     end
%     if (tmp==nPiSubj)
% 	bAproxTst = 0;
%     else
% 	nPiSubj=tmp;
%     end
% else
%     bAproxTst = 0;
% end

%-Compute permutations of subjects
%=======================================================================
%-All possible labelings correspond to the binary representation of
% numbers {1...2^nSubj}.
if nSubj<=35
  if (bAproxTst)
    tmp = randperm(2^nSubj)-1;
    tmp = tmp(1:nPiSubj)';
  else
    tmp = (0:(2^nSubj-1))';
  end
  %-Generate labelings of subjects as +/-1
  PiSubj=[];
  for i=(nSubj-1):-1:0
    PiSubj = [PiSubj,2*(tmp>=2^i)-1];
    tmp = tmp - (tmp>=2^i)*2^i;
  end
  % Look for correct labeling
  d = find(all((PiSubj==meshgrid(iSubjC,1:size(PiSubj,1)))'));
  if (length(d)~=1 && ~bAproxTst)
    error('SnPM:CorrectLabelMissing', 'Internal error: Correct labeling is not in the perms');
  elseif (length(d)~=1)
    % Correct labeling randomly removed, insert at top
    PiSubj(1,:) = iSubjC;
  else
    % Swap correct labeling to top
    PiSubj(d,:) = PiSubj(1,:);
    PiSubj(1,:) = iSubjC(1:nSubj);
  end
else
  % Here we are always approximate 
  PiSubj = [...
      iSubjC
      2*randi(2,nPiSubj-1,nSubj)-3];
end

%-If not approximate then we can just calc half
%----------------------------------------------------------------------
if ~bAproxTst
  PiSubj = PiSubj(PiSubj(:,1)==1,:);
  bhPerms=1;
else
  bhPerms=0;
end

%-Build conditions perumutions - PiCond
%=======================================================================
% is there a better way?
tmp1=zeros(nSubj*nScan,nSubj);
tmp2=meshgrid(1:nScan,1:nSubj)';
tmp3=meshgrid(1:nSubj,1:nScan); tmp3 = tmp3(:)';
tmp1((tmp3-1)*(nSubj+1)*nScan + tmp2(:)') = iCond(tmp2(:));
PiCond = PiSubj*tmp1';
clear tmp1 tmp2 tmp3

%-Change Cond def; Old: [1=A, -1=B]; New: [1=A, 2=B]
%=======================================================================
iCond(iCond==-1) = 2;
PiCond(PiCond==-1) = 2;

%-Build correct perm design matrix, partition B
%=======================================================================
%-Use implicit SumToZero constraints via relative block effects & pinv.
%-See spm_DesMtx for more information on this.
[B Bnames] = spm_DesMtx(iSUBJ,'+0m','Subj');

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
sDesign = sprintf('MultiSub: Paired T test; 2 conditions, 1 scan per condition: %d(subj) %d(cond)x%d(repl)%s', ...
	nSubj,nCond,nRepl);
sPiCond = sprintf('Permutations of conditions by subject, %d permutations, bhPerms=%d', ...
	size(PiCond,1)*(bhPerms+1),bhPerms);
