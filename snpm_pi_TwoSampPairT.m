% Mfile snpm_pi_TwoSampPairT.m
% SnPM PlugIn design module - 2 group, 2 condition, 1 scan per condition
% 2 Groups: Test diff of response; 2 conditions, 1 scan per condition
% FORMAT snpm_pi_TwoSampPairT.m
%
% See body of snpm_ui for definition of PlugIn interface.
%_______________________________________________________________________
%
% snpm_pi_TwoSampPairT.m is a PlugIn for the SnPM design set-up program,
% creating design and permutation matrix appropriate for two group, two
% condition interaction analyses where there is just *one* scan per
% condition.  Only the interaction between activation and group is examined.
% Use snpm_MSA2x to examine intragroup (main) effects.
%
%
%-Number of permutations
%=======================================================================
%
% There are nTot-choose-nGrp possible permutations, where 
% nTot is the total number of subjects (and scans) and nGrp is the
% size of one of the groups.
%
%	prod(1:nTot)/prod(1:nGrp)^2
% 
%
%-Prompts
%=======================================================================
%
% '# of subjects in group <num> ?': Enter the number of subjects in the
% specified group.
%
% 'Select scans, group <num>, subj <num>':  Select the two scans for the
% specified subject.  Enter the scans in time order.
%
% 'Enter scan index: (AB|BA)':  Assuming you consistently label one
% condition 'A', one condition 'B', enter AB or BA as the labels for the two
% scans just entered.  This PlugIn computes the interaction as
% (A1-B1)-(A2-B2) (where A1 is group 1 condition A, etc) though you can get
% both positive and negative effects in % the results module.
%
% '# of confounding covariates' & '[<len>] - Covariate <num>': Use these
% prompts to specify a covariate of no interest.
%
% '<nPerms> Perms.  Use approx. test':  This prompt will inform you of the
% number of possible permutations, that is, the number of ways the group
% labels can be arranged under the assumption that there is no group
% effect.  Fewer than 200 permutations is undesirable; more than 10,000
% is unnecessary.  If the number of permutations is much greater than 10,000
% you should use an approximate test.  Answering 'y' will produce another
% prompt... 
% '# perms. to use? (Max <MaxnPerms>)': 10,000 permutations is regarded as
% a sufficient number to characterize the permutation distribution well.
%
% 
%-Exchangeability
%=======================================================================
% 
% There are two exhangeability issues:  Subject exchangeability and scan
% exchagebility. 
% 
% If the "AB" order is not randomized then one must assume that the data are
% exchangeable under the null hypothesis of no main effect.  Essentially, if
% the order is not randomized you must assume there is no temporal efffect.  
%
% If subjects are not randomized to group, one must assume that that the
% subjects are exchangeable under the null hypothesis, that is that none of
% them are a priori "special". 
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
% iGrp          - Group indicator vector
% nSubGrp       - 2-vector of group counts (subjects per group)
% 
%_______________________________________________________________________
% Based on snpm_MG2x.m v1.5
% @(#)snpm_MG2i.m	3.3 Thomas Nichols & Andrew Holmes 04/06/08
%	$Id: snpm_pi_TwoSampPairT.m,v 8.1 2009/01/29 15:02:57 nichols Exp $	

%-----------------------------functions-called------------------------
% spm_DesMtx
% spm_select
% spm_input
%-----------------------------functions-called------------------------


%-Initialisation
%-----------------------------------------------------------------------
nCond    = 2;			% Number of conditions
nStud    = 2;			% Number of groups
iGloNorm = '123';		% Allowable Global norm. codes
sDesSave = 'iStud iCond nSubj'; % PlugIn variables to save in cfg file
global TEST;
if isempty(TEST) || ~TEST % When testing we need a fixed seed
    rand('seed',sum(100*clock));	% Initialise random number generator
end
iStudC   = [];			% +1/-1 version of iStud

%-Get filenames and iStud, the subject group labels
%=======================================================================
P = '';
for stud=1:nStud
    tmp = sprintf('# of subjects in group %d ?',stud);
    nSubj = spm_input(tmp,'+0');
    for subj=1:nSubj    
	tmp = sprintf('Select scans, group %d, subj %d ',stud,subj);
	tP = spm_select(2,'image',tmp);
	% get condition label	
	Cond = [];    
	while(isempty(Cond))
	    tmp='Enter scan index: (AB|BA)';
	    tmpCond = upper(spm_input(tmp,'+0','s'));
	    if (strcmp(tmpCond,'AB'))
		Cond = [+1 -1];
	    elseif (strcmp(tmpCond,'BA'))
		Cond = [-1 +1];
	    else	    
		fprintf(2,'%cEnter either AB or BA',7);
	    end	    
	end
	P = str2mat(P,tP);
	% update indicators	
	iCond = [iCond, Cond];
	iSubj = [iSubj, subj*ones(1,nCond)];
	iStudC = [iStudC, 3-2*stud];
	iStud  = [iStud, stud*ones(1,nCond)];
    end
end
P(1,:)=[];

nFlip   = sum(iStudC==-1);
nSubj   = iSubj([diff(iStud),1]~=0);                    %-#subject per study
nSUBJ   = sum(nSubj);                                   %-#subjects in total
tmp     = cumsum([0,nSubj]);
iSUBJ   = iSubj+tmp(cumsum([1,diff(iStud)]));           %-Index to subjects
nScan   = nSubj*nCond;

%-Get confounding covariates
%-----------------------------------------------------------------------
G = []; Gnames = ''; Gc = []; Gcnames = ''; q = nScan;
g = spm_input('# of confounding covariates','+1','0|1|2|3|4|5|>',0:6,1);
if (g == 6), g = spm_input('# of confounding covariates','+1'); end
while size(Gc,2) < g
  nGcs = size(Gc,2);
  d = spm_input(sprintf('[%d] - Covariate %d',[q,nGcs + 1]),'0');
  if (size(d,1) == 1), d = d'; end
  if size(d,1) == q
    %-Save raw covariates for printing later on
    Gc = [Gc,d];
    %-Always Centre the covariate
    bCntr = 1;	    
    if bCntr, d  = d - ones(q,1)*mean(d); str=''; else, str='r'; end
    G = [G, d];
    dnames = [str,'ConfCov#',int2str(nGcs+1)];
    for i = nGcs+1:nGcs+size(d,1)
      dnames = str2mat(dnames,['ConfCov#',int2str(i)]); end
    Gcnames = str2mat(Gcnames,dnames);
  end
end
%-Strip off blank line from str2mat concatenations
if size(Gc,2), Gcnames(1,:)=[]; end
%-Since no FxC interactions these are the same
Gnames = Gcnames;


%-Compute permutations of subjects
%=======================================================================
%-Compute permutations for a single exchangability block
%-----------------------------------------------------------------------
%-NB: m-Choose-n = exp(gammaln(m+1)-gammaln(m-n+1)-gammaln(n+1))
nPiStud = round(exp(gammaln(nSUBJ+1)-gammaln(nSUBJ-nFlip+1)-gammaln(nFlip+1)));
bAproxTst = spm_input(sprintf('%d Perms. Use approx. test?',nPiStud),...
							'+1','y/n')=='y';
if (bAproxTst)
  tmp = 0;
  while ((tmp>nPiStud) | (tmp==0) )
    tmp = spm_input(sprintf('# perms. to use? (Max %d)',nPiStud),'+0');
    tmp = floor(max([0,tmp]));
  end
  if (tmp==nPiStud), bAproxTst=0; else, nPiStud=tmp; end
end

%-Two methods for computing permutations, random and exact; exact
% is efficient, but a memory hog; Random is slow but requires little
% memory.
%-We use the exact one when the nSUBJ is small enough; for nSUBJ=12,
% PiStud will initially take 384KB RAM, for nSUBJ=14, 1.75MB, so we 
% use 12 as a cut off. (2^nSUBJ*nSUBJ * 8bytes/element).  
%-If user wants all perms, then random method would seem to take an
% absurdly long time, so exact is used.

if nSUBJ<=12 | ~bAproxTst                    % exact method

    %-Generate all labellings of nSUBJ subjects as +/- 1
    PiStud=[];
    for i=0:nSUBJ-1
	PiStud=[ones(2^i,1),PiStud;-ones(2^i,1),PiStud];
    end
    %-Trim to labellings with correct group numbers
    PiStud=PiStud(sum(PiStud'==-1)==nFlip,:);

    %-Only do half the work, if possible
    bhPerms=0;
    if ~bAproxTst & (nFlip==nSUBJ/2) % balanced group numbers
	% Here, PiStud should *always* satisfy:
	% all(all(PiStud(PiStud(:,1)==1,:)==flipud(-PiStud(PiStud(:,1)==-1,:))))
	PiStud=PiStud(PiStud(:,1)==1,:);
	bhPerms=1;
    elseif bAproxTst                 % pick random supsample of perms
	tmp=randperm(size(PiStud,1));
	PiStud=PiStud(tmp(1:nPiStud),:);
        % Note we may have missed iStudC!  We catch this below.	
    end	

else                                          % random method
    
    % Allocate final result
    PiStud = zeros(nPiStud,nSUBJ);

    % Fill first row  
    PiStud(1,:) = iStudC;
    % Fill subsequent rows, checking that we're not repeating  
    for i=2:nPiStud
      tmp=PiStud(i-1,randperm(nSUBJ));
      while any(all(PiStud(1:(i-1),:)'==meshgrid(tmp,1:(i-1))'))
	tmp=PiStud(i-1,randperm(nSUBJ));
      end
      PiStud(i,:)=tmp;
    end      

    bhPerms=0;    
end

%-Check each perm in PiStuds sums to nStud1-nStud2
if ~all(all(PiStud*ones(nSUBJ,1)==nSubj*[1 -1]'))
	error('Invalid PiStud computed!'), end

%-Find (maybe) iStudC in PiStud, move iStudC to 1st; negate if neccesary
%-----------------------------------------------------------------------
perm = find(all((meshgrid(iStudC,1:size(PiStud,1))==PiStud)'));
if (bhPerms)
    perm=[perm,-find(all((meshgrid(iStudC,1:size(PiStud,1))==-PiStud)'))];
end
if length(perm)==1
    if (perm<0), PiStud=-PiStud; perm=-perm; end
    %-Actual labelling must be at top of PiStud
    if (perm~=1)
	PiStud(perm,:)=[];
	PiStud=[iStudC;PiStud];
    end
    if ~bAproxTst    
        %-Randomise order of PiStuds, unless already randomized
        % Allows interim analysis	
	PiStud=[PiStud(1,:);PiStud(randperm(size(PiStud,1)-1)+1,:)];
    end	
elseif length(perm)==0 & (nScan<=12) & bAproxTst
    % Special case where we missed iStud; order of perms is random 
    % so can we can just replace first perm.
    PiStud(1,:) = iStudC;
    perm = 1;
else
    error(['Bad PiStud (' num2str(perm) ')'])
end    

%-Turn PiStud into PiCond
% Expand subject labels into scan labels, convert [1,-1]->[1 2]
% Note snpm_cp uses PiCond to determine number of permutations
tmp=[ones(2,nSUBJ); zeros(2*(nSUBJ),nSUBJ)];
tmp=reshape(tmp(1:2*nSUBJ*nSUBJ),2*nSUBJ,nSUBJ)';
PiCond = 1.5-PiStud*tmp/2;

%-Build correct perm design matrix, partition B
%=======================================================================
%-Use implicit SumToZero constraints via relative block effects & pinv.
%-See spm_DesMtx for more information on this.
[B Bnames] = spm_DesMtx(iSUBJ,'+0m','Subj');


%-Form non-null design matrix partitions (Globals handled later)
%=======================================================================
%-Form for HC computation at permutation perm
sHCform = 'spm_DesMtx([PiCond(perm,:)'' iCond''],''-'',[''Group'';''Cond ''])';
%-Condition partition
[H,Hnames] = spm_DesMtx([iStud' iCond'],'-',['Group';'Cond ']);
%-Contrast of condition effects
% (spm_DesMtx puts condition effects in index order)
CONT       = [1,-1,-1,1];


%-Design description
%-----------------------------------------------------------------------
sDesign = sprintf('2 Groups: Test diff of response; 2 conditions, 1 scan per condition: %d(GrpA),%d(GrpB)',nSubj);
sPiCond = sprintf('%d permutations of conditions, bhPerms=%d',size(PiCond,1)*(bhPerms+1),bhPerms);
