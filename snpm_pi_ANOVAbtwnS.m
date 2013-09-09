% Mfile snpm_pi_ANOVAbtwnS
% SnPM PlugIn design module - Between Subject ANOVA, k groups, 1 scan per subject
% FORMAT snpm_pi_ANOVAbtwnS 
%
% See body of snpm_ui for definition of PlugIn interface.
%_______________________________________________________________________
%
% snpm_pi_ANOVAbtwnS is a PlugIn for the SnPM design set-up program,
% creating design and permutation matrix appropriate for k groups
% analyses where there is just *one* scan per subject.
%
% This PlugIn can be regarded as a generalization of a two-sample t-test to
% k>2 groups.  Instead of producing a t-statistic, it produces a F statistic.  
%
% The PlugIn can test for two different types of effects.  It can test
% for the presence of any *non-zero* effect among the k groups; that is,
% it tests the null hypothesis that all of the groups are mean zero.  Or
% it can test for the presence of any differences *between* the k groups;
% that is, it tests the null hypothesis that all of the groups have some 
% (possibly non-zero) common mean.
%
%-Number of permutations
%=======================================================================
%
% There are (nScan)!/(GrpCnt[1]!*GrpCnt[2]!*...*GrpCnt[nCond]!)
% possible permutations, where nScan is the total number of scans and 
% GrpCnt(i) is the size of the ith group.
%
%     round(exp(gammaln(nScan+1)-sum(gammaln(GrpCnt+1))))
% 
%
%-Prompts
%=======================================================================
%
% 'Select all scans':  Enter the scans to be analyzed; the order 
% is not important as the specification of which scans belong to which
% groups will be specified subsequently.
%
% '# of confounding covariates' & '[<len>] - Covariate <num>': Use these
% prompts to specify a covariate of no interest.  As mentioned above,
% fitting a confounding covariate of age may be desirable.
%
% 'Enter Subject index: (A/B/...)':  Use A's, B's and et.c. to indicate which 
% scans belong to which group. 
% You must enter one letter for each scan entered above.
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
% GrpCnt        - A vector of group counts
%
%_______________________________________________________________________
% Based on snpm_MG2x.m v1.7
% $Id: snpm_pi_ANOVAbtwnS.m,v 8.1 2009/01/29 15:02:57 nichols Exp $	

%-----------------------------functions-called------------------------
% spm_DesMtx
% spm_select
% spm_input
%-----------------------------functions-called------------------------

% Programmer's note
%
% If the null hypothesis is that all means are zero, then an alternative
% permutation scheme is to flip the signs of the individual's data (or,
% possibly, flipping the sign of an entire group together?).  Need
% further evaluation to see which permutation scheme is best.
% 

%-Initialisation
%-----------------------------------------------------------------------
%%% nCond    = 2;			% Number of conditions (groups)
iGloNorm = '123';		% Allowable Global norm. codes
sDesSave = 'iCond GrpCnt';	% PlugIn variables to save in cfg file
rand('seed',sum(100*clock));	% Initialise random number generator

%-Get filenames and iCond, the condition labels
%=======================================================================
P = spm_select(Inf,'image','Select all scans');
nScan = size(P,1);

%-Get the condition (group) labels
%=======================================================================
while(1)
    nCond = spm_input('Number of groups k=','+0','w',3,1);
    if (nCond <= 2)
         fprintf(2,'%cNumber of groups should be greater than 2.',7)
    else
         break
    end
end

if nCond>255, error('Can''t support more than 255 groups'); end

tmp0='A/B/...';

while(1)
    tmp=['Enter subject index: (',tmp0, ')[',int2str(nScan),']'];
    iCond = spm_input(tmp,'+1','s');
    %-Convert A/B/C notation to 1,2,...,k vector - assume A-B is of interest
    iCond = abs(upper(iCond(~isspace(iCond))));
    iCond = iCond-min(iCond)+1;
    
    unique_iCond = unique(iCond);
    
    %-Check validity of iCond
    if length(iCond)~= nScan
        fprintf(2,'%cEnter indicies for exactly %d scans',7,nScan)
    elseif length(unique_iCond) ~= nCond
        fprintf(2,'%cEnter indicies for exactly %d groups',7,nCond) 	
    else
        % Deal with the 'Skip' situation, e.g. if users input 'A C A C D D',
        % Then iCond = [1 2 1 2 3 3];
        for i = 1:length(unique_iCond)
           iCond(iCond==unique_iCond(i)) = i;
        end
        
	GrpCnt = zeros(1,nCond);
        for (i = 1:nCond)
           GrpCnt(i) = sum(iCond==i);
        end
        break	
    end
end

%-Get the F contrasts
%-----------------------------------------------------------------------
b_all_zero = spm_input('Null Hypothesis: Groups are','+1','b','all zero|all equal',[1,0],1);

%-Get confounding covariates
%-----------------------------------------------------------------------
G = []; Gnames = ''; Gc = []; Gcnames = ''; q = nScan;
% g = spm_input('# of confounding covariates','+1','0|1|2|3|4|5|>',0:6,1);
% if (g == 6), g = spm_input('# of confounding covariates','+1'); end
% while size(Gc,2) < g
%   nGcs = size(Gc,2);
%   d = spm_input(sprintf('[%d] - Covariate %d',[q,nGcs + 1]),'0');
%   if (size(d,1) == 1), d = d'; end
%   if size(d,1) == q
%     %-Save raw covariates for printing later on
%     Gc = [Gc,d];
%     %-Always Centre the covariate
%     bCntr = 1;	    
%     if bCntr, d  = d - ones(q,1)*mean(d); str=''; else, str='r'; end
%     G = [G, d];
%     dnames = [str,'ConfCov#',int2str(nGcs+1)];
%     for i = nGcs+1:nGcs+size(d,1)
%       dnames = str2mat(dnames,['ConfCov#',int2str(i)]); end
%     Gcnames = str2mat(Gcnames,dnames);
%   end
% end
% %-Strip off blank line from str2mat concatenations
% if size(Gc,2), Gcnames(1,:)=[]; end
% %-Since no FxC interactions these are the same
% Gnames = Gcnames;


%-Compute permutations of conditions
%=======================================================================
%-Compute permutations for a single exchangability block
%-----------------------------------------------------------------------
%-NB: m-Choose-n = exp(gammaln(m+1)-gammaln(m-n+1)-gammaln(n+1))
%-NB: a! = exp(gammaln(a+1))
%-NB: nPiCond =
%(nScan)!/(GrpCnt[1]!*GrpCnt[2]!*...*GrpCnt[nCond]!)

nPiCond_mx = round(exp(gammaln(nScan+1)-sum(gammaln(GrpCnt+1))));

bAproxTst = spm_input(sprintf('%d Perms. Use approx. test?',nPiCond_mx),...
							'+1','y/n')=='y';
if (bAproxTst)
  tmp = 0;
  Defperm = min(10000,nPiCond);
  while ((tmp>nPiCond) | (tmp==0) )
    tmp = spm_input(sprintf('# perms. to use? (Max %d)',nPiCond),'+0','w',Defperm,1);
    tmp = floor(max([0,tmp]));
  end
  nPiCond=tmp; 
  if (tmp==nPiCond), bAproxTst=0; end
else
  nPiCond=nPiCond_mx;
end
snpm_check_nperm(nPiCond,nPiCond_mx);


%-Two methods for computing permutations, random and exact; exact
% is efficient, but a memory hog; Random is slow but requires little
% memory.
%-We use the exact one when the nScan is small enough; [previously, when
% we have two groups, for nScan=12,
% PiCond will initially take 384KB RAM, for nScan=14, 1.75MB, so we 
% use 12 as a cut off. (2^nScan*nScan * 8bytes/element)]. Now, if we
% assume nCond (number of groups) =3, for nScan=9, 1.35MB; for nScan=10,
% 4.51MB. If nCond = 4, for nScan=9, 18MB; for nScan=8, 4MB. 
%-If user wants all perms, then random method would seem to take an
% absurdly long time, so exact is used.

hash = 1:nScan;
hash = nCond.^(hash-1);
true_hash = iCond * hash';
% Basically the idea here is that we can use one number to uniquely
% identify one permutation. For example, if we nCond=3, the hash sequence
% will be 1, 3, 9, 27, ..., the inner product of iCond and hash will be
% unique for each iCond. This is like regard iCond as a 3-base number and
% it will be unique both on base 3 (i.e. original iCond) and base 10 (the
% inner product of iCond and hash).

if nScan<=10 | ~bAproxTst                    % exact method

    %-Generate all labellings of nScan scans as 1,2,3,...
    PiCond=uint8([]);
    for i=0:nScan-1
      tmp = uint8([]);
      for j=1:nCond
	tmp = [tmp;uint8(j*ones(nCond^i,1)),PiCond];
      end
      PiCond = tmp;
    end
    
    %-Trim to labellings with correct group numbers
    for i=1:(nCond-1)
      PiCond=PiCond(sum(PiCond'==i)==GrpCnt(i),:);
    end

    if bAproxTst                % pick random supsample of perms
	  tmp=randperm(size(PiCond,1));
	  PiCond=PiCond(tmp(1:nPiCond),:);
        % Note we may have missed iCond!  We catch this below.	
    end	

else %(nScan>10 & bAproxTst)                   % random method
    
    % Allocate final result
    PiCond = zeros(nPiCond,nScan);

    % Fill first row  
    PiCond(1,:) = iCond;
    % Fill subsequent rows, checking that we're not repeating  
    for i=2:nPiCond
      tmp=PiCond(i-1,randperm(nScan));
      while any(abs(PiCond(1:(i-1),:)*hash' - tmp*hash')/(tmp*hash') < sqrt(eps))
         tmp=PiCond(i-1,randperm(nScan));
      end
      PiCond(i,:)=tmp;
    end      
    PiCond = uint8(PiCond);
end


%-Check PiConds sum to 1*GrpCnt(1)+2*GrpCnt(2)+...
tmp = 1:nCond;
row_total = tmp*GrpCnt';

if ~all(all(double(PiCond)*ones(nScan,1)==row_total))
	error('Invalid PiCond computed!'), end

%-Find (maybe) iCond in PiCond, move iCond to 1st; 
%-----------------------------------------------------------------------
perm_hash = double(PiCond) * hash';
perm = find(abs(true_hash-perm_hash)/true_hash < sqrt(eps));

if length(perm) > 1
    perm = find(all((meshgrid(iCond,1:size(PiCond(perm,:),1))==PiCond(perm,:))'));
    if (perm~=1)
	PiCond(perm,:)=[];
	PiCond=[iCond;PiCond];
    end
elseif length(perm)==1
    %-Actual labelling must be at top of PiCond
    if (perm~=1)
	PiCond(perm,:)=[];
	PiCond=[iCond;PiCond];
    end
elseif length(perm)==0 & (nScan<=10) & bAproxTst
    % If ~bAproxTst, we won't miss iCond;
    % If (nScan>10)& bAproxTst, we use random method, iCond is guaranteed to be there.
    % So the only way in which we miss iCond is: (nScan<=10) & bAproxTst.
    % Special case where we missed iCond; order of perms is random 
    % so we can just replace first perm.
    PiCond(1,:) = iCond;
    perm = 1;
else    
    error(['Bad PiCond (' num2str(perm) ')'])
end    


%-Form non-null design matrix partitions (Globals handled later)
%=======================================================================
%-Form for HC computation at permutation perm
sHCform    = 'spm_DesMtx(double(PiCond(perm,:)),''-'',''Cond'')';
%-Condition partition
[H,Hnames] = spm_DesMtx(iCond,'-','Cond');
%-Contrast of condition effects
% (spm_DesMtx puts condition effects in index order)
if b_all_zero
    CONT       = eye(nCond);
else
    CONT       = spm_DesMtx(1:nCond, '+0')';
end
%-No block/constant
B=[]; Bnames='';

%-F statistic's numerator degree of freedom
%-Use df1 to denote it.
%Note: in snpm_cp, F statistic's denominator degree of freedom is denoted by df.
if b_all_zero
    df1       = nCond;
else
    df1       = nCond-1;
end


%-Design description
%-----------------------------------------------------------------------
%%%GrpCnt = [nScan-nFlip nFlip];
tmp = [];
for (i=1:nCond)
    tmp0 = ['%d(Grp',char(i+64),') '];
    tmp  = [tmp,tmp0];
end
tmp = ['%d Groups, between group ANOVA, 1 scan per subj: ', tmp];  
sDesign = sprintf(tmp,nCond,GrpCnt);
sPiCond = sprintf('%d permutations of conditions',size(PiCond,1));
