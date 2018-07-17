% Mfile snpm_pi_ANOVAwithinS
% SnPM PlugIn design module - Within Subject ANOVA, k diffs/contrasts per subject
% FORMAT snpm_pi_ANOVAwithinS 
%
% See body of snpm_ui for definition of PlugIn interface.
%_______________________________________________________________________
%
% snpm_pi_ANOVAwithinS is a PlugIn for the SnPM design set-up program,
% creating design and permutation matrix appropriate for one group
% analyses where there are multiple scans per subject, and where each
% scan is itself a difference image or contrast image.  This plug in
% effects a within subject ANOVA.
%
% A common use of this PlugIn is for an F test for a set of contrasts.
% For each subject, we have k contrasts jointly expressing some effect of
% interest.  Under the null hypothesis we assume that the data for each
% subject is unpeturbed by multiplication by -1.  That is, under the null
% hypothesis the multivariate measurements are all mean zero and
% symmetrically distributed.  We assume exchangeability between subjects
% (just as we usually assume independent subjects) but *do* *not* assume
% that the k values for each subject are independent.
%
% The PlugIn tests for the presence of *any* effect among the k
% contrasts.  That is, it tests the null hypothesis that all of the
% effects are mean zero.
%
%-Number of permutations
%=======================================================================
%
% There are 2^(nSubj-1) possible permutations, where nSubj is the total
% number of subjects. Intuitively, each subject can be assigned to +1 or
% -1, so we should have 2^nSubj possible permutations. However, since we
% are doing an F test and all +1's and all -1's would give us the same F
% statistic. To avoid the redundance, therefore we explicitly assign the
% first subject to +1 group.    
% 
% It is recommended that at least 7 or 8 subjects are used; with only 6
% subjects, the permutation distribution will only have 2^5 = 32 elements
% and the smallest p-value will be 1/32=0.03125.
%
%
%-Prompts
%=======================================================================
% '# subjects': Input the number of subjects.
%
% '# scans per subject': Input the number of scans per subject.
%
% 'Subject x: Select scans in time order': For each subject x, enter the
% scans to be analyzed in time order. Note: the order should be the same
% for each subject. 
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
% CONT          - a contrast for examination, a square matrix (nRepl*nRepl) 
% sDesign       - String defining the design
% sDesSave      - String of PlugIn variables to save to cfg file
%
% - design -
% H,Hnames      - Condition partition of design matrix, & effect names
% B,Bnames      - Block partition (constant term), & effect names
%
% - extra -
% iCond         - Condition indicator vector
% nRepl         - Number of scans per subject
% iRepl         - A vector indicating the order of the scans. Basically
%                 it is [1:nRepl, 1:nRepl, ...] (1:nRepl repeated by
%                 nSubj times).
% sHCform_Mtx   - A matrix that will be used when sHCform is called.  
%
%_______________________________________________________________________
% Copyright (C) 2013 The University of Warwick
% Id: snpm_pi_ANOVAwithinS.m  SnPM13 2013/10/12
% Thomas Nichols, Camille Maumet
% Based on snpm_MS1.m, V3.2 04/06/08

%-----------------------------functions-called------------------------
% spm_DesMtx
% spm_select
% spm_input
%-----------------------------functions-called------------------------

% 
% Note:  For a multisubject, no-replication design, exchagiblity is
% guaranteed for all observations by random selection of subjects from
% the populations of interest.  Hence, Xblk is all scans, and does not
% need to be accounted for.
%

%-Initialisation
%-----------------------------------------------------------------------
iGloNorm = '123';		% Allowable Global norm. codes
sDesSave = 'iRepl sHCform_Mtx';		        % PlugIn variables to save in cfg file

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
nSubj    = size(job.fsubject,2);%spm_input('# subjects','+1');
if (nSubj==1), error('SnPM:SingleSubj', 'Use single subject plug for single subjects'); end    

%-Get number of scans per subject - nSubj x nRepl design
nRepl    =  unique(arrayfun(@(x) numel(x.scans), job.fsubject));%spm_input('# scans per subject','+1');
if numel(nRepl) > 1
    error('SnPM:DifferentReplications', 'All subjects must have the same number of replications')
end


%-Get filenames and iCond, the condition labels
%=======================================================================
P     = [];
iRepl = [];
iSubj = [];
for subj=1:nSubj
    %tmp = ['Subject ',int2str(subj),': Select scans in time order'];
    P = str2mat(P, str2mat(job.fsubject(subj).scans)); %str2mat(P,spm_select(nRepl,'image',tmp));
    iRepl = [iRepl, 1:nRepl];
    iSubj = [iSubj, subj*ones(1,nRepl)];
end
P(1,:) = [];

iCond = ones(1,nSubj);

%-Get confounding covariates
%-----------------------------------------------------------------------
G = []; Gnames = ''; Gc = []; Gcnames = ''; q = nSubj*nRepl;
if numel(job.cov) > 0 %isfield(job.covariate,'cov_Val')
    for i = 1:numel(job.cov)
        d = job.cov(i).c;
        if (size(d,1) == 1)
            d = d';
        end
        nGcs = size(Gc,2);
        if size(d,1) ~= q
            error('SnPM:InvalidCovariate', sprintf('Covariate [%d,1] does not match number of scans [%d]',...
                size(job.cov(i).c,1),q))
        else
            %-Save raw covariates for printing later on
            Gc = [Gc,d];
            % Center
            d  = d - ones(q,1)*mean(d); str='';
            G = [G, d];
            dnames = job.cov(i).cname;
            Gcnames = str2mat(Gcnames,dnames);
        end
    end
    %-Strip off blank line from str2mat concatenations
    if size(Gc,2)
        Gcnames(1,:)=[];
    end
end
%-Since no FxC interactions these are the same
Gnames = Gcnames;


%-Compute permutations of subjects (we'll call them scans)
%=======================================================================
%-Compute permutations for a single exchangability block
%-----------------------------------------------------------------------
nPiCond_mx = 2^(nSubj-1);
% Note: here nPiCond is half of its usual value. The reason is we are
% calculating F stat.
nPiCond = job.nPerm;
if job.nPerm >= nPiCond_mx
    bAproxTst=0;
    if job.nPerm > nPiCond_mx
        fprintf('NOTE: %d permutations requested, only %d possible.\n',job.nPerm, nPiCond_mx)
        nPiCond = nPiCond_mx;
    end
else
    bAproxTst=1;
end

%-Two methods for computing permutations, random and exact; exact
% is efficient, but a memory hog; Random is slow but requires little
% memory.
%-We use the exact one when the nSubj is small enough; for nSubj=12,
% PiCond will initially take 384KB RAM, for nSubj=14, 1.75MB, so we 
% use 12 as a cut off. (2^nSubj*nSubj * 8bytes/element).  
%-If user wants all perms, then random method would seem to take an
% absurdly long time, so exact is used.

if nSubj<=12 || ~bAproxTst                    % exact method

    %-Generate all labellings of nSubj scans as +/- 1
    PiCond=[];
    for i=0:nSubj-2
	PiCond=[ones(2^i,1),PiCond;-ones(2^i,1),PiCond];
    end
    
    a = ones(size(PiCond,1),1);
    PiCond =[a,PiCond];
    
    if bAproxTst                 % pick random supsample of perms
	tmp=randperm(size(PiCond,1)-1);
	PiCond=PiCond([1 tmp(1:nPiCond)],:);
    end	
    
    % Set bhPerms=0. The reason is this:
    % the permutations with all +1's or all -1's will give the same F.
    % So we just want to count half of all possible permutations.
    % Another way to think about it is to always keep first subject as +1.
    bhPerms=0;
    
elseif nSubj<=53      % random method, using integer indexing
    
    d       = nPiCond-1;
    tmp     = pow2(0:nSubj-2)*iCond(1:(nSubj-1))';  % Include correctly labeled iCond

    while (d>0)
      tmp = union(tmp,floor(rand(1,d)*2^(nSubj-1)));
      tmp(tmp==2^(nSubj-1)) = [];  % This will almost never happen
      d   = nPiCond-length(tmp);
    end
    
    % randomize tmp before it is used to get PiCond
    rand_tmp=randperm(length(tmp));
    tmp=tmp(rand_tmp);
    
    PiCond = 2*rem(floor(tmp(:)*pow2(-(nSubj-1-1):0)),2)-1;
    
    a = ones(size(PiCond,1),1);
    PiCond =[a,PiCond]; 
    
    bhPerms=0;    

else    % random method, for nSubj>=54, when exceeding
        % double-precision's significand's 53 bit precision
        % For now, don't check for duplicates
    
    d       = nPiCond-1;
    PiCond  = [iCond;
	       2*(rand(nPiCond-1,nSubj)>0.5)-1];
    
    bhPerms=0;    

end

%-Find (maybe) iCond in PiCond, move iCond to 1st; negate if neccesary
%-----------------------------------------------------------------------
perm = find(all((meshgrid(iCond,1:size(PiCond,1))==PiCond)'));
if (bhPerms)
    perm=[perm,-find(all((meshgrid(iCond,1:size(PiCond,1))==-PiCond)'))];
end
if length(perm)==1
    if (perm<0), PiCond=-PiCond; perm=-perm; end
    %-Actual labelling must be at top of PiCond
    if (perm~=1)
	PiCond(perm,:)=[];
	PiCond=[iCond;PiCond];
    end
    if ~bAproxTst    
	%-Randomise order of PiConds, unless already randomized
	% Allows interim analysis	
	PiCond=[PiCond(1,:);PiCond(randperm(size(PiCond,1)-1)+1,:)];
    end	
else    
    error('SnPM:InvalidPiCond', ['Bad PiCond (' num2str(perm) ')'])
end    


%-Form non-null design matrix partitions (Globals handled later)
%=======================================================================
%-Form for HC computation at permutation perm
sHCform_Mtx = spm_DesMtx(iSubj)';

sHCform = 'diag(PiCond(perm,:)*sHCform_Mtx)*spm_DesMtx(iRepl,''-'',''Scan'')';
%-Condition partition
[H,Hnames] = spm_DesMtx(iRepl,'-','Scan');
%-Contrast of condition effects
% (spm_DesMtx puts condition effects in index order)
CONT       = eye(nRepl);
%-No block/constant
B=[]; Bnames='';

% clear iCond, because we don't want to keep it in snpm_ui.
iCond=[];

%-Calculate df1 (the numerator df for F stat)
df1 = nRepl;

%-Design description
%-----------------------------------------------------------------------
sDesign = sprintf('Multisubject, Within Subject ANOVA, multiple scans per subj: %d(subj)',nSubj);
sPiCond = sprintf('%d permutations of conditions, bhPerms=%d',size(PiCond,1)*(bhPerms+1),bhPerms);
