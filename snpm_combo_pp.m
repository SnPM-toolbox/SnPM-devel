function snpm_combo_pp(CWD, job)
% SnPM post processing and results display
% FORMAT snpm_combo_pp(CWD)
%
% CWD -	Directory containing SnPM results files
%
% If CWD is not specified then user is prompted to locate results file SnPM.mat
%
%___________________________________________________________________
% 
% snpm_combo_pp.m is the post-processing function for combined
% voxel-cluster size inference. This program is a modified version
% of snpm_pp.m used in SnPM. 
%
% Before using this program, first you should have run the "setup" and
% "compute" steps of SnPM. In the setup step, you should specify 
% to assess cluster extent. Otherwise you will not be able to 
% perform a cluster size test, and thus no combined test either.
%
% Based on SnPM_ST.mat, the supra-threshold portion of all the
% realizations, this program calculates, (a) partial test p-values
% or voxel test p-values and cluster size test p-values, then
% (b) combined functions based on p-values from (a), then (c)
% calculates p-values for combined tests, and finally (d)
% calculates the meta-combining function and perform the meta-
% combined test. This program stores information from all the
% clusters in all the permutations. So for the steps (c)-(d),
% this recorded cluster information is used without doing
% further permutations.
%
% theta = 0.5 corresponds to the equally weighted statistics,  For theta > 0.5,
% the combining function becomes more sensitive to high intensity peaks;
% whereas for theta < 0.5, the combining function becomes more sensitive to
% large clusters. For theta = 0, the test becomes a cluster size test and
% for theta = 1, the test becomes a peak intensity test.
%
% Details regarding the combined and meta-combined tests
% are found in
%
% Hayasaka S, and Nichols TE. (2004)
% Combining Voxel Intensity and Cluster Extent with 
% Permutation Test Framework.
% 
% THIS IS A BETA VERSION. THE PROGRAM IS STILL BEING TESTED!
% 
%_______________________________________________________________________
% Copyright (C) 2013 The University of Warwick
% Id: snpm_combo_pp.m  SnPM13 2013/10/12
% Satoru Hayasaka, Darren Gitelman, Camille Maumet
% Based on snpm_pp.m 

%-Setup
%=======================================================================
fprintf('\nSnPM: snpm_combo_pp\n'),fprintf('%c','='*ones(1,72)),fprintf('\n')
MLver=version;MLver=MLver(1);

%-Initialise variables & constants
%-----------------------------------------------------------------------
tol = 1e-4;	% Tolerance for comparing real numbers
		% Two reals with abs(a-b)<tol are considered equal
		% ( Reals have to be compared for equality when        )
		% ( computing adjusted p-values                        )

%-SetUp figure window
%-----------------------------------------------------------------------
Finter = spm_figure('FindWin','Interactive');
Fgraph = spm_figure('FindWin','Graphics');
if isempty(Fgraph), Fgraph=spm_figure('Create','Graphics'); end
spm_clf(Finter), spm_clf(Fgraph)
set(Finter,'Name','SnPM PostProcess');


%-Get Data
%=======================================================================
% Get analysis directory
if nargin==0
    tmp = spm_select(1,'SnPM.mat','Select SnPM.mat for analysis...');
    CWD  = spm_str_manip(tmp,'hd');
end
load(fullfile(CWD,'SnPM'))

%-Load Config file & SnPM permutation data
if exist(fullfile(CWD,'SnPMcfg.mat'),'file')
   load(fullfile(CWD,'SnPMcfg'))
else
   fprintf('Error!! Cannot find SnPMcfg.mat file.\n');
   fprintf('Make sure to run Setup and Compute in SnPM.\n');
   fprintf('Check if\n    %s \nis the right directory.\n',CWD);
   fprintf('\nEnding the program .........\n\n');
   return;
end

%-Ask whether positive or negative effects be analysed
%-----------------------------------------------------------------------
if STAT == 'T'
    bNeg = job.Tsign==-1;
else
    bNeg = 0;
end

%-Form full Tmax distribution
%-----------------------------------------------------------------------
%-Tmin are in second column of MaxT, stored with *+ve* values
if bhPerms
	MaxT   = [ MaxT; flipud(fliplr(MaxT)) ];
	PiCond = [PiCond; -flipud(PiCond)];
end
%-Take MaxT for increases or decreases according to bNeg
MaxT = MaxT(:,bNeg+1);
nPerm = size(MaxT,1);
[StMaxT, iStMaxT] = sort(MaxT);

%-Load statistic image
%-----------------------------------------------------------------------
load(fullfile(CWD,'SnPMt'))
load(fullfile(CWD,'XYZ'))

%-Negate if looking at negative contrast
%-----------------------------------------------------------------------
if bNeg
	SnPMt    = -SnPMt;
	CONT     = -CONT;
end

%-Get ORIGIN, etc
DIM    = [V(1).dim(1) V(1).dim(2) V(1).dim(3)];
VOX    = [V(1).mat(1,1) V(1).mat(2,2) V(1).mat(3,3)];
MAT    = V(1).mat;
IMAT   = inv(MAT);
ORIGIN = IMAT(1:3,4);

% Template vol structure
Vs0 =  V(1);
Vs0.dt = [spm_type('float64'), spm_platform('bigend')];
% Vs0 = struct('fname',	'',...
% 	     'dim',	[DIM,spm_type('float')],...
% 	     'mat',	MAT,...
% 	     'pinfo',	[1 0 0]',...
% 	     'descrip',	'');

%-Write out filtered statistic image?  (Get's done later)
%-----------------------------------------------------------------------
WrtFlt = isfield(job.WriteFiltImg, 'name'); %spm_input('Write filtered statistic img?','+1','y/n',[1,0],2);
if WrtFlt
	%WrtFltFn = 'SnPMt_filtered';
	%WrtFltFn=spm_input('Filename ?','+1','s',WrtFltFn);
    %    WrtFltFn = [WrtFltFn, '.img'];
    WrtFltFn = job.WriteFiltImg.name;
end



%-Get inference parameters
%=======================================================================

%-Get corrected threshold
%-----------------------------------------------------------------------
alpha = job.Thr.Clus.ClusMass.PFilt; %spm_input('Corrected p value for filtering','+1','e',0.05);

%-Compute critical threshold for level alpha test
%-----------------------------------------------------------------------
if alpha < 1
	c=ceil((1-alpha)*nPerm);
	C_MaxT=StMaxT(c);
else
	%-Just use voxels with +ve valued SnPMt
	C_MaxT=0;
end

%-Ask whether SupraThreshold cluster size test required
%-----------------------------------------------------------------------
%-If chosen alpha specifies a critical threshold less than the threshold
% ST_Ut used to collect suprathreshold data in snpm_cp, then it makes
% no sense to analyse by spatial extent since the voxels are individually
% significant.
%-bST flags whether spatial extent information was collected.)
bSpatEx = bST & exist(fullfile(CWD,'SnPM_ST.mat'))==2;
if bSpatEx && (C_MaxT <= ST_Ut) && (alpha ~= 1)
  str = 'Voxelwise corrected threshold = %g, which is smaller ';
  str = [str 'than minimum saved suprathreshold information (%g)'];
  str = [str '\nAll results significant voxelwise.'];
  warning('SnPM:VoxelwiseCorrThreshSmaller', sprintf(str,C_MaxT,ST_Ut))
end
if bSpatEx~=1
  error('SnPM:SupraStatMissing', 'Suprathreshold stats not collected!  Cannot do cluster-combining!')
end

%-Get primary threshold for STC analysis if requested
%-----------------------------------------------------------------------
if bSpatEx
    % Save original ST_Ut
    ST_Ut_0 = ST_Ut;
    %-Threshold must be greater or equal to that (ST_Ut) used to collect
    % suprathreshold data in snpm_cp
    %-If a test level alpha has been set, then it there's no sense in having
    % the threshold greater than C_MaxT, above which voxels are individually 
    % significant
    %tmp = 0;
    primaryThresh = job.Thr.Clus.ClusMass.PrimThresh;
    if bVarSm
        %-If using pseudo-statistics then can't use (uncorrected) 
        % upper tail p-values to specify primary threshold
        if alpha == 1	% Not filtering on significance
            if ~(primaryThresh>=ST_Ut)
                error('SnPM:PseudoStatWithPValueThreshold', ...
                    ['Using pseudo-statistics you can''t use (uncorrected)'... 
                        'upper tail p-values to specify primary threshold']);
            end   
        else
            if ~(primaryThresh>=ST_Ut && primaryThresh<C_MaxT)
                        error('SnPM:PseudoStatWithPValueThreshold', ...
                            ['Using pseudo-statistics you can''t use (uncorrected)'... 
                                'upper tail p-values to specify primary threshold']);
            end  
        end
    else
        %-Statistic image is t with df degrees of freedom
        pU_ST_Ut  = 1-spm_Tcdf(ST_Ut,df);
        if alpha==1	% Not filtering on significance
            if ~( primaryThresh>=ST_Ut || (primaryThresh>0 && primaryThresh<=pU_ST_Ut))
                error('SnPM:InvalidPrimaryThresh', ['Primary threshold must be >=' num2str(ST_Ut) ...
                    ' and >0 and <=' num2str(pU_ST_Ut) ]);
            end  
        else
            pU_C_MaxT = 1-spm_Tcdf(C_MaxT,df);
            if ~((primaryThresh>=ST_Ut && primaryThresh<C_MaxT) || ...
                    (primaryThresh>pU_C_MaxT && primaryThresh<=pU_ST_Ut))
                error('SnPM:InvalidPrimaryThresh', ['Primary threshold must be >=' num2str(ST_Ut) ...
                    ' and <' num2str(C_MaxT) ' or >' num2str(pU_C_MaxT) ...
                    ' and <= ' num2str(pU_ST_Ut)]);
            end
            clear pU_C_MaxT
        end
        clear pU_ST_Ut
        if (primaryThresh < 1) 
            primaryThresh = spm_invTcdf(1-primaryThresh,df); 
        end
    end
    ST_Ut = primaryThresh;

    %
    % Getting combined test parameters
    %
    Theta = job.Thr.Clus.ClusMass.Theta;
    if (Theta<0+tol || Theta>1-tol)
        error('SnPM:InvalidTheta', 'Theta should be between 0 and 1');
    end
    mTheta = Theta/(1-Theta);  %-Weight for mass combining
    
    %-Picking which combining test
    % Only mass combining for now
    iW    = 3;%    = spm_input('Choose combining function',1,'b', ...
%                       'Fisher|Tippet|Mass|All',[1:4],1);

end


%-Show permutation distributions?
%-----------------------------------------------------------------------
%ShwDst = spm_input('Display permutation distribution[s]?','+1','y/n',[1,0],1);
ShwDst = 1;

%=======================================================================
%- C O M P U T A T I O N
%=======================================================================
set(Finter,'Pointer','Watch')

%-Calculate distribution of Maximum Suprathreshold Cluster size
%-Calculate critical Suprathreshold Cluster Size
%=======================================================================
if bSpatEx
	fprintf('Working on spatial extent...\n');

	%-Compute suprathreshold voxels - check there are some
	%---------------------------------------------------------------
	fprintf('\tComputing suprathreshold voxels...');
	Q     = find(SnPMt > ST_Ut);
	SnPMt = SnPMt(Q);
	XYZ   = XYZ(:,Q);
	if isempty(Q)
		set(Finter,'Pointer','Arrow')
		figure(Fgraph)
		axis off
		text(0,0.97,CWD,'Fontsize',16,'FontWeight','Bold');
		text(0,0.93,sprintf('No voxels above threshold %4.2f',ST_Ut));
		ShowDist(MaxT,C_MaxT);
		return
	end
	fprintf('done\n')

	%-Load & condition statistics
	%---------------------------------------------------------------
	fprintf('\tLoading & conditioning SupraThreshold statistics...');
	load(fullfile(CWD,'SnPM_ST'))
	%-SnPM_ST stores columns of [x;y;z;abs(t);perm] with perm negative
	% where the exceedence was t < -ST_Ut_0
	%-Trim statistics according to threshold ST_Ut, if ST_Ut > ST_Ut_0
	tmp = find(SnPM_ST(4,:)>ST_Ut);
	SnPM_ST = SnPM_ST(:,tmp);
	clear tmp;	    

	%-Negate perm numbers if looking at negative contrast
	if bNeg
		SnPM_ST(5,:) = -SnPM_ST(5,:);
	end
	if bhPerms
		%-Renumber negative perms according to -flipud PiCond
		tQ = SnPM_ST(5,:)<0;
		SnPM_ST(5,tQ) = nPerm +1 +SnPM_ST(5,tQ);
	else
		%-Not bhPerms: Lose entries for negative excursions
		SnPM_ST = SnPM_ST(:,SnPM_ST(5,:)>0);
	end
	fprintf('done\n')

	%-Calculate distribution of Maximum SupraThreshold Cluster size
	%---------------------------------------------------------------
        % ClInfo = [x;y;z;ClusterPeak;ClusterSize;ClusterMass;ClusterInd;Perm]
        ClInfo = [];

	fprintf('\tComputing dist. of max SupraThreshold cluster size: ');
	MaxSTCS = zeros(nPerm,1);
	SetLvl  = zeros(nPerm,1);
	fprintf('\nPerms left:     ');
    
	% preallocate ClInfo to improve speed. This will be far larger than it
	% needs to be but it should be OK.
	ClInfo = zeros(8,size(SnPM_ST,2));
	ClInfoEnd = 0;
	for i = nPerm:-1:1
        if (rem(i,10)==0)
		  fprintf('\b\b\b\b%-4u',i)
		  drawnow
		end


		tQ = (SnPM_ST(5,:)==i);
		if any(tQ)
                        %-Filter out data for this perm
                        subSnPM_ST = SnPM_ST(:,tQ);

         		%-Compute cluster labellings for this perm
			%===== SnPM99 change =================
         		Locs_mm = SnPM_ST(1:3,tQ);
			Locs_mm (4,:) = 1;
			Locs_vox = IMAT * Locs_mm;
			%===== SnPM99 change =================
			tmp = spm_clusters(Locs_vox(1:3,:));
			%-Work out maximum cluster size (honest!)
			SetLvl(i)  = max(tmp);
                        tmpCS = zeros(1,SetLvl(i));
                        tmpCS = diff(find([diff([0,sort(tmp)]),1]));
                        MaxSTCS(i) = max(tmpCS);
                        for j=1:SetLvl(i)
	                   %-Recording peak, its location, and excess mass
                           % then creates a big vector ClInfo with 
                           % all cluster info
                           tmpCLoc  = (tmp==j);
                           subsubST = subSnPM_ST(:,tmpCLoc);
                           [tmpPeak, tmpPLoc] = max(subsubST(4,:));
                           tmpXYZ   = subsubST(1:3,tmpPLoc);
                           tmpMass  = sum((subsubST(4,:)-ST_Ut).^mTheta);
                           tmpVec   = [tmpXYZ; tmpPeak; tmpCS(j); tmpMass; j; i];
                           ClInfo(:,ClInfoEnd+j) = tmpVec;
                        end
			ClInfoEnd = ClInfoEnd+SetLvl(i);
		end
	end
    
	% get rid of any unused entries in ClInfo
	ClInfo = ClInfo(:,1:ClInfoEnd);
    
        clear tmpCS tmpCLoc tmpPeak tmpPLoc tmpXYZ tmpMass tmpVec;

	fprintf('\b\b\b\bdone\n');
	%-Save perm 1 stats for use later - [X;Y;Z;T;perm;STCno]
	STCstats = [ SnPM_ST(:,tQ); tmp];

	%-Compute critical SupraThreshold Cluster size
	[StMaxSTCS, iStMaxSTCS] = sort(MaxSTCS);
	if alpha < 1
		C_STCS = StMaxSTCS(c);
	else
		C_STCS = 0;
	end

	%-Check XYZ for points > ST_Ut in perm 1 matches
	% XYZ computed above for SnPMt > ST_Ut
	if ~all(all( SnPM_ST(1:3,SnPM_ST(5,:)==1) == XYZ ))
		error('SnPM:InvalidSTXYZ', 'ST XYZ don''t match between STCS & thresh')
	end

        %-- VOXEL CLUSTER COMBINED TEST STATISTICS ----
        %-- Assigning each voxel corrected voxel and cluster p-values
        %   then calculate the combining functions
        % CorrPs = [Vox p; Cl p; Fisher p; Tippet p; ExcessMass p; Meta p]
        % ComboF = [Fisher Wf; Tippet Wt; Excess mass Wm; Meta Wa]
        CorrPs = ones(6,size(ClInfo,2));
        ComboF = zeros(4,size(ClInfo,2));

	fprintf('Calculating corr p-values for each voxel / cluster\n');
	for ip = 1:(nPerm-1)
        if rem(ip,50)==0, fprintf('.'), end
        Qvox = find(ClInfo(4,:)>StMaxT(ip));
        CorrPs(1,Qvox) = (nPerm-ip)/nPerm;
        Qcl  = find(ClInfo(5,:)>StMaxSTCS(ip));
        CorrPs(2,Qcl) = (nPerm-ip)/nPerm;
    end
	fprintf('.Done!\n');

	fprintf('Calculating voxel-cluster combining functions\n');
        ComboF(1,:) = -2*(2*Theta*log(CorrPs(1,:)) + ...
                          2*(1-Theta)*log(CorrPs(2,:)));
        ComboF(2,:) = 1 - min(2*Theta*log(CorrPs(1,:)), ...
                          2*(1-Theta)*log(CorrPs(2,:)));
	ComboF(3,:) = ClInfo(6,:);
	fprintf('.Done!\n');

        %-- Corrected p-value for combo function
        MaxWf = zeros(nPerm,1);
        MaxWt = zeros(nPerm,1);
        MaxWm = zeros(nPerm,1);
        for i=1:nPerm
           QPerm  = find(ClInfo(8,:) == i);
           if ~isempty(QPerm)
              tmpWf    = ComboF(1,QPerm);
              MaxWf(i) = max(tmpWf);
              tmpWt    = ComboF(2,QPerm);
              MaxWt(i) = max(tmpWt);
              tmpWm    = ComboF(3,QPerm);
              MaxWm(i) = max(tmpWm);
           end
        end

        [StMaxWf, iStMaxWf] = sort(MaxWf);
        [StMaxWt, iStMaxWt] = sort(MaxWt);
        [StMaxWm, iStMaxWm] = sort(MaxWm);
	C_Wcomb = zeros(4,1);
	if alpha < 1
		C_Wcomb(1) = StMaxWf(c);
		C_Wcomb(2) = StMaxWt(c);
		C_Wcomb(3) = StMaxWm(c);
	else
		C_Wcomb(1) = 0;
		C_Wcomb(2) = 0;
		C_Wcomb(3) = 0;
	end

	fprintf('Calculating corr p-values for combining functions\n');
	for ip = 1:(nPerm-1)
           if rem(ip,50)==0, fprintf('.'), end
           Qcomb = find(ComboF(1,:)>StMaxWf(ip));
           CorrPs(3,Qcomb) = (nPerm-ip)/nPerm;
           Qcomb = find(ComboF(2,:)>StMaxWt(ip));
           CorrPs(4,Qcomb) = (nPerm-ip)/nPerm;
           Qcomb = find(ComboF(3,:)>StMaxWm(ip));
           CorrPs(5,Qcomb) = (nPerm-ip)/nPerm;
        end
	fprintf('.Done!\n');


        %
        % Doing Meta-combining
        % 
 
	%-First, calculate the meta statistic
        fprintf('Calculating meta-combining function\n');
        ComboF(4,:) = 1 - min([log(CorrPs(3,:))' ...
                      log(CorrPs(4,:))' log(CorrPs(5,:))'],[],2)';
  

        % Then calculate the max meta-combining function
        MaxWa = zeros(nPerm,1);
        for ip=1:nPerm
           QPerm  = find(ClInfo(8,:) == ip);
           if ~isempty(QPerm)
              tmpWa    = ComboF(4,QPerm);
              MaxWa(ip) = max(tmpWa);
           end
        end


        % Corrected critical meta-statistic
        [StMaxWa, iStMaxWa] = sort(MaxWa);
	if alpha < 1
		C_Wcomb(4) = StMaxWa(c);
	else
		C_Wcomb(4) = 0;
	end

        % P-value calculation for meta-combining
	fprintf('Calculating corr p-values for meta-combining \n');
	for ip = 1:(nPerm-1)
           if rem(ip,50)==0, fprintf('.'), end
           Qcomb = find(ComboF(4,:)>StMaxWa(ip));
           CorrPs(6,Qcomb) = (nPerm-ip)/nPerm;
        end
	fprintf('.Done!\n');

        % Combining all the MaxW distributions to be plotted later
        MaxW  = [MaxWf; MaxWt; MaxWm; MaxWa];
       
end


%-Save some time consuming results
%-----------------------------------------------------------------------
if bSpatEx, save SnPM_pp.mat STCstats MaxSTCS SetLvl, end
if bSpatEx, save SnPM_combo.mat CorrPs ComboF ClInfo MaxWf MaxWt MaxWm MaxWa, end

%-Filter data at specified corrected p-value alpha
%=======================================================================
if bSpatEx
    %-Analysing spatial extent
    %-NB:alpha==1 implies C_MaxT==C_STCS==0.
    % Since ST_Ut>0 filtering has no effect if alpha==1, so skip it.
    if alpha<1
	%-Filter on significance of cluster size
	%---------------------------------------------------------------
	fprintf('Filtering on cor.sig. at suprathreshold cluster level...');
	nSTC     = max(STCstats(6,:));
	STCS     = diff(find([diff([0,sort(STCstats(6,:))]),1]));
	Q        = [];
	for i = 1:nSTC
		tQ = find(STCstats(6,:)==i);
	     ttQ= find(ClInfo(7,:)==i & ClInfo(8,:)==1);
		if ( STCS(i) > C_STCS || max(STCstats(4,tQ)) > C_MaxT || ...
                     max(ComboF(iW,ttQ)) > C_Wcomb(iW))
			Q        = [Q tQ];
		end
	end
	if ~isempty(Q)
		SnPMt    = SnPMt(Q);
		XYZ      = XYZ(:,Q);
		STCstats = STCstats(:,Q);
	end
	fprintf('done\n')
    end
else
	%-Truncate at critical threshold for level alpha test
	% NB if alpha==1 then C_MaxT is set to 0, and filter on +ve SnPMt
	fprintf('Filtering on cor.sig. at voxel level...');
	Q     = find(SnPMt > C_MaxT);
	if length(Q)
		SnPMt = SnPMt(Q);
		XYZ   = XYZ(:,Q);
	end
	fprintf('done\n')
end



%-Return if there are no voxels
%-----------------------------------------------------------------------
if isempty(Q)
	set(Finter,'Pointer','Arrow')
	figure(Fgraph)
	axis off
	text(0,0.97,CWD,'Fontsize',16,'FontWeight','Bold');
	tmp='voxels'; if bSpatEx, tmp='suprathreshold clusters'; end
	text(0,0.93,sprintf(...
		'No %s significant at alpha=%6.4f (corrected)',tmp,alpha));
	if bSpatEx
	  ShowDist(MaxT,C_MaxT,MaxSTCS,C_STCS);
	else	   
	  ShowDist(MaxT,C_MaxT);
	end
	return
end

%-Characterize local excursions in terms of maxima:
% #voxels STC_N; MaxTs STC_SnPMt; locations STC_XYZ, & region# STC_r
%-----------------------------------------------------------------------
%===== SnPM99 change =============================================
TempXYZmm = XYZ;
TempXYZmm(4,:) = 1;
TempXYZvoxel = IMAT*TempXYZmm;
TempXYZvoxel= TempXYZvoxel(1:3,:);

[STC_N, STC_SnPMt, STC_XYZ, STC_r] = spm_max(SnPMt,TempXYZvoxel);

TempXYZvoxel = STC_XYZ;
TempXYZvoxel(4,:) = 1;
TempXYZmm = MAT * TempXYZvoxel;
STC_XYZ = TempXYZmm(1:3,:);
%===== SnPM99 change =============================================

%-Compute adjusted significances for local maxima, & regions (if required)
%-----------------------------------------------------------------------
Pt = ones(size(STC_r));
for i = 1:length(STC_r)
	%-Use a > b -tol rather than a >= b to avoid comparing reals
	Pt(i) = sum(MaxT > STC_SnPMt(i) -tol) / nPerm;
end
%if ~bVarSm
%	Pu    = 1 - spm_Tcdf(STC_SnPMt,df);
% end
if bSpatEx
	%-Compute single step adjusted p-values for region size: pSTSC_SS
	Pn    = ones(1,length(STC_r));
	Ww    = ones(4,length(STC_r));
	Pw    = ones(4,length(STC_r));
	for i = 1:length(STC_r)
	   Pn(i) = sum(MaxSTCS>=STC_N(i)) / nPerm;
           Qvox  = find((ClInfo(1,:)==STC_XYZ(1,i)) & ...
                        (ClInfo(2,:)==STC_XYZ(2,i)) & ...
                        (ClInfo(3,:)==STC_XYZ(3,i)) & ...
                        (ClInfo(8,:)==1));
	   if ~isempty(Qvox)
	      Ww(:,i) = ComboF(:,Qvox);
              Pw(:,i) = CorrPs(3:6,Qvox);
           end
	end
	save('SnPM_combo','Ww','Pw','STC_XYZ','-append');
end




%=======================================================================
%-D I S P L A Y
%=======================================================================
figure(Fgraph)


if (ShwDst)
  axis off
  if (bSpatEx)
    text(0,0.97,'Permutation Distribution','Fontsize',16,'FontWeight','Bold');
    ShowDist(MaxT,C_MaxT,MaxSTCS,C_STCS);
  else	     
    text(0,0.97,'Permutation Distributions','Fontsize',16,'FontWeight','Bold');
    ShowDist(MaxT,C_MaxT);
  end
end

% Change the style of pressing ' return' to clicking on a new button.
%spm_print
%disp('Press <RETURN> to continue'); pause
if false
    if spm_input('Review permutation distributions.',1,'bd',...
                  'Print & Continue|Continue',[1,0],1)
        spm_print
    end
end

spm_clf(Fgraph)


%-Maximium intenisty projection of SPM{Z}
%=======================================================================
hmip = axes('Position',[0.05 0.5 0.5 0.5]);
snpm_mip(SnPMt,XYZ,MAT,DIM); axis image
if bVarSm
	title('SnPM{Pseudo-t}','FontSize',16,'Fontweight','Bold')
else
	title('SnPM{t}','FontSize',16,'Fontweight','Bold')
end

%-Design matrix and contrast
%=======================================================================
hDesMtx = axes('Position',[0.65 0.6 0.2 0.2]);
imagesc((spm_DesMtx('Sca', [H,C,B,G],HCBGnames) + 1)*32)
xlabel 'Design Matrix'
set(hDesMtx,'XTick',[],'XTickLabel','')
hConAxes = axes('Position',[0.65 0.8 0.2 0.1]);

h = bar(CONT(1,:), 'FaceColor',[1 1 1]*.8, 'BarWidth', 1);  
hold on
tX = get(h,'XData'); tY = get(h,'YData');
bar_width = get(h, 'BarWidth');
set(gca,'Xlim',[min(tX(:))-bar_width/2 max(tX(:))+bar_width/2]) 
title 'contrast'; axis off; hold off




%-Table of regional effects
%=======================================================================
%-Table headings
%-----------------------------------------------------------------------
hTable = axes('Position',[0.1 0.1 0.8 0.46],...
	      'YLim',[0,27],'YLimMode','manual',...
	      'DefaultTextInterpreter','Tex',...
	      'DefaultTextVerticalAlignment','Baseline',...
	      'Visible','off');
%	      'DefaultHorizontalAlignment','right',...
y = 26;
dy=1;
text(0,y,['P values & statistics:   ',spm_str_manip(CWD,'a40')],...
	'FontSize',12,'FontWeight','Bold','Interpreter','none');
y  = y -dy;
line([0 1],[y y],'LineWidth',3,'Color','r')
y = y -dy;

tCol       = [  0.00      0.16  ...			%-Cluster
	        0.25      0.43  ...                     %-combo cluster
		0.55      0.71   ...                    %-Voxel
                0.86      0.93    1.00];        	%-XYZ

PF    = spm_platform('fonts');   %-Font names (for this platform)

%-Construct table header
%-----------------------------------------------------------------------
set(gca,'DefaultTextFontName',PF.helvetica,'DefaultTextFontSize',10)

Hp = [];
h  = text(0.10,y,	'cluster-level','FontSize',10,'HorizontalAlignment','Center' );		
h  = line([tCol(1),0.20],[1,1]*(y-dy/4),'LineWidth',0.5,'Color','r');	
h  = text(tCol(1),y-9*dy/8,	'\itp_{corrected}');        Hp = [Hp,h];
h  = text(tCol(2),y-9*dy/8,	'\itk ');		


h  = text(0.37,y,	'combo cluster-level','FontSize',10, 'HorizontalAlignment','Center');   Hp=[Hp,h];		
h  = line([tCol(3),0.50],[1,1]*(y-dy/4),'LineWidth',0.5,'Color','r');	Hp=[Hp,h];		
h  = text(tCol(3),y-9*dy/8,	'\itp_{corrected}');        Hp = [Hp,h];
h  = text(tCol(4),y-9*dy/8,	'\itw ');	Hp=[Hp,h];	

								
text(0.67,y,		'voxel-level','FontSize',10, 'HorizontalAlignment','Center');
line([tCol(5),0.80],[1,1]*(y-dy/4),'LineWidth',0.5,'Color','r');
h  = text(tCol(5),y-9*dy/8,	'\itp_{FWE-corr}');	
if ~bVarSm
	h  = text(tCol(6),y-9*dy/8,	'\itt');
else
	h  = text(tCol(6)-0.02,y-9*dy/8,	'Pseudo-t');
end


text(tCol(7),y-dy/2,'{x,y,z} mm','FontSize',10);

y     = y - 7*dy/4;
line([0 1],[y y],'LineWidth',1,'Color','r')
y     = y - 5*dy/4;


%    text(0.00,y,'region','FontSize',10);
%text(0.00,y,'size {k}','FontSize',10);
%if bSpatEx
%	text(0.11,y,'P(K_{max}>= k)','FontSize',10);
%	text(0.62,y,'w','FontSize',10);
%	text(0.68,y,'P(W_{max}>= w)','FontSize',10);
%	
%end
%
%text(0.35,y,'P(T_{max}>= u)','FontSize',10);
%
%text(0.84,y,'{x,y,z} mm','FontSize',10);
%y = y -0.8;
%line([0 1],[y y],'LineWidth',3,'Color',[0 0 0])
%y  = y -1;

Fmtst = {	'%0.4f','%0.0f', ...	        %-cluster
                '%0.4f','%6.2f', ...            %-combo cluster
                '%0.4f','%6.2f', ...            %-voxel 	
		'%3.0f','%3.0f', '%3.0f'};      %-XYZ

%-Column Locations
%-----------------------------------------------------------------------
%tCol       = [  0.08     0.18     0.33  ...			
%	        0.44     0.66     0.76           ...  
%                0.86     0.93     1.00];			%-XYZ

%-List of maxima
%-----------------------------------------------------------------------
r = 1;
bUsed = zeros(size(STC_SnPMt));
while max(STC_SnPMt.*(~bUsed)) && (y > 3)

	[null, i] = max(STC_SnPMt.*(~bUsed));	% Largest t value
	j         = find(STC_r == STC_r(i));	% Maxima in same region


	%-Print region and largest maximum
	%-------------------------------------------------------------------
	
	StrAttr = {'Fontsize',10,'ButtonDownFcn','get(gcbo, ''UserData'')',...
		  'HorizontalAlignment','right'};
	StrAttrB = {StrAttr{:},'FontWeight','Bold'};
%	text(0.00,y,sprintf('%0.0f',r),'UserData',r,StrAttrB{:})
	
        if bSpatEx
		text(tCol(1)+0.09,y,sprintf(Fmtst{1},Pn(i)),...
				'UserData',Pn(i),StrAttrB{:})
                text(tCol(3)+0.08,y,sprintf(Fmtst{3},Pw(iW,i)),...
				'UserData',Pw(iW,i),StrAttrB{:})
                text(tCol(4)+0.04,y,sprintf(Fmtst{4},Ww(iW,i)),...
				'UserData',Ww(iW,i),StrAttrB{:})
	else
	  set(Hp, 'Visible','off')
	end


        text(tCol(2)+0.03,y,sprintf(Fmtst{2},STC_N(i)),'UserData',STC_N(i),StrAttrB{:})
	text(tCol(5)+0.08,y,sprintf(Fmtst{5},Pt(i)),...
	                       'UserData',Pt(i),StrAttrB{:})
	text(tCol(6)+0.04,y,sprintf(Fmtst{6},STC_SnPMt(i)),...
	                       'UserData',STC_SnPMt(i),StrAttrB{:})
	
	text(tCol(7),y,sprintf(Fmtst{7},STC_XYZ(1,i)),...
	                       'UserData',STC_XYZ(:,i),StrAttrB{:})
	text(tCol(8),y,sprintf(Fmtst{8},STC_XYZ(2,i)),...
	                       'UserData',STC_XYZ(:,i),StrAttrB{:})
	text(tCol(9),y,sprintf(Fmtst{9},STC_XYZ(3,i)),...
	                       'UserData',STC_XYZ(:,i),StrAttrB{:})
	y = y -1;

	%-Print up to 3 secondary maxima (>8mm apart)
	%-------------------------------------------------------------------
	[null, k] = sort(-STC_SnPMt(j));	% Sort on t value
	D         = i;
	for i = 1:length(k)
	    d     = j(k(i));
	    if min( sqrt( sum((STC_XYZ(:,D) - ...
			STC_XYZ(:,d)*ones(1,size(D,2))).^2) ) ) > 8
		if length(D) < 3
       	
	text(tCol(5)+0.08,y,sprintf(Fmtst{5},Pt(d)),...
	                       'UserData',Pt(d),StrAttr{:})
	text(tCol(6)+0.04,y,sprintf(Fmtst{6},STC_SnPMt(d)),...
	                       'UserData',STC_SnPMt(d),StrAttr{:})
	
	
	text(tCol(7),y,sprintf(Fmtst{7},STC_XYZ(1,d)),...
	                       'UserData',STC_XYZ(:,d),StrAttr{:})
	text(tCol(8),y,sprintf(Fmtst{8},STC_XYZ(2,d)),...
	                       'UserData',STC_XYZ(:,d),StrAttr{:})
	text(tCol(9),y,sprintf(Fmtst{9},STC_XYZ(3,d)),...
	                       'UserData',STC_XYZ(:,d),StrAttr{:})
		  
		    D = [D d];
		    y = y -1;
		end
	    end
	end

	bUsed(j) = (bUsed(j) | 1 );		%-Mark maxima as "used"
	r = r + 1;				% Next region
end
clear i j k D d r


%-Footnote with SnPM parameters
%=======================================================================
line([0,1],[0.5,0.5],'LineWidth',1,'Color','r')
y = 0;
if bSpatEx
    tmp = sprintf('Threshold = %7.4f',ST_Ut);
    if ~bVarSm
	tmp=[tmp,sprintf(' (p = %6.4f)',spm_Tcdf(-ST_Ut,df))];
    end
    text(0,y,tmp,'FontSize',8)
    text(0.7,y,sprintf('Critical STCS = %d voxels',C_STCS),'FontSize',8)
    y = y -0.8;
end
text(0,y,sprintf('alpha = %6.4f, df = %d',alpha,df),'FontSize',8)
text(0.7,y,sprintf('Critical threshold = %7.4f',C_MaxT),'FontSize',8)
y = y -0.8;
text(0,y,sprintf('Volume = %d %5.2fx%5.2fx%5.2f mm voxels',...
	S,VOX(1),VOX(2),VOX(3)),'FontSize',8);
y = y -0.8;
text(0,y,sprintf('Design: %s',sDesign),'FontSize',8);
y = y -0.8;
text(0,y,sprintf('Perms: %s',sPiCond),'FontSize',8);
if bVarSm
	y = y -0.8;
	text(0,y,sVarSm,'FontSize',8)
end

y = -1.6;
NameW = {'Fisher','Tippet','Excess mass','Meta-combining'};
if bSpatEx
   text(0.7,y,sprintf('Combining function W: %s',NameW{iW}),'FontSize',8);
   y = y - 0.8;
   text(0.7,y,sprintf('Critical W = %7.4f',C_Wcomb(iW)),'FontSize',8);
   y = y - 0.8;
   text(0.7,y,sprintf('Theta = %4.2f',Theta),'FontSize',8);
end

%spm_print
%Set a button, so the user can decide whether to print the page of results to spm2.ps.
if false % interactive display inactive
    if spm_input('Review results.',1,'bd','Print|Done',[1,0],1)
        spm_print
    end 
end

set(Finter,'Pointer','Arrow')

%- Image output?
%=======================================================================
%-Write out filtered SnPMt?
if WrtFlt

	Fname = WrtFltFn;

	%-Dont ask about t2z conversion
	%---------------------------------------------------------------
	bt2z = 0;
	if ~bVarSm
% 	    bt2z = spm_input('Convert t -> z prior to writing?',...
% 	    	'+1','y/n')=='y';
	    tmp = sprintf('SPMt - %d df',df);
	else
	    tmp = 'SnPMt - pseudo t';
	end

	%-Reconstruct filtered image from XYZ & SnPMt
	%---------------------------------------------------------------
	t = zeros(1,prod(DIM));
	if ~bt2z
		t(spm_xyz2e(XYZ,V)) = SnPMt;
	else
		t(spm_xyz2e(XYZ,V)) = spm_t2z(SnPMt,df);
		tmp = [tmp,' (Gaussianised)'];
	end
	if ~bSpatEx
		tmp=sprintf('%s p<%10g corrected @ voxel level',tmp,alpha);
	elseif bt2z
		tmp=sprintf('%s p<%10g corrected @ cluster level, u=%4.2',...
	    		tmp,alpha,spm_t2z(ST_Ut,df));
	else
		tmp=sprintf('%s p<%10g corrected @ cluster level, u=%4.2',...
	    		tmp,alpha,ST_Ut,df);
	end

	%-Write out to analyze file
	%---------------------------------------------------------------
% 	Vs = Vs0; 
% 	Vs.fname = Fname; Vs.descrip = tmp;
% 	Vs.dim   = Vs.dim(1:3);
% 	Vs = sf_create_vol(Vs);
% 	t = reshape(t,DIM);
% 	for p=1:Vs.dim(3)
% 	  Vs = spm_write_plane(Vs,t(:,:,p),p);
% 	end
% 	Vs = sf_close_vol(Vs);
% 	clear t
    %-Write out to image file
    %---------------------------------------------------------------
    Vs = snpm_clone_vol(Vs0, Fname, tmp); 
    Vs =  spm_create_vol(Vs);
    t = reshape(t,DIM);
    for p=1:Vs.dim(3)
        Vs = spm_write_plane(Vs,t(:,:,p),p);
    end
    Vs =  sf_close_vol(Vs);
    clear t
end

%-Reset Interactive Window
%-----------------------------------------------------------------------
spm_figure('Clear','Interactive')



function ShowDist(T,cT,C,cC)
%
% Display permutation distributions on current figure, using position
% pos.
%
% We assume that the observed (aka correctly labeled data) are the first
% element in the distribution.
%

pos1 = [0.125 0.50 0.75 0.3];
pos2 = [0.125 0.08 0.75 0.3];

%
% Display intensity perm dist
%
% Bin width rule from Scott, "Multivariate Density Estimation", 1992, pg 55.
%
axes('position',pos1)
BinWd  = 3.5*std(T)*length(T)^(-1/3);  
nBin   = floor((max(T)-min(T))/BinWd)+1;
nBin   = min(max(nBin,10),50);

hist(T,nBin);
h = findobj(gca,'Type','patch');
set(h,'FaceColor',[.5 .5 .5]);
title('Permutation Distribution:  Maximum Statistic','FontSize',14)
Ylim = get(gca,'ylim'); Xlim = get(gca,'xlim');
line(T(1)*[1 1],Ylim.*[1 0.95],'LineStyle',':');
text(T(1)+diff(Xlim)*0.01,Ylim(2)*0.95,'Observed','FontSize',10)
line(cT*[1 1],Ylim.*[1 0.85],'LineStyle','-');
text(cT+diff(Xlim)*0.01,Ylim(2)*0.85,'Threshold','FontSize',10);

if (nargin>2)
    
    %
    % Display cluster size perm dist
    %
    
    axes('position',pos2)
    BinWd  = 3.5*std(C)*length(C)^(-1/3);
    nBin   = floor((max(C)-min(C))/BinWd)+1;
    nBin   = min(max(nBin,10),50);
    
    hist(C,nBin);
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor',[.5 .5 .5]);
    title('Permutation Distribution: Maximum Cluster Size','FontSize',14)
    Ylim = get(gca,'ylim'); Xlim = get(gca,'xlim');
    line(C(1)*[1 1],Ylim.*[1 0.95],'LineStyle',':');
    text(C(1)+diff(Xlim)*0.01,Ylim(2)*0.95,'Observed','FontSize',10)
    line(cC*[1 1],Ylim.*[1 0.85],'LineStyle','-');
    text(cC+diff(Xlim)*0.01,Ylim(2)*0.85,'Threshold','FontSize',10)
    
end

return



function Vs = sf_create_vol(Vs)
switch spm('ver')
 case 'SPM99'
  Vs = spm_create_image(Vs);
 case 'SPM2'
  Vs = spm_create_vol(Vs);
 case 'SPM5'
  Vs = spm_create_vol(Vs);
end         
return


function Vs = sf_close_vol(Vs)
switch spm('ver')
 case 'SPM99'
  % n/a; File closed by spm_write_plane
 case 'SPM2'
  % Vs = spm_close_vol(Vs);
  % Don't need to close images in SPM5?
end         
return
