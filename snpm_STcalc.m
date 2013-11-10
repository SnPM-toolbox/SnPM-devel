function STCS = snpm_STcalc(varargin)
% Utility function for spatial statistics
%
% FORMAT STCS = snpm_STcalc('init',nPerm)
% Top-level initialization
%
% FORMAT STCS = snpm_STcalc('update',STCS,ST,XYZ,isPos,perm,pU_ST_ut,DF)
% STCS - Data structure
% ST   - Suprathreshold statistic values  (1 x N)
% XYZ  - Voxel locations of corresponding (3 x N)
% isPos - 1 for positive; 2 for negative
% perm  - Permutation No.
% ST_Ut - Cluster-defining threshold
% DF    - Degrees of freedom, used for t-2-z transformation
%
% FORMAT STCS = snpm_STcalc('double', STCS)
% if bhPerms, STCS needs to be doubled. 
% 
% For internal use only 
%
% FORMAT STCS = snpm_STcalc('initK',STCS,C,isPos,perm)
% Initialization for permutation K
%
%_______________________________________________________________________
% Copyright (C) 2013 The University of Warwick
% Id: snpm_STcalc.m  SnPM13 2013/10/12
% Thomas Nichols, Jun Ding & Hui Zhang
% Optimizations by Darren Gitelman

switch (lower(varargin{1}))
 case 'init'
  nPerm = varargin{2};
  % Notation for documentation below
  %   i - Permutation number
  %   n - isPos; 1=Pos, 2=Neg
  %   k - cluster index
  % For example, LMxT{i,n}{k} contains information about the i-th perm's
  % kth cluster.
  %
  STCS = struct('nPerm',nPerm,...       %-Number of permutations
		             ... 
		'LMxT', [],  ...        %-Local maximum intensity; 
		             ...        % LMxT{i,n}{k} is a vector of a
                             ...        % cluster's local maxima 
		             ...
		'nLMxT',[],  ...        %-Number of local maximum;
                             ...        % nLMxT{i,n}(k) is a cluster's
                             ...        % number of local maxima
		             ...
		'K',    [],  ...        %-Cluster sizes; 
		             ...        % K{i,n}(k) is a clusters size
		             ...        %
		'MxT',  [],  ...        %-Cluster maximum intensity;
		             ...        % MxT{i,n}(k) the max within
                             ...        % a cluster
		             ...        %
		'MxK',  zeros(nPerm,2),...
                             ...        %-Maximum cluster size
		             ...        % MxK(i,n) is global max cluster
                             ...        % size for perm (i,n)
		             ...        %
		'C',    zeros(nPerm,2),...
                             ...        %-Number of clusters
		             ...        % C(i,n) is the number of
                             ...        % clusters in perm (i,n)
                             ...
		'AvgT', [],  ...        %-Mean intensity minus threshold
                             ...        % AvgT{i,n}(k) is a number of a
                             ...        % cluster's average - threshold
                             ...
         	'CMssT', [],  ...       %-Cluster Mass (T_bar-thershold)*cluster size
                             ...        % CMssT(i,n) is a number of a
                             ...        % cluster's mass
                             ...
         	'MxAvgT',zeros(nPerm,2), ...   
                             ...        %-Maximum AvgT in perm(i,n)
                             ...
         	'MxCMssT',[],  ...      %-Maximum cluster mass
                             ...        % MxCMssT(i,n) is max cluster mass of perm (i,n)
		             ...
		'AvgZ', [],  ...        %-Avg intensity for Gaussianized data
                             ...
         	'CMssZ', [], ...        %-Cluster Mass for Gaussianized data
                             ...
         	'MxAvgZ',zeros(nPerm,2), ...   
                             ...        %-Maximum AvgZ in perm(i,n)
                             ...
         	'MxCMssZ',[] ...        %-Maximum CMssZ
                             ...        % MxCMssZ(i,n) is max Gaussanized
                             ...        % cluster mass of perm (i,n)
          );
  return

 case 'initk'
  STCS = varargin{2};
  C = varargin{3};
  isPos = varargin{4};
  perm = varargin{5};
  STCS.LMxT{perm,isPos}   = cell(C,1);      %- local maximum intensity
  STCS.nLMxT{perm,isPos}  = zeros(C,1);     %- number of local maximum intensity
  
  STCS.K{perm,isPos}      = zeros(C,1);     %- Cluster size
  STCS.MxT{perm,isPos}    = zeros(C,1);     %- Cluster maximum intensity;
  STCS.AvgT{perm,isPos}   = zeros(C,1);     %- Mean intensity minus threshold
  STCS.CMssT{perm,isPos}  = zeros(C,1);     %- Cluster Mass (T_bar-thershold)*cluster
                                            %- size
  STCS.AvgZ{perm,isPos}   = zeros(C,1);     %- Avg intensity for
                                            %- Gaussianized data
  STCS.CMssZ{perm, isPos} = zeros(C,1);     %- Cluster mass for
                                            %- Gaussianized data
	

  return
  
 case 'update'
  STCS = varargin{2};
  ST = varargin{3};
  XYZ = varargin{4};
  isPos = varargin{5};
  perm = varargin{6};
  ST_Ut = varargin{7};
  df    = varargin{8};
  
  %- determine whether the user input is a p value or real threshold
  if (ST_Ut < 1),
    ST_Ut = spm_invTcdf(1-ST_Ut,df);
  end 
  
  [N,Z,M,A] = spm_max(ST,XYZ);
  Aindex = spm_clusters(XYZ); %- cluster indexes
    
  STCS = snpm_STcalc('initK',STCS,max(A),isPos,perm);
  
  % Compute the stuff!
  STCS.MxK(perm,isPos) = max(N);
  STCS.C(perm,isPos) = max(A);
  for i=1:max(A)
     d = find(A==i);
     Zd = Z(d);

     STCS.K{perm,isPos}(i)      = N(d(1));
     STCS.MxT{perm,isPos}(i)    = max(Zd);
     STCS.LMxT{perm,isPos}{i}   = Zd';
     STCS.nLMxT{perm,isPos}(i)  = numel(d);
     
     indexd = find(Aindex==i);    %- find the indexes for corresponding cluster
     mST_ST_Ut = mean(ST(indexd))-ST_Ut;
     STCS.AvgT{perm,isPos}(i) = mST_ST_Ut;
     STCS.CMssT{perm,isPos}(i) = N(d(1))*(mST_ST_Ut);
     
     %- for Gaussianized t image

     mt2zST_t2zST_Ut = mean(spm_t2z(ST(indexd),df) -  spm_t2z(ST_Ut,df));
     STCS.AvgZ{perm,isPos}(i) = mt2zST_t2zST_Ut;
     STCS.CMssZ{perm,isPos}(i) = N(d(1))*(mt2zST_t2zST_Ut);

  end
 

  if (~isempty(STCS.AvgT{perm,isPos}))
    STCS.MxAvgT(perm,isPos) = max(STCS.AvgT{perm,isPos});
    STCS.MxAvgZ(perm,isPos) = max(STCS.AvgZ{perm,isPos}); %- for Gaussianized t image
  end
  if (~isempty(STCS.CMssT{perm,isPos}))
    STCS.MxCMssT(perm,isPos) = max(STCS.CMssT{perm,isPos});
    STCS.MxCMssZ(perm,isPos) = max(STCS.CMssZ{perm,isPos}); %- for Gaussianized t image
  end
  return
 
 case 'double'
  STCS = varargin{2};
  
  STCS.MxK=[STCS.MxK; flipud(fliplr(STCS.MxK))];
  STCS.C=[STCS.C; flipud(fliplr(STCS.C))];
  STCS.LMxT=[STCS.LMxT; flipud(fliplr(STCS.LMxT))];
  STCS.nLMxT=[STCS.nLMxT; flipud(fliplr(STCS.nLMxT))];
  STCS.K=[STCS.K; flipud(fliplr(STCS.K))];
  STCS.MxT=[STCS.MxT; flipud(fliplr(STCS.MxT))];


  STCS.AvgT    = [STCS.AvgT;flipud(fliplr(STCS.AvgT))];
  STCS.CMssT   = [STCS.CMssT;flipud(fliplr(STCS.CMssT))];
  STCS.MxAvgT  = [STCS.MxAvgT; flipud(fliplr(STCS.MxAvgT))];
  STCS.MxCMssT = [STCS.MxCMssT; flipud(fliplr(STCS.MxCMssT))];
  
  
  %- for Gaussianized t image
  STCS.AvgZ    = [STCS.AvgZ;flipud(fliplr(STCS.AvgZ))];
  STCS.CMssZ   = [STCS.CMssZ;flipud(fliplr(STCS.CMssZ))];
  STCS.MxAvgZ  = [STCS.MxAvgZ;flipud(fliplr(STCS.MxAvgZ))];
  STCS.MxCMssZ = [STCS.MxCMssZ;flipud(fliplr(STCS.MxCMssZ))];

  STCS.nPermReal = STCS.nPerm*2;  

  return
  
end


  











