function snpmui = snpm_bch_ui_Corr
% Design-specific menu items
%
%_______________________________________________________________________
% Based on Volkmar Glauche's MatlabBatch example code (Rev 1716)
% $Id$ Thomas Nichols, Emma Thomas

DesNm = 'MultiSub: Simple Regression; 1 covariate of interest';
DesFile = mfilename;
DesHelp = {'',...
	  'This analysis is appropriate for a correlation of a continuous variable with the imaging data.',...
	  '',...
	  'Number of permutations. The number of possible permutations with this design is n! (read "n factorial") where n is the number of subjects.  In Matlab, use the GAMMA command to check the number of permutations, as in n=10;gamma(n+1)',...
	  '',...
      'At least 4 subjects are needed for a P-value of 0.05 to be obtained; 7 subjects gives 5040 permutations and is perhaps a practical lower limit.',...
	  };

% Set up menu items specific to this plug in

%Enter Covariate Value
CovInt            = cfg_entry;
CovInt.name       = 'Covariate';
CovInt.tag        = 'CovInt';
CovInt.strtype    = 'e';
CovInt.val        = {};
CovInt.num        = [1 Inf];
CovInt.help       = {'This is the variable to correlate with the imaging data.'}; 

% ---------------------------------------------------------------------
% c Vector
% ---------------------------------------------------------------------
c         = cfg_entry;
c.tag     = 'c';
c.name    = 'Vector';
c.help    = {'Vector of covariate values'};
c.strtype = 'e';
c.num     = [Inf 1];
% ---------------------------------------------------------------------
% cname Name
% ---------------------------------------------------------------------
cname         = cfg_entry;
cname.tag     = 'cname';
cname.name    = 'Name';
cname.help    = {'Name of covariate'};
cname.strtype = 's';
cname.num     = [1 Inf];

% ---------------------------------------------------------------------
% mcov Covariate
% ---------------------------------------------------------------------
cov         = cfg_branch;
cov.tag     = 'cov';
cov.name    = 'Covariate';
cov.val     = {c cname };
cov.help    = {'Add a new covariate to your experimental design'};
% ---------------------------------------------------------------------
% generic Covariates
% ---------------------------------------------------------------------
generic         = cfg_repeat;
generic.tag     = 'generic';
generic.name    = 'Covariates of no interest';
generic.help    = {'Covariates of no interest'};
generic.values  = {cov };
generic.num     = [0 Inf];

%% Executable Branch
snpmui = snpm_bch_ui(DesNm,DesFile,DesHelp,{CovInt generic});
