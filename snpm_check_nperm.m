function snpm_check_perm(nPerm,TotPerm)
%
% Issue consistent warning messages about number of permutations
%
%
%_______________________________________________________________________
% Thomas Nichols
% $Id$

% Use very arbitrary heuristic... if fewer than 100 permutations used
% *and* less than 90% of possible permutations considered, then issue a
% specific warning about variability of Monte Carlo variation.

str='';
if nPerm<100 & nPerm>=0.9*TotPerm
  str = sprintf(['Very few (%d) permutations used, nonparametric P-values are very coarse (but exact).\n'...
		 'Smallest possible P-value is %0.4f.'],nPerm,1/nPerm);
elseif nPerm<100 & nPerm<0.9*TotPerm
  str = sprintf(['Very few (%d) permutations used, nonparametric P-values are very coarse will vary \n'...
		 'substantially over repeated re-analyses.'...
		 'Smallest possible P-value is %0.4f.'],nPerm,1/nPerm);
elseif nPerm>10000
  str = sprintf(['Many (%d) permutations used, analysis may take a very long time\n'...
		 'to complete.  Consiser running fewer permutations.'],nPerm,1/nPerm);
elseif nPerm>25000
  str = sprintf(['%cInsane nummber of permutations requested (%d)!\n'...
		 'Are you sure!?'],7,nPerm);
end
if ~isempty(str)
  warning(str)
end
