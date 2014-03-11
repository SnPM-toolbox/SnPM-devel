function V = snpm_clone_vol(Vtemplate, fname, descrip)
% Clone volume based on a template
% FORMAT V = snpm_clone_vol(Vtemplate, fname, descrip)
%
% Vtemplate ...
%_______________________________________________________________________
% $Id: snpm_clone_vol.m,v 8.1 2009/01/29 15:02:57 nichols Exp $


V        = Vtemplate;
V        = rmfield(V,{'fname','descrip','n','private'});
V.fname  = fname;
V.dim    = V.dim(1:3);
V.dt     = [spm_type('float64'), spm_platform('bigend')];
V.pinfo  = [1 0 0]';
V.descrip= descrip;

return



