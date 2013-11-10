function V = snpm_clone_vol(Vtemplate, fname, descrip)
% Clone volume based on a template
% FORMAT V = snpm_clone_vol(Vtemplate, fname, descrip)
%
% Vtemplate ...
%_______________________________________________________________________
% Copyright (C) 2013 The University of Warwick
% Id: snpm_clone_vol.m  SnPM13 2013/10/12
% Thomas Nichols

V        = Vtemplate;
V        = rmfield(V,{'fname','descrip','n','private'});
V.fname  = fname;
V.dim    = V.dim(1:3);
V.dt     = [spm_type('float64'), spm_platform('bigend')];
V.pinfo  = [1 0 0]';
V.descrip= descrip;

return



