%
% FORMAT spm_defaults
%_______________________________________________________________________
%
% This file is intended to be customised for the site.
% Individual users can make copies which can be stored in their own
% matlab subdirectories. If ~/matlab is ahead of the SnPM directory
% in the MATLABPATH, then the users own personal defaults are used.
%
% Care must be taken when modifying this file
%_______________________________________________________________________
% @(#)snpm_defaults.m	3.2 Thomas Nichols 04/08/24
%	$Id: snpm_defaults.m,v 8.1 2009/01/29 15:02:57 nichols Exp $	

global SnPMdefs

% Spatial extent threshold parameters
%------------------------------------------------------------------------
SnPMdefs.STalpha = 0.01; % T values above this sig are save for ST analysis
SnPMdefs.STprop  = 0.10; % 100*(1-STprop)%ile of observed Psuedo T values saved

% When to work volumetrically?
%------------------------------------------------------------------------
SnPMdefs.nMax4DefVol = 16; % Default to volumetric if less than this many scans
