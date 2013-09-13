function varargout=snpm(Action)
% SnPM: nonParametric statistical analysis toolbox for SPM
%_______________________________________________________________________
%  ___       ____  __  __ 
% / __) ___ (  _ \(  \/  )  Statistical nonParametric Mapping toolbox
% \__ \( _ ) )___/ )    (   The Wellcome Department of Cognitive Neurology
% (___/(_)_)(__)  (_/\/\_)  SnPM8
%_______________________________________________________________________
%
% This is SnPM8, the beta release of the nonParametric toolbox for SPM
% SnPM8 is written for SPM8    -    released Friday March 19th, 2010
%
% Please refer to this version as "SnPM8" in papers and
% communications.
%
%                           ----------------
%
% Please report bugs to spm-bugs@fil.ion.ucl.ac.uk. Peculiarities and
% general queries should be raised on the SPM discussion list. We will
% report patches to the list and make updates available if necessary.
%
%                           ----------------
%
% See the "About SnPM" topic (snpm.man) for information on using SnPM
%
% Additional information may be found at the SPMweb site:
%                  http://www.fil.ion.ucl.ac.uk/spm/snpm
% ...where details of the SPM email discussion list can be found:
%                  http://www.fil.ion.ucl.ac.uk/spm/help
%
%                           ----------------
%
%                SnPM, by Andrew Holmes and Thomas Nichols
%                   <snpm-authors@fil.ion.ucl.ac.uk>
%_______________________________________________________________________
% 
%	$Id$	

%-----------------------------functions-called------------------------
% spm
% spm_figure
% spm_help
% snpm_cp
% snpm_pp
% snpm_ui
%-----------------------------functions-called------------------------

%-Parameters
%-----------------------------------------------------------------------
SnPMver    = 'SnPM8';

%-Format arguments
%-----------------------------------------------------------------------
if nargin == 0, Action='Init'; end


if strcmp(lower(Action),lower('Init'))
%=======================================================================
global MODALITY
if isempty(spm_figure('FindWin','Menu'))
	spm('PET')
	clc
else
	clc
end
snpm_defaults;
snpm('AsciiWelcome')
spm_help('!Disp','snpm.m','','Graphics',snpm('Ver'))
snpm_init; % Initialize Matlab Batch system
% snpm('CreateMenuWin')


elseif strcmp(lower(Action),lower('AsciiWelcome'))
%=======================================================================
disp( ' ___       ____  __  __                                                  ')
disp( '/ __) ___ (  _ \(  \/  )  Statistical nonParametric Mapping toolbox      ')
disp( '\__ \( _ ) )___/ )    (   University of Warwick')
disp(['(___/(_)_)(__)  (_/\/\_)  Version: ',snpm('Ver')])
fprintf('\n')


elseif strcmp(lower(Action),lower('Ver'))
%=======================================================================
% snpm('Ver')
varargout = {SnPMver};

elseif strcmp(lower(Action),lower('Colour'))
%=======================================================================
% snpm('Colour')
%-----------------------------------------------------------------------
% %-Developmental livery
% varargout = {[0.7,1.0,0.7], 'Lime Green'};
%-Distribution livery
varargout = {[0.8 0.8 1.0], 'Diluted Blackcurrent Purple'};

else
%=======================================================================
error('Unknown action string')

%=======================================================================
end
