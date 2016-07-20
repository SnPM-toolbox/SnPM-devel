function varargout=snpm(Action)
% SnPM: nonParametric statistical analysis toolbox for SPM
%_______________________________________________________________________
%  ___       ____  __  __   Statistical nonParametric Mapping toolbox
% / __) ___ (  _ \(  \/  )  Dpt of Statistics & Warwick Manufacturing Group
% \__ \( _ ) )___/ )    (   The University of Warwick
% (___/(_)_)(__)  (_/\/\_)  SnPM13 
%_______________________________________________________________________
%
% This is SnPM13, the SPM batch-compatible release of the nonParametric 
% toolbox for SPM
% SnPM13 is written for SPM8 and SPM12b -released Friday October 11th, 2013
%
% Please refer to this version as "SnPM13" in papers and
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
% Copyright (C) 2013-14 The University of Warwick
% Id: snpm.m  SnPM13.01 2014/01/31
% Thomas Nichols

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
SnPMver    = 'SnPM13.1.03';

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
[Finter,Fgraph,CmdLine] = spm('FnUIsetup','Statistical non-Parametric Mapping (SnPM)');
snpm_init; % Initialize Matlab Batch system
snpm('AsciiWelcome')
spm_help('!Disp','snpm.m','','Graphics',snpm('Ver'))
snpm('CreateMenuWin')

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

elseif strcmpi(Action,'CreateMenuWin')
%=======================================================================
% snpm('CreateMenuWin')
% Initialise batch system
spm_jobman('initcfg');
% Open batch window
spm_jobman('interactive');
disp('SnPM13 tools are available in the SPM batch window under SPM -> Tools -> SnPM')

else
%=======================================================================
error('SnPM:UnknownAction', 'Unknown action string')

%=======================================================================
end
