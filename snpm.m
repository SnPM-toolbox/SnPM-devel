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


elseif strcmp(lower(Action),lower('CreateMenuWin'))
%=======================================================================
close(findobj(get(0,'Children'),'Tag','SnPM Menu'))

%-Open SnPM menu window
%-----------------------------------------------------------------------
S = get(0,'ScreenSize');
F = figure('Color',[1 1 1]*.8,...
	'Name',snpm('Ver'),...
	'NumberTitle','off',...
	'Position',[S(3)/2-200,S(4)/2-140,300,320],...
	'Resize','off',...
	'Tag','SnPM Menu',...
	'Pointer','Watch',...
	'MenuBar','none',...
	'Visible','off');

%-Frames and text
%-----------------------------------------------------------------------
axes('Position',[0 0 80/300 320/320],'Visible','Off')
text(0.5,0.475,'SnPM',...
	'FontName','Times','FontSize',72,...
	'Rotation',90,...
	'VerticalAlignment','middle','HorizontalAlignment','center',...
	'Color',[1 1 1]*.6);

text(0.2,0.96,'Statistical nonParametric Mapping',...
	'FontName','Times','FontSize',16,'FontAngle','Italic',...
	'FontWeight','Bold',...
	'Color',[1 1 1]*.6);

uicontrol(F,'Style','Frame','Position',[095 005 200 290],...
	'BackgroundColor',snpm('Colour'));
uicontrol(F,'Style','Frame','Position',[105 015 180 270]);

%-Buttons to launch SnPM functions
%-----------------------------------------------------------------------
uicontrol(F,'Style','Text',...
	'String','Analysis',...
	'HorizontalAlignment','Center',...
	'Position',[115 250 160 030],...
	'ForegroundColor','w');


str = [...
    'if exist(fullfile(''.'',''SnPM.mat''),''file'')==2 & ',...
    'spm_input({''Current directory contains existing SnPM results files.'',',...
        '[''(pwd = '',pwd,'')''],'' '',',...
        '''Continuing will overwrite existing results!''},1,''bd'',',...
        '''stop|continue'',[1,0],1), tmp=0; else, tmp=1; end, ',...
    'if tmp, snpm_ui;',...
    'end'];

uicontrol(F,'String','Setup',...
	'Position',[115 220-(1-1)*35 130 025],...
	'CallBack',str,...
	'Interruptible','on',...
	'ForegroundColor','k');
uicontrol(F,'String','?',...
	'Position',[250 220-(1-1)*35 025 025],...
	'CallBack','spm_help(''snpm_ui.m'')',...
	'Interruptible','on',...
	'ForegroundColor','g');

uicontrol(F,'String','Compute',...
	'Position',[115 220-(2-1)*35 130 025],...
	'CallBack','snpm_cp',...
	'Interruptible','on',...
	'ForegroundColor','k');
uicontrol(F,'String','?',...
	'Position',[250 220-(2-1)*35 025 025],...
	'CallBack','spm_help(''snpm_cp.m'')',...
	'Interruptible','on',...
	'ForegroundColor','g');

uicontrol(F,'String','Results',...
	'Position',[115 220-(3-1)*35 130 025],...
	'CallBack','snpm_pp',...
	'Interruptible','on',...
	'ForegroundColor','k');
uicontrol(F,'String','?',...
	'Position',[250 220-(3-1)*35 025 025],...
	'CallBack','spm_help(''snpm_pp.m'')',...
	'Interruptible','on',...
	'ForegroundColor','g');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

uicontrol(F,'String','Voxel-cluster Results new!',...
	'Position',[115 220-(4-1)*35 130 025],...
	'CallBack','snpm_combo_pp',...
	'Interruptible','on',...
	'ForegroundColor','k');
uicontrol(F,'String','?',...
	'Position',[250 220-(4-1)*35 025 025],...
	'CallBack','spm_help(''snpm_combo_pp.m'')',...
	'Interruptible','on',...
	'ForegroundColor','g');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




uicontrol(F,'String','About SnPM',...
	'Position',[115 060 160 030],...
	'CallBack','spm_help(''snpm.man'')',...
	'ForegroundColor','g');

uicontrol(F,'String','Close SnPM',...
	'Position',[115 025 160 030],...
	'CallBack','close(gcf)',...
	'ForegroundColor','r');

set(F,'Pointer','Arrow','Visible','on')



elseif strcmp(lower(Action),lower('Colour'))
%=======================================================================
% snpm('Colour')
%-----------------------------------------------------------------------
% %-Developmental livery
% varargout = {[0.7,1.0,0.7], 'Lime Green'};
%-Distribution livery
varargout = {[0.8 0.8 1.0], 'Diluted Blackcurrent Purple'};


elseif strcmp(lower(Action),lower('Quit'))
%=======================================================================
close(findobj(get(0,'Children'),'Tag','SnPM Menu'))

else
%=======================================================================
error('Unknown action string')

%=======================================================================
end
