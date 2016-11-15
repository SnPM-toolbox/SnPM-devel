function varargout=snpm_combo(Action)
% SnPM Combo: An SnPM toolbox for combined cluster-voxel tests
%
% This is a beta version with ABSOLUTELY NO WORRANTY, still
% being tested.
%
% For more detailed descriptions of the combined tests, please
% refer to our manuscript
% http://www.sph.umich.edu/fni-stat/Docs/Combo.pdf
%
% This toolbox requires SnPM99 in addition to SPM99. To use this
% toolbox, you should run the "Setup" and "Compute" steps for
% SnPM first. Then use this toolbox's "ComboResults" button
% instead of SnPM's "Results" button to assess the results
% with combined tests.
%
%------------------------------------------------------------------
% Based on snpm.m v2.4 Andrew Holmes & Thomas Nichols 
% @(#)snpm_combo.m	3.5 Satoru Hayasaka 04/07/02
%	$Id: snpm_combo.m,v 8.1 2009/01/29 15:02:57 nichols Exp $	

%-Format arguments
%-----------------------------------------------------------------------
if nargin == 0, Action='Init'; end


if strcmp(lower(Action),lower('Init'))
%=======================================================================
global MODALITY
if isempty(spm_figure('FindWin','Menu'))
	snpm
	clc
else
	clc
end
snpm_defaults;
snpm_combo('AsciiWelcome')
spm_help('!Disp','snpm_combo.m','','Graphics')
snpm_combo('CreateMenuWin')




elseif strcmp(lower(Action),lower('AsciiWelcome'))
%=======================================================================
disp( ' SnPM Combo Toolbox                                   ')
disp( ' -----------------------------------------------------')
disp( ' A combined test toolbox for SnPM                     ')
disp( ' Starting up .......................                  ')
fprintf('\n')




elseif strcmp(lower(Action),lower('CreateMenuWin'))
%=======================================================================
%close(findobj(get(0,'Children'),'Tag','SnPM Menu'))

%-Open SnPM menu window
%-----------------------------------------------------------------------
S = get(0,'ScreenSize');
F = figure('Color',[1 1 1]*.8,...
	'Name','SnPM Combo',...
	'NumberTitle','off',...
	'Position',[S(3)/2-200,S(4)/2-140-230,300,200],...
	'Resize','off',...
	'Tag','SnPM Combo',...
	'Pointer','Watch',...
	'MenuBar','none',...
	'Visible','off');

%-Frames and text
%-----------------------------------------------------------------------
axes('Position',[0 0 80/300 280/280],'Visible','Off')
text(0.20,0.535,'SnPM',...
     'FontName','Times','FontSize',24,...
     'Rotation',90,...
	'Color',[1 1 1]*.6);
text(0.4,0.8,'Combo',...
	'FontName','Times','FontSize',64,...
	'Color',[1 1 1]*.6);
text(0.45,0.575,'Combined Test Toolbox',...
	'FontName','Times','FontSize',16,'FontAngle','Italic',...
	'FontWeight','Bold',...
	'Color',[1 1 1]*.6);

uicontrol(F,'Style','Frame','Position',[010 010 280 080],...
	'BackgroundColor',snpm_combo('Colour'));
uicontrol(F,'Style','Frame','Position',[020 020 260 060]);


%-Buttons to launch SnPM functions
%-----------------------------------------------------------------------
uicontrol(F,'String','Combined Results',...
	'Position',[030 035 130 030],...
	'CallBack','snpm_combo_pp',...
	'Interruptible','on',...
	'ForegroundColor','k');
uicontrol(F,'String','?',...
	'Position',[165 035 030 030],...
	'CallBack','spm_help(''snpm_combo.man'')',...
	'Interruptible','on',...
	'ForegroundColor','g');
uicontrol(F,'String','Close',...
	'Position',[200 035 070 030],...
	'CallBack','close(gcf)',...
	'ForegroundColor','r');




set(F,'Pointer','Arrow','Visible','on')



elseif strcmp(lower(Action),lower('Colour'))
%=======================================================================
% snpm_combo('Colour')
%-----------------------------------------------------------------------
% %-Developmental livery
varargout = {[0.7,1.0,0.7], 'Lime Green'};
%-Distribution livery
%varargout = {[0.8 0.8 1.0], 'Diluted Blackcurrent Purple'};



else
%=======================================================================
error('SnPM:UnknownAction', 'Unknown action string')

%=======================================================================
end
