% =========================================================================
% This function can be called with no input arguments to check whether or
% not the installed version of SnPM is up to date. It does this by looking
% on GitHub to find the latest version release number and comparing it to
% the version number stored in the local version of SnPM.
%
% Authors: Thomas Nichols, Tom Maullin (23/07/2018).
% =========================================================================

function snpm_update()

    % Obtain SnPM version number.
    vsnpm=snpm('ver');
    
    % Look to github for a newer version.
    url='https://github.com/SnPM-toolbox/SnPM-devel/releases/latest';
    [s,stat]=urlread(url);
    
    % Error if we couldn't contact github.
    if ~stat 
      error('Can''t contact GitHub');
    end
    
    % Look for latest SnPM version number.
    [tok,x]=regexp(s,'SnPM [0-9.]*','match','tokens');
    tok=tok{1};
    tok=strrep(tok,' ','');
    
    % Tell the user whether their version is the newest.
    if strcmp(tok,vsnpm)
      msg = 'Your version of SnPM is current.';
    else
      msg = ['A new version (' tok ') is available; please visit ' url ' to update your installation.'];
    end

    if ~nargout, fprintf([blanks(9) msg '\n']); end

end
