function [] = AddPaths(RapidPTPath)
%AddPaths Add all paths needed

    % Add path mex path and grasta paths
    PATH = RapidPTPath;
    addpath(PATH);
    addpath(strcat(PATH, '/util'));
    addpath(strcat(PATH, '/include/grasta.1.2.0'));
    addpath(strcat(PATH, '/include/grasta.1.2.0/mex'));

end

