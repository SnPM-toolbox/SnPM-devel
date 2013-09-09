function varargout = snpm_get_defaults(defstr, varargin)
% Get/set the defaults values associated with an identifier
% FORMAT defval = snpm_get_defaults(defstr)
% Return the defaults value associated with identifier "defstr". 
% Currently, this is a '.' subscript reference into the global  
% "defaults" variable defined in spm_defaults.m.
%
% FORMAT snpm_get_defaults(defstr, defval)
% Sets the defaults value associated with identifier "defstr". The new
% defaults value applies immediately to:
% * new modules in batch jobs
% * modules in batch jobs that have not been saved yet
% This value will not be saved for future sessions of SPM. To make
% persistent changes, edit spm_defaults.m.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% T. Nichols
% Based on spm_get_defaults.m 2696 2009-02-05 20:29:48Z guillaume, Glauche

global SnPMdefs;
if isempty(SnPMdefs)
    snpm_defaults;
end

% % construct subscript reference struct from dot delimited tag string
% tags = textscan(defstr,'%s', 'delimiter','.');
% subs = struct('type','.','subs',tags{1}');
subs = defstr;

if nargin == 1
    varargout{1} = getfield(SnPMdefs, subs);
else
    SnPMdefs = setfield(SnPMdefs, subs, varargin{1});
end
