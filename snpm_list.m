function varargout = spm_list(varargin)
% Support function for snpm_pp.m
% FORMAT spm_list('TxtList',TabDat,c)
% FORMAT spm_list('CSVList',TabDat,ofile)
% FORMAT spm_list('XLSList',TabDat,ofile)
%
% 
%_______________________________________________________________________
% Copyright (C) 2016 The University of Warwick
% snnpm_list.m, based on spm_list.m 6903
% Thomas Nichols


switch lower(varargin{1})

    %======================================================================
    case 'csvlist'            %-Export table to comma-separated values file
    %======================================================================
    % FORMAT snpm_list('CSVList',TabDat,ofile)

        if nargin<2, error('Not enough input arguments.'); end
        TabDat = varargin{2};
        if nargin == 3, ofile = varargin{3};
        else            ofile = [tempname '.csv']; end
        
        fid = fopen(ofile,'wt');
	ncol=size(TabDat.hdr,2);
        fprintf(fid,[repmat('%s,',1,ncol-1) '%d,,,\n'],TabDat.hdr{1,:});
        fprintf(fid,[repmat('%s,',1,ncol) '\n'],TabDat.hdr{2,:});
        fmt = TabDat.fmt;
        [fmt{2,:}] = deal(','); fmt = [fmt{:}];
        fmt(end:end+1) = '\n'; fmt = strrep(fmt,' ',',');
        for i=1:size(TabDat.dat,1)
            fprintf(fid,fmt,TabDat.dat{i,:});
        end
        fclose(fid);
        if nargin == 2, open(ofile); end
    
    %======================================================================
    otherwise                                       %-Unknown action string
    %======================================================================
        error('Unknown action string')
end
