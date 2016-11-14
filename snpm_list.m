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
    case 'txtlist'                                 %-Print ASCII text table
    %======================================================================
    % FORMAT snpm_list('TxtList',TabDat,c)

        if nargin<2, error('Not enough input arguments.'); end
        if nargin<3, c = 1; else c = varargin{3}; end
        TabDat = varargin{2};

        %-Table Title
        %------------------------------------------------------------------
        fprintf('\n\nStatistics: %s\n',TabDat.tit)
        fprintf('%c',repmat('=',1,80)), fprintf('\n')

        %-Table header
        %------------------------------------------------------------------
        fprintf('%s\t',TabDat.hdr{1,c:end-1}), fprintf('%s\n',TabDat.hdr{1,end})
        fprintf('%s\t',TabDat.hdr{2,c:end-1}), fprintf('%s\n',TabDat.hdr{2,end})
        fprintf('%c',repmat('-',1,80)), fprintf('\n')

        %-Table data
        %------------------------------------------------------------------
        for i = 1:size(TabDat.dat,1)
            for j=c:size(TabDat.dat,2)
                fprintf(TabDat.fmt{j},TabDat.dat{i,j});
                fprintf('\t')
            end
            fprintf('\n')
        end
        for i=1:max(1,12-size(TabDat.dat,1)), fprintf('\n'), end
        fprintf('%s\n',TabDat.str)
        fprintf('%c',repmat('-',1,80)), fprintf('\n')

        %-Table footer
        %------------------------------------------------------------------
        for i=1:size(TabDat.ftr,1)
            fprintf([TabDat.ftr{i,1} '\n'],TabDat.ftr{i,2});
        end
        fprintf('%c',repmat('=',1,80)), fprintf('\n\n')

        
    %======================================================================
    case 'xlslist'                                  %-Export table to Excel
    %======================================================================
    % FORMAT snpm_list('XLSList',TabDat,ofile)

        if nargin<2, error('Not enough input arguments.'); end
        TabDat = varargin{2};
        if nargin == 3, ofile = varargin{3};
        else            ofile = [tempname '.xls']; end
        
        d          = [TabDat.hdr(1:2,:);TabDat.dat];
        xyz        = d(3:end,end);
        xyz        = num2cell([xyz{:}]');
        d(:,end+1) = d(:,end);
        d(:,end+1) = d(:,end);
        d(3:end,end-2:end) = xyz;
        xlswrite(ofile, d);
        if nargin == 2, winopen(ofile); end
    
    %======================================================================
    case 'csvlist'            %-Export table to comma-separated values file
    %======================================================================
    % FORMAT snpm_list('CSVList',TabDat,ofile)

        if nargin<2, error('Not enough input arguments.'); end
        TabDat = varargin{2};
        if nargin == 3, ofile = varargin{3};
        else            ofile = [tempname '.csv']; end
        
        fid = fopen(ofile,'wt');
        fprintf(fid,[repmat('%s,',1,7) ',,\n'],TabDat.hdr{1,1:7});
        fprintf(fid,[repmat('%s,',1,7) '%s\n'],TabDat.hdr{2,1:8});
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
