function spm_append_96(MAT,X,SizeWarn)
% appends a matrix to a MAT-file (Level 1 format)
% FORMAT spm_append(MAT,X,SizeWarn);
% MAT  - name of .mat file
% X    - matrix 
% SizeWarn - Warning string to be issued if file size approaches 1GB.
%____________________________________________________________________________
%
% spm_append is used to append variables to a matrix in a MAT-
% file without loading the matrix into working memory.  This saves
% memory and time.  spm_append is used by spm_spm as the latter
% cycles over planes saving data for subsequent analysis
%
% The name of the matrix is set to MAT ensuring only one matrix per
% MAT-file.  MAT.mat is created if necessary.
%
% NOTE: this routine will require updating for Level 2 MAT-file format and
% forwards compatibility with new versions of MATLAB
% spm_append will only append to MAT-files in pwd
%
% This is SPM96's spm_append, intended for use with SnPM.
%
%_______________________________________________________________________
% Copyright (C) 2013 The University of Warwick
% Id: spm_append_96.m  SnPM13 2013/10/12
% Thomas Nichols, Andrew Holmes	

global LAST_MAT_SIZE_WARNING % Used to avoid spewing repeated warnings about file size

%----------------------------------------------------------------------------
if ~length(X); return; end
X     = real(X);


% create if necessary
%----------------------------------------------------------------------------
fid   = fopen(fullfile(pwd,[MAT '.mat']),'r');
if fid < 0
        eval([MAT ' = X(:,1);']);
        eval(['save ' MAT '.mat ' MAT ' -V4']);
        spm_append_96(MAT,X(:,[2:size(X,2)]));
	LAST_MAT_SIZE_WARNING = '';
        return
end

% check for number of rows compatibility
%----------------------------------------------------------------------------
[m n] = size(X);
HDR   = fread(fid,5,'int32');
if m ~= HDR(2);  error('SnPM:InvalidSizes', '  Incompatible sizes');  end
fclose(fid);

% append and update nunmber of columns
%----------------------------------------------------------------------------
fid   = fopen([MAT '.mat'],'a');
fwrite(fid,X(:),'double');
Len   = ftell(fid);
if Len>2^(10*3)*0.99  % 99% of 1GB
  if ~strcmp(LAST_MAT_SIZE_WARNING,MAT)
    warning('SnPM:VeryLargeFile', ['Very large file!  (' MAT ')\n' SizeWarn])
    LAST_MAT_SIZE_WARNING = MAT;
  end
end
fclose(fid);
fid   = fopen([MAT '.mat'],'r+');
fseek(fid,8,'bof');
fwrite(fid,(HDR(3) + n),'int32');
fclose(fid);
