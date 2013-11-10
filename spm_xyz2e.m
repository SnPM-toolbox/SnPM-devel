function e = spm_xyz2e(XYZ,V)
% Converts co-ordinates to element number
% FORMAT e = spm_xyz2e(xyz,V)
% XYZ    - Coordinates, with xyz triples in columns
% V      - a vector of structures containing image volume information
% e      - element number
%__________________________________________________________________________
%
% Given coordinates, origin and voxel dimensions, spm_xyz2e returns
% the element number corresponding to these coordinates.
%
% This is useful for indexing 3d volumes held MatLab matrices such that x
% varies fastest, then y, and then z.
%
% This function is obsolete in MatLab5, where nD sparse matrix commands
% can be used.
%
% In particular, point lists of coordinates (XYZ) and values (X) can be 
% reconstructed into (3d) images in the above format as follows:
%	>> Ximg = zeros(1,prod(ImDim));
%	>> Ximg(spm_xyz2e(XYZ,[ImDim,VoxDim,Origin])) = X;
%
%_______________________________________________________________________
% Copyright (C) 2013 The University of Warwick
% Id: spm_xyz2e.m  SnPM13 2013/10/12
% Thomas Nichols, Andrew Holmes

%-Condition arguments
%-----------------------------------------------------------------------
if (nargin<2), error('Insufficient arguments'), end

ImDim  = V(1).dim';
M      = V(1).mat(1:3, 1:3);
VoxDim = sqrt(diag(M'*M));
MAT    = V(1).mat;
IMAT   = inv(V(1).mat);
Origin = IMAT(1:3,4);

%-Computation
%-----------------------------------------------------------------------
n       = size(XYZ,2);
tmp     = IMAT * [XYZ; ones(1,n)];
rcp     = round(tmp(1:3,:));
DimMult = cumprod([1,ImDim(1:2)']);
OffSets = [0;1;1]*ones(1,n);
e       = ((rcp-OffSets)'*DimMult')';
