function [tfced] = snpm_tfce(img,H,E,C,dh)
% Threshold-Free Cluster Enhancement as per Smith & Nichols (2009).
% FORMAT [tfced] = snpm_tfce(img,H,E,C,dh);
% img   - the 3D image to be transformed
% H     - height exponent, default = 2
% E     - extent exponent, default = 0.5
% C     - connectivity, default 18 (6 = surface, 18 = edge, 26 = corner)
% df    - size of steps for cluster formation, default = 0.1
%
%   More steps will be more precise but will require more time and memory.
%   The H & E default parameter settings match FSL's randomise/fslmaths.
%   The C default setting matches SPM's default cluster forming. 
%   To match FSL's ramdomise default setting, use 26 instead
%   
%__________________________________________________________________________
% 
% The maximal statistic technique is combined
% with the threshold free cluster enhancement (TFCE) transformation due to
% Smith & Nichols (2009), which obviates the need for arbitrary voxelwise
% cluster-forming thresholds and instead produces continuous correct
% p-values for all voxels. Although some spatial specifity is lost
% relative to purely voxelwise approach, this approach, like cluster
% corrections, is substantially less conservative due to the fact that
% it capitalizes on spatial dependency in the data. 
%
% This function is taken from Mark Thornton's MatlabTFCE toolbox
% see https://github.com/markallenthornton/MatlabTFCE
%
%_______________________________________________________________________
%
% Kyoung whan Choe (kywch@uchicago.edu)

% set cluster thresholds
threshs = 0:dh:max(img(:));
threshs = threshs(2:end);
ndh = length(threshs);

% find positive voxels (greater than first threshold)
nvox = length(img(:));

% find connected components
vals = zeros(nvox,1);
cc = arrayfun(@(x) bwconncomp(bsxfun(@ge,img,x),C), threshs);
for h = 1:ndh
    clustsize = zeros(nvox,1);
    ccc = cc(h);
    voxpercc = cellfun(@numel,ccc.PixelIdxList);
    for c = 1:ccc.NumObjects
        clustsize(ccc.PixelIdxList{c}) = voxpercc(c);
    end
    % calculate transform
    curvals = (clustsize.^E).*(threshs(h)^H);
    vals = vals + curvals;
end
tfced = NaN(size(img));
tfced(:) = vals.*dh;

end

