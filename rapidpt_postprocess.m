function [ ] = rapidpt_postprocess(MaxT, SnPMt, XYZ, brain, alpha, params)
%rapidpt_postprocess Summary of this function goes here
%   Detailed explanation goes here

    xSize = params.xdim;
    ySize = params.ydim;
    zSize = params.zdim;

    [ tThresh ] = getTThreshold(MaxT, alpha);
    sigVoxelsIndices = find(SnPMt >= tThresh);
    numSigVoxels = size(sigVoxelsIndices,2);

    coords = floor(repmat(params.origin,1,numSigVoxels) - XYZ(:,sigVoxelsIndices));

    newBrain = zeros(xSize,ySize,zSize);
    for i=1:numSigVoxels
        newBrain(coords(1,i),coords(2,i),coords(3,i)) = 1;
    end

    brain_bin = brain; 
    brain_bin.img = double(newBrain);
    imgStr = 'activeBrain';
    
    save_nii(brain_bin,strcat('outputs/',imgStr,'_',num2str(alpha),'.nii'));
    save(strcat('outputs/coords_',imgStr,'_',num2str(alpha),'.mat'),'coords');


end

