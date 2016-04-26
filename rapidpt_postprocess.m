function [ ] = rapidpt_postprocess(MaxT, SnPMt, XYZ, brain, alpha, params)
%rapidpt_postprocess Summary of this function goes here
%   Detailed explanation goes here

    xSize = params.xdim;
    ySize = params.ydim;
    zSize = params.zdim;
    numAlphas = size(alpha,2);
    tThresh = zeros(1,numAlphas);
    for i=1:numAlphas
        currAlpha = alpha(1,i);
        tThresh(1,i) = getTThreshold(MaxT, currAlpha);
        sigVoxelsIndices = find(SnPMt >= tThresh(1,i));
        numSigVoxels = size(sigVoxelsIndices,2);

        coords = floor(repmat(params.origin,1,numSigVoxels) - XYZ(:,sigVoxelsIndices));

        newBrain = zeros(xSize,ySize,zSize);
        for j=1:numSigVoxels
            newBrain(coords(1,j),coords(2,j),coords(3,j)) = 1;
        end

        brain_bin = brain; 
        brain_bin.img = double(newBrain);
        imgStr = 'activeBrain';
    
        save_nii(brain_bin,strcat('outputs/',imgStr,'_',num2str(currAlpha),'.nii'));
        save(strcat('outputs/coords_',imgStr,'_',num2str(currAlpha),'.mat'),'coords');
    end
    



end

