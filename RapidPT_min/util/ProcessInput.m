function [ sub, numPermutations, maxRank, trainNum, maxCycles, iter, write ] = ProcessInput( inputs, N )
%ProcessInput Get all input variables to run permutation testing
% %     %%% INPUTS
% %     A structure filed with following arguments
% %     inputs.datapath    :       path to mat file containing the Data matrix 
% %                                (REQUIRED) Two fields : Data and labeling
% %                                Data - a matrix of size N X V
% %                                labels - a vector of length N (2 groups)
% %                                CONTENTS SHOULD BE NAMED "Data" and "labels"                                
% %                                (N : number of instances , V : Data dimension) 
% %     inputs.sub         :       sub-sampling rate (0 to 1) (DEFAULT = 0.05)
% %     inputs.T           :       number of permutations (DEFAULT = 10^4)
% %     inputs.maxrank     :       rank for estimating the low rank subspace (DEFAULT = N)
% %     inputs.trainNum   :       number of permutations for training (DEFAULT = 100)
% %     inputs.maxCycles   :       number of cycles for training (DEFAULT = 3)
% %     inputs.iter        :       number of iterations for matrix completion (DEFAULT = 30)  
% %     inputs.writing     :       if 0 - outputs only maxnull (SEE BELOW) 
% %                                if 1 - outputs maxnull, U and W (DEFAULT = 0)
% %     inputs.saveDir        :       path to save the outputs (DEFAULT : working folder)  
% %     inputs.timingDir        :       path to save the timings (DEFAULT : working folder)  
    
    fprintf('Rapid Permutation Testing Parameters: \n')
    % Check if each value was assigned in inputs. If they were assigned
    % take the assigned valued, if not assign them a default value.
    if isfield(inputs,'sub') sub = inputs.sub; else sub = 0.05; end
    if isfield(inputs,'T') numPermutations = inputs.T; else numPermutations = 10^4; end
    if isfield(inputs,'maxRank') maxRank = inputs.maxRank; else maxRank = N; end
    if isfield(inputs,'trainNum') trainNum = inputs.trainNum; else trainNum = 100; end
    if isfield(inputs,'maxCycles') maxCycles = inputs.maxCycles; else maxCycles = 3; end
    if isfield(inputs,'iter') iter = inputs.iter; else iter = 30; end    
    if isfield(inputs,'write') write = inputs.write; else write = 0; end    
%     if isfield(inputs,'saveDir') saveDir = inputs.saveDir; else saveDir = pwd; end
%     if isfield(inputs,'timingDir') timingDir = inputs.timingDir; else timingDir = pwd; end

    fprintf('    sub = %d \n', sub)
    fprintf('    numPermutations = %d \n', numPermutations)
    fprintf('    maxRank = %d \n', maxRank)
    fprintf('    trainNum = %d \n', trainNum)
    fprintf('    maxCycles = %d \n', maxCycles)
    fprintf('    iter = %d \n', iter)
    fprintf('    write = %d \n', write)
%     fprintf('    saveDir = %s \n', saveDir)
%     fprintf('    timingDir = %s \n', timingDir)
end

