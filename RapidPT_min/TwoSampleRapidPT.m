%% Efficient permutation testing using Matrix completion
% % the following function computes the max Null statistic distribution 
% % in its current format, the code only uses t-statistics

%%% Original Library:
% % RapidPT: https://github.com/felipegb94/RapidPT

%%% Corresponding paper :
% % Speeding up Permutation Testing in Neuroimaging 
% % C Hinrichs, VK Ithapu, Q Sun, SC Johnson, V Singh
% % NIPS 2013
% % Accelerating Permutation Testing in Neuroimaging through Subspace Tracking
% % F Gutierrez-Barragan VK Ithapu, C Hinrichs, T Nichols, SC Johnson, V Singh
% % Under Preparation


%%% Inputs: 
%       - Data: NxV matrix - N: number of subject, V: number of voxels per subject.
%       - numPermutations: Number of permutation to perform
%       - nGroup1: Number of subjects in group 1
%       - writingVal: Can be 1 or 0. 0 means that we only save the maxnull
%         distribution produced by the algorithm. 1 means that we save the maxnull
%         distribution AND the low-rank basis matrix and coefficient matrix with
%         which you can recover the Permutation testing matrix
%%% Outputs
%       - outputs.MaxT - Maximum t-statistics at each permutation.
%       - outputs.MaxNull - Binned MaxT statistics.
%       - outputs.tThreshold - Calculated tThreshold for the input p-value.
%       - timings.tTraining - Time it took for the first part of RapidPT
%       procedure
%       - timings.tRecovery - Time it took to recover the whole Permutation
%       matrix.
%       - timings.tTotal - Total time it took.

function [ outputs, timings ] = TwoSampleRapidPT(Data, numPermutations, nGroup1, writingVal, RapidPTLibraryPath)
% RapidPT 
%   Modified permutation testing algorithm described in the following paper
%   Speeding up Permutation Testing in Neuroimaging  C Hinrichs, VK Ithapu, Q Sun, SC Johnson, V Singh, NIPS 2013
    
    % N subjects, V voxels (or statistics)
    [N,V] = size(Data);

    % Set keys for input struct
    dataKey = 'data';
    %labelsKey = 'labels';
    nGroup1Key = 'nGroup1';
    %nGroup2Key = 'nGroup2';
    subKey = 'sub';
    TKey = 'T';
    maxRankKey = 'maxRank';
    trainNumKey = 'trainNum';
    maxCyclesKey = 'maxCycles';
    iterKey = 'iter';
    writingKey = 'writing';

%% Fixed values for RapidPT
% The following parameters have been extensively tested and have produced
% good results, independent of the dataset size and number of permutations

    
    if(numPermutations < 10000)      
        subVal = {0.005};
    else
        subVal = {0.0035};
    end
    
    maxCyclesVal = {3}; % Number of cycles for training.
    iterVal = {30}; % Number of iterations for matrix completion.

%% Variable inputs for RapidPT
% These inputs depend on the dataset size.

    nGroup1Val = nGroup1; % Size of group 1
    assert(N > nGroup1, 'nGroup1 cannot be larger than the number of subjects in the data');
    TVal = {numPermutations}; % Number of Permutations.
    maxRankVal = {N}; % Rank for estimating the low rank subspace
    trainNumVal = N; % Default number of training samples


    inputs = struct(dataKey, Data,...
                    nGroup1Key, nGroup1Val,...
                    subKey, subVal,...
                    TKey, TVal,...
                    maxRankKey, maxRankVal,...
                    trainNumKey, trainNumVal,...
                    maxCyclesKey, maxCyclesVal,...
                    iterKey, iterVal,...
                    writingKey, writingVal);
            
    [outputs, timings] = RapidPT(inputs, RapidPTLibraryPath);

    
end



