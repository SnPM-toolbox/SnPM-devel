function [ indexMatrix, permutationMatrix1, permutationMatrix2 ] = TwoSampleGetPermutationMatrices(numPermutations, N, nGroup1 )
%GetPermutationMatrices 
%   * indexMatrix: Each column contains the indeces that will be used to
%   make a permutation. Rows 1-nGroup1 are the indeces for group 1 and
%   rows nGroup1-size(labels,1) are the indeces for group 2. This matrix is
%   used when doing serial permutation testing where each permutation is
%   calculated one by one.
%
%   * permutationMatrix1: Matrix composed of 1's and 0's. At each row the
%   columns that contain 1's refer to the corresponding subject in the data
%   matrix at that index and this subject will be in group 1 for that
%   permutations.
%
%   * permutationMatrix2: Same as permutationMatrix2 but for group 2
%   instead of 1.
    rng('default');
    rng('shuffle');
    indexMatrix = zeros(numPermutations, N); 
    permutationMatrix1 = zeros(numPermutations, N);
    permutationMatrix2 = zeros(numPermutations, N);

    for t = 1:1:numPermutations
        currLabels = randperm(N);
        currLabels1 = currLabels(1:(nGroup1));
        currLabels2 = currLabels(nGroup1+1:end);
        permutationMatrix1(t,currLabels1) = 1;
        permutationMatrix2(t,currLabels2) = 1;
        indexMatrix(t,:) = currLabels;
    end

end





