function [ subV ] = CheckSamplingRate( N, V, sub)
% min_subV Validate that we are sampling enough entries from the matrix
%   * Condition: 
%`      We want V*subV > r*log(V) 

minSampling = round(N * log(V));
subV = round(sub * V); 
if(minSampling > subV)
    warning('Warning: The input sub-sampling rate is too low! RapidPT requires that V*subV > N*log(V) to be able to accurately recover the Maxnull distribution. See paper reference (Sec 3.1) for more information on the theory behind this requirement.');
    warning('We are automatically changing the sub-sampling rate such that V*sub > 2*N*log(V).');
    minSub = 2*N*log(V)/V;
    subV = round(minSub*V);
    fprintf( 'Old sub-sampling rate, sub = %d \n', sub);
    fprintf( 'New sub-sampling rate, sub = %d, subV = %d \n', minSub, subV);
end
    
    
end

