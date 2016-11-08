
%% function for generating histogram

function [hist bins] = gen_hist(X, binsize)

N = numel(X); X = sort(X(:));

if numel(binsize) == 1, bins = min(X):binsize:max(X); else
    bins = binsize; end

hist = zeros(numel(bins), 1);

head = 1;
for i=1:numel(bins)
    while head <= numel(X) && X(head) <= bins(i)
        hist(i) = hist(i) + 1/N;
        head = head + 1;
    end
end

end