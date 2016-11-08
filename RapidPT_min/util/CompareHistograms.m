function [ kldivergence ] = CompareHistograms( distribution1, distribution2 )
%CompareHistograms Summary of this function goes here
%   Detailed explanation goes here

    binRes = 0.05;
    Tbins = -9:binRes:9;
    
    normalizedDistribution1 = gen_hist(distribution1, Tbins);
    toAdd_MaxT = min(normalizedDistribution1(normalizedDistribution1 > 0))/100;
    normalizedDistribution1 = normalizedDistribution1 + toAdd_MaxT;
    normalizedDistribution1 = normalizedDistribution1./sum(normalizedDistribution1);
    
    normalizedDistribution2 = gen_hist(distribution2, Tbins);
    toAdd_MaxT = min(normalizedDistribution2(normalizedDistribution2 > 0))/100;
    normalizedDistribution2 = normalizedDistribution2 + toAdd_MaxT;
    normalizedDistribution2 = normalizedDistribution2./sum(normalizedDistribution2);

    kldivergence = kldiv2(Tbins', normalizedDistribution1, normalizedDistribution2, 'sym');                 

end