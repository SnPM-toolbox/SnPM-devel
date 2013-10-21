function t = snpm_bch_check_TwoSampT(job)
    t={};
    % If covariate is defined, must have exactly one value per subject
    if numel(job.cov)>0
        numCov = numel(job.cov);
        numScans = numel(job.scans1)+numel(job.scans2);
        for i = 1:numCov
            if numel(job.cov(i).c) ~= numScans
                t = {['Wrong number of values for covariates ''' job.cov(i).cname ''': ' ...
                    num2str(numel(job.cov(i).c)) ' instead of ' num2str(numScans) ' (number of scans).' ] };
            end
        end
    end
end