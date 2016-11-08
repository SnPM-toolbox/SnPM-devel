function [ ] = CompareHistogramsVisually(distribution1,distribution2,name1,name2)
%CompareHistogramsVisually Summary of this function goes here
%   Detailed explanation goes here
    
    if ~exist('name1','var')
        name1 = 'distribution1';
    end
    if ~exist('name2','var')
        name2 = 'distribution2';
    end
    
    [counts1,centers] = hist(distribution1, 100);
    [counts2] = hist(distribution2, centers);
    
    plot(centers, counts1, '*r', centers, counts2, '+b');
    title('Maxnull distribution');
    xlabel('Value');
    ylabel('Count');
    
    legend(name1,name2);
        

end