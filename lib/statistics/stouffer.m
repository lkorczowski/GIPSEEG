function pcomb = stouffer(p,options)
% Stouffer et al's (1949) unweighted method for combination of 
% independent p-values via z's 
if nargin<2
    options='twosided';
end
    

        
    switch options
        case 'twosided'
        pcomb = (1-erf(sum(sqrt(2) * erfinv(1-2*p))/sqrt(2*length(p))))/2;
        case 'onesided'
           pcomb= normcdf(sum(norminv(p))/sqrt(length(p)));
            
           %pcomb = (1-erf(sum(erfinv(1-p))/sqrt(length(p))));

    end    
