function p = stat_fisher_pcomb(p)
% Fisher's (1925) method for combination of independent p-values [1]
% Code adapted from Bailey and Gribskov (1998) [2]
%
% [1] Fisher RA (1925). Statistical methods for research workers (13th edition). London: Oliver and Boyd.
% [2] Bailey TL, Gribskov M (1998). Combining evidence using p-values: application to sequence homology searches. Bioinformatics, 14 (1) 48-54.
%
% Author: Peter Watson, http://imaging.mrc-cbu.cam.ac.uk/statswiki/FAQ/CombiningPvalues

product=prod(p);
n=length(p);
if n<=0
    error('pfast was passed an empty array of p-values')
elseif n==1
    p = product;
    return
elseif product == 0
    p = 0;
    return
else
    x = -log(product);
    t=product;
    p=product;
    for i = 1:n-1
        t = t * x / i;
        p = p + t;
    end
end