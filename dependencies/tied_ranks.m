function ranks = tied_ranks(x,tol)
% like tiedrank(), but faster and without the bells & whistles.
%
%                               Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                               2013-07-22

if ~exist('tol','var')
    tol = 10*eps; end
ranks = 1:length(x);
[sorted,order] = sort(x);
tie_ranges = reshape(find(diff([false abs(diff(sorted(:)'))<tol false])),2,[])';
for r=1:size(tie_ranges,1)
    range = tie_ranges(r,1):tie_ranges(r,2);
    ranks(range) = sum(ranks(range))/(range(end)-range(1)+1);
end
ranks(order) = ranks;
ranks = reshape(ranks,size(x));
