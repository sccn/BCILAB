function d = gmmbvl_sqdist(a,b)
% gmmbvl_sqdist - computes pairwise squared Euclidean distances between points

% original version by Roland Bunschoten, 1999
%
% $Name:  $

if size(a,1)==1
  d = repmat(a',1,length(b)) - repmat(b,length(a),1); 
  d = d.^2;
else
  aa = sum(a.*a); bb = sum(b.*b); ab = a'*b; 
  d = abs(repmat(aa',[1 size(bb,2)]) + repmat(bb,[size(aa,2) 1]) - 2*ab);
end
