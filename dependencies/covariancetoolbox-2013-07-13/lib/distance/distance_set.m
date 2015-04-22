function d = distance_set(COV,C,method_dist,arg_dist)

if (nargin<3)||(isempty(method_dist))
    method_dist = 'euclid';
end
if (nargin<4)
    arg_dist = {};
end

Ntrial = size(COV,3);
d = zeros(Ntrial,1);
for i=1:Ntrial
    d(i) = distance(COV(:,:,i),C,method_dist,arg_dist);
end