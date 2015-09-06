function D = midist(u,x)
% compute distance metrix based on mutual information
% compatible with pdist function

[MI,vMI,D] = minfo(u,x);

D = 1-D';
