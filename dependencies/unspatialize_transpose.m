function X = unspatialize_transpose(X,dims)
% Inverse operation of spatialize_transpose()
X = permute(reshape(X,[size(X,1) dims([1 3:end])]),[2 1 3:length(dims)]);