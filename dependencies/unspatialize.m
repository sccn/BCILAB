function X = unspatialize(X,dims)
% Inverse operation of spatialize()
X = permute(reshape(X',dims([2 1 3:end])),[2 1 3:length(dims)]);