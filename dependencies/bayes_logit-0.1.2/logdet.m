function out = logdet(A)
%% out = logdet(A)
%
% computes the log(det(A)) of a positive definite A.
%
% This function is more accurate than log(det(A))

out = 2 * sum(log(diag(chol(A))), 1);