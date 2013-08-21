function g=lcvplot(alpha,varargin)
%
% Computes and plots the Likelihood Cross-Validation score (LCV)
% for local fits with different smoothing parameters.
%
% The first argument to lcvplot(), alpha, should be a matrix with one
% or two columns (first column = nearest neighbor component, second
% column = constant component). Each row of this matrix is, in turn,
% passed as the 'alpha' argument to lcv() (and locfit()). The results
% are stored in a matrix, and LCV score ploted against the degrees of
% freedom.

k = size(alpha,1);
z = zeros(k,4);

for i=1:k
  z(i,:) = lcv(varargin{:},'alpha',alpha(i,:));
end;

plot(z(:,3),z(:,4));
xlabel('Fitted DF');
ylabel('LCV');

g = [alpha z];
return;
