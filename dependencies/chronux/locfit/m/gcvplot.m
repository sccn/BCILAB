function g=gcvplot(alpha,varargin)
%
% Computes and plots the Generalized Cross-Validation score (GCV)
% for local fits with different smoothing parameters.
%
% The first argument to gcvplot(), alpha, should be a matrix with one
% or two columns (first column = nearest neighbor component, second
% column = constant component). Each row of this matrix is, in turn,
% passed as the 'alpha' argument to gcv() (and locfit()). The results
% are stored in a matrix, and GCV score ploted against the degrees of
% freedom.

k = size(alpha,1);
z = zeros(k,4);

for i=1:k
  z(i,:) = gcv(varargin{:},'alpha',alpha(i,:));
end;

plot(z(:,3),z(:,4));
xlabel('Fitted DF');
ylabel('GCV');

g = [alpha z];
return;
