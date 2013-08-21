function g=aicplot(alpha,varargin)
%
% Computes and plots the -2*AIC
% for local fits with different smoothing parameters.
%
% The first argument to aicplot(), alpha, should be a matrix with one
% or two columns (first column = nearest neighbor component, second
% column = constant component). Each row of this matrix is, in turn,
% passed as the 'alpha' argument to aic() (and locfit()). The results
% are stored in a matrix, and aic score ploted against the degrees of
% freedom.

k = size(alpha,1);
z = zeros(k,4);

for i=1:k
  z(i,:) = aic(varargin{:},'alpha',alpha(i,:));
end;

plot(z(:,3),z(:,4));
xlabel('Fitted DF');
ylabel('AIC');

g = [alpha z];
return;
