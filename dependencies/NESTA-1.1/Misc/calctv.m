% Calculates the TV of N1xN2 image
% x is the vector version of the image
function [tv,val] = calctv(N1,N2,x)
X = reshape(x,N1,N2);
tv = sum(sum(sqrt([diff(X,1,2) zeros(N1,1)].^2 + [diff(X,1,1); zeros(1,N2)].^2 )));
val=max(max(sqrt([diff(X,1,2) zeros(N1,1)].^2 + [diff(X,1,1); zeros(1,N2)].^2 )));
end
