function [A,sumV] = quadcof(N,NW,order)
% Helper function to calculate the nonstationary quadratic inverse matrix
% Usage: [A,sumV] = quadcof(N,NW,order)
% N             (number of samples)
% NW: Time bandwidth product
% order: order (number of coefficients, upto 4NW)
%
% Outputs: 
%
% A: quadratic inverse coefficient matrix
% sumV: sum of the quadratic inverse eigenvectors

A = zeros(2*NW-2,2*NW-2,order);
V = quadinv(N,NW);
[P,alpha] = dpss(N,NW,'calc');

for ii = 1:order

  for jj = 1:2*NW
    for kk = 1:2*NW
	A(jj,kk,ii) = sqrt(alpha(jj)*alpha(kk))*...
                      sum(P(:,jj).*P(:,kk).*V(:,ii));
    end;
  end;
end;
sumV=sum(V)/N;

  
