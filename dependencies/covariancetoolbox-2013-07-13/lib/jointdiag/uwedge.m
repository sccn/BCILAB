function [W_est Ms crit]=uwedge(M,W_est0)
%
% Uniformly Weighted Exhaustive Diagonalization using Gauss itEration (U-WEDGE)
%
% Coded by Petr Tichavsky, March 2008, updated July 2008
%
% Please cite
%
%  P. Tichavsky, A. Yeredor and J. Nielsen,
%     "A Fast Approximate Joint Diagonalization Algorithm
%      Using a Criterion with a Block Diagonal Weight Matrix",
%      ICASSP 2008, Las Vegas
% or                                                                            
%  P. Tichavsky and A. Yeredor,                                                
%     "Fast Approximate Joint Diagonalization Incorporating Weight Matrices"   
%     IEEE Transactions of Signal Processing, to appear, 2009.                 
%
% Input: M .... the matrices to be diagonalized, stored as [M1 M2 ... ML]
%        West0 ... initial estimate of the demixing matrix, if available
%
% Output: W_est .... estimated demixing matrix
%                    such that W_est * M_k * W_est' are roughly diagonal
%         Ms .... diagonalized matrices composed of W_est*M_k*W_est'
%         crit ... stores values of the diagonalization criterion at each
%                  iteration
%
[d Md]=size(M);
L=floor(Md/d);
Md=L*d;
iter=0;
eps=1e-7;
improve=10;
if nargin<2
   [H E]=eig(M(:,1:d));
   W_est=diag(1./sqrt(abs(diag(E))))*H';
else
   W_est=W_est0;
end 
Ms=M; 
Rs=zeros(d,L);
for k=1:L
      ini=(k-1)*d;
      M(:,ini+1:ini+d)=0.5*(M(:,ini+1:ini+d)+M(:,ini+1:ini+d)');
      Ms(:,ini+1:ini+d)=W_est*M(:,ini+1:ini+d)*W_est';
      Rs(:,k)=diag(Ms(:,ini+1:ini+d));
end
crit=sum(Ms(:).^2)-sum(Rs(:).^2);
while improve>eps && iter<100
  B=Rs*Rs';
  for id=1:d
      C1(:,id)=sum(Ms(:,id:d:Md).*Rs,2);
  end
  D0=B.*B'-diag(B)*diag(B)';
  A0=eye(d)+(C1.*B-diag(diag(B))*C1')./(D0+eye(d));
  W_est=A0\W_est;
  Raux=W_est*M(:,1:d)*W_est';
  aux=1./sqrt(abs(diag(Raux)));
  W_est=diag(aux)*W_est;  % normalize the result
  for k=1:L
     ini=(k-1)*d;
     Ms(:,ini+1:ini+d) = W_est*M(:,ini+1:ini+d)*W_est';
     Rs(:,k)=diag(Ms(:,ini+1:ini+d));
  end
  critic=sum(Ms(:).^2)-sum(Rs(:).^2);
  improve=abs(critic-crit(end));
  crit=[crit critic];
  iter=iter+1;
end 
 


