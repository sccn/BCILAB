function [W_est Ms crit]=wedge(M,H,W_est0,maxnumit)
%
% Weighted Exhaustive Diagonalization using Gauss itEration (WEDGE)
%
% Coded by Petr Tichavsky, March 2008
%
% Please cite
%
%  P. Tichavsky, A. Yeredor and J. Nielsen,
%     "A Fast Approximate Joint Diagonalization Algorithm
%      Using a Criterion with a Block Diagonal Weight Matrix",
%      ICASSP 2008, Las Vegas
%
% Input: M .... the matrices to be diagonalized, stored as [M1 M2 ... ML]
%        H .... diagonal blocks of the weight matrix stored similarly
%                     as M, but there is dd2 blocks, each of the size L x L
%        West0 ... initial estimate of the demixing matrix, if available
%        maxnumit ... maximum number of iterations
%
% Output: W_est .... estimated demixing matrix
%                    such that W_est * M_k * W_est' are roughly diagonal
%         Ms .... diagonalized matrices composed of W_est*M_k*W_est'
%         crit ... monitors the convergence
%
%
[d Md]=size(M);
L=floor(Md/d);
dd2=d*(d-1)/2;
Md=L*d;
iter=0;
eps=1e-7;
improve=10;
if nargin<4
   maxnumit=100;
end   
if nargin<3
   [T E]=eig(M(:,1:d));
   W_est=diag(1./sqrt(diag(E)))*T';
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
critic=10;   
crit=[];  
while critic>eps  && iter<maxnumit
 b11=zeros(dd2,1); b12=b11; b22=b11; c1=b11; c2=c1;
 m=0;
 for id=2:d        
    for id2=1:id-1
        m=m+1; im=(m-1)*L;
        Wm=H(:,im+1:im+L);
        Yim=Ms(id,id2:d:Md);
        Wlam1=Wm*Rs(id,:)';
        Wlam2=Wm*Rs(id2,:)';
        b11(m,1)=Rs(id2,:)*Wlam2;
        b12(m,1)=Rs(id,:)*Wlam2;
        b22(m,1)=Rs(id,:)*Wlam1;
        c1(m,1)=c1(m,1)+Wlam2'*Yim';
        c2(m,1)=c2(m,1)+Wlam1'*Yim';
     end
  end
  det0=b11.*b22-b12.^2;
  d1=(c1.*b22-b12.*c2)./det0;
  d2=(b11.*c2-b12.*c1)./det0;
  m=0;
  A0=eye(d);
  for id=2:d
      A0(id,1:id-1)=d1(m+1:m+id-1,1)';
      A0(1:id-1,id)=d2(m+1:m+id-1,1);
      m=m+id-1;
  end  
  Ainv=inv(A0);  
  W_est=Ainv*W_est;
  Raux=W_est*M(:,1:d)*W_est';
  aux=1./sqrt(diag(Raux));
  W_est=diag(aux)*W_est;  % normalize the result
  for k=1:L
     ini=(k-1)*d;
     Ms(:,ini+1:ini+d) = W_est*M(:,ini+1:ini+d)*W_est';
     Rs(:,k)=diag(Ms(:,ini+1:ini+d));
  end
  critic=(sum(d1.^2)+sum(d2.^2))/dd2;
  crit=[crit critic];
end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% of wedge

