function [W Ms crit]=uwedge_c(M,W_est0);                                       
%                                                                              
% Uniformly Weighted Exhaustive Diagonalization using Gauss itEration (U-WEDGE)
% for complex-valued matrices (not necessarily symmetric).                     
%                                                                              
% For non-commertial use only.                                                 
%                                                                              
% Coded by Petr Tichavsky, October 2008                                        
%                                                                              
% Please cite                                                                  
%                                                                              
%  P. Tichavsky, A. Yeredor and J. Nielsen,                                    
%     "A Fast Approximate Joint Diagonalization Algorithm                      
%      Using a Criterion with a Block Diagonal Weight Matrix",                 
%      ICASSP 2008, Las Vegas                                                  
%                                                                              
%  P. Tichavsky and A. Yeredor,                                                
%     "Fast Approximate Joint Diagonalization Incorporating Weight Matrices"   
%     IEEE Transactions of Signal Processing, to appear, 2009.                 
%                                                                              
% Input: M .... the matrices to be diagonalized, stored as [M1 M2 ... ML]      
%        West0 ... initial estimate of the demixing matrix, if available    %                                                                              
% Output: W  .... estimated demixing matrix                                    
%                    such that W * M_k * W' are roughly diagonal               
%         Ms .... diagonalized matrices composed of W*M_k*W'                   
%         crit ... final value of the diagonalization criterion                
%                                                                              
[d Md]=size(M);                                                                
L=floor(Md/d);                        % number of matrices                                                               
Md=L*d;                                                                        
iter=0;                                                                        
eps=1e-4;                                                                      
M0=M(:,1:d);                                                             
if nargin<2                                                                    
   [H E]=eig(M0);                                                              
   W=diag(1./sqrt(abs(diag(E))))*H';                                           
else                                                                           
   W=W_est0;                                                                   
end                                                                            
Ms=M;  Ms2=M;                                                                  
Rs=zeros(d,L);                                                                 
for k=1:L                                                                      
      ini=(k-1)*d;                      
      Ms(:,ini+1:ini+d)=W*M(:,ini+1:ini+d)*W';                                 
      Rs(:,k)=diag(Ms(:,ini+1:ini+d));                                         
      Ms2(:,ini+1:ini+d)=conj(Ms(:,ini+1:ini+d)');                             
end                                                                            
crit=sum(abs(Ms(:)).^2)-sum(abs(Rs(:)).^2);                                    
while crit>eps && iter<100                                                     
  B=real(Rs*Rs');                                                              
  Q=diag(B);                                                                   
  C1=zeros(d,d);                                                               
  for id=1:d                                                                   
    C1(:,id)=sum(Ms(:,id:d:Md).*(ones(d,1)*conj(Rs(id,:))),2)...               
      +conj(sum(Ms2(:,id:d:Md).*(ones(d,1)*conj(Rs(id,:))),2));                
  end                                                                          
  D0=2*(B.*B-Q*Q')+eye(d);                                                     
  A0=eye(d)+(C1'.*B-diag(Q)*C1)./D0;                                           
  W=A0\W;                                                                      
  Raux=W*M0*W';                                                                
  aux=1./sqrt(abs(diag(Raux)));                                                
  W=diag(aux)*W;  % normalize the result                                       
  Ms=W*M;                                                                      
  for k=1:L                                                                    
     ini=(k-1)*d;  
     Ms(:,ini+1:ini+d) = Ms(:,ini+1:ini+d)*W';                                  
     Rs(:,k)=diag(Ms(:,ini+1:ini+d));                                          
     Ms2(:,ini+1:ini+d)=conj(Ms(:,ini+1:ini+d)');                              
  end                                                                          
  crit=(sum(abs(A0(:)))-d)/d^2;                                                
  iter=iter+1;                                                                 
end                                                                            
crit=sum(abs(Ms(:)).^2)-sum(abs(Rs(:)).^2);                                    
end %%%%%%%%%%%%%   of UWEDGE_C                                      
