% A = riemann_mean(B,epsilon,tol)
%
% Calcul du barycentre des matrice de covariances.
% A : baricentre des K matrices NxN
%
% B : Matrice NxNxK
% epsilon : Pas de la descente de gradient
% tol : arret de la descente si le crit�re < tol


function [A critere niter] = riemann_diag_mean(B,args)

N_itermax = 10;
if (nargin<2)||(isempty(args))
    tol = 10^-10;
    A = mean(B,3);
else
    tol = args{1};
    A = args{2};
end

niter = 0;
fc = 0;
Ntrial = size(B,3);
Tn = zeros(size(B,1),Ntrial);
B2 = zeros(size(B));
while (niter<N_itermax)
    niter = niter+1;
    % Tangent space mapping
    P = A^-0.5;
    for i=1:Ntrial
        B2(:,:,i) = P*B(:,:,i)*P;
    end
    W = uwedge([eye(size(A)) B2(:,:)]);
    for i=1:Ntrial 
        Tn(:,i) = (log(diag(W*B2(:,:,i)*W')));
    end
    % sum of the squared distance
    fcn = sum(sum(Tn.^2));
    % improvement
    conv = abs((fcn-fc)/fc);
    if conv<tol % break if the improvement is below the tolerance
       break; 
    end
    % arithmetic mean in tangent space
    V = inv(W);
    %diag(mean(Tn,2))
    %TA = V*diag(mean(Tn,2))*V';
    % back to the manifold
    P = A^0.5;
    A = P*V*diag(exp(mean(Tn,2)))*V'*P;
    %A = P*expm(TA)*P;
    fc = fcn
end

if niter==N_itermax
    disp('Warning : Nombre d''it�rations maximum atteint');
end

critere = fc;
