% A = LogaffineMeanCOV(B,epsilon,tol)
%
% Calcul de la moyenne des matrice de covariances selon la distance log affine.
% A : médiane des K matrices NxN
%
% B : Matrice NxNxK
% epsilon : Pas de la descente de gradient
% tol : arret de la descente si le critère < tol


function A = logeuclid_mean(B)

K = size(B,3); % Nombre de matrices
% Initialisation avec la matrice moyenne arythmetique
fc = zeros(size(B,1));
for i=1:K
    fc = fc + logm(B(:,:,i));
end
A = expm( (1/K)*fc);

end

