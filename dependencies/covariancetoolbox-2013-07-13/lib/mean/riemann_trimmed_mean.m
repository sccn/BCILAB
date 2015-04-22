% A = riemann_trimmed_mean(B,epsilon,tol)
%
% Calcul du barycentre des matrice de covariances.
% A : baricentre des K matrices NxN
%
% B : Matrice NxNxK
% epsilon : Pas de la descente de gradient
% tol : arret de la descente si le critère < tol


function [A critere niter] = riemann_trimmed_mean(B,args)

N_itermax = 100;
if (nargin<2)||(isempty(args))
    tol = 10^-10;
    trim = 10;
    A = trimmean(B,trim,3);
else
    tol = args{1};
    A = args{2};
    trim = args{3};
end

niter = 0;
fc = 0;

while (niter<N_itermax)
    niter = niter+1;
    % Tangent space mapping
    T = Tangent_space(B,A);
    % sum of the squared distance
    fcn = sum(sum(T.^2));
    % improvement
    conv = abs((fcn-fc)/fc);
    if conv<tol % break if the improvement is below the tolerance
       break; 
    end
    % arithmetic mean in tangent space excluding outliers
    TA =  trimmean(T,trim,2);
    % back to the manifold
    A = UnTangent_space(TA,A);
    fc = fcn;
end

if niter==N_itermax
    disp('Warning : Nombre d''itérations maximum atteint');
end

critere = fc;
