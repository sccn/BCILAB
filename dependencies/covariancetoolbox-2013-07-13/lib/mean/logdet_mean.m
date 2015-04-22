% from Means of Hermitian positive-definite matrices based on the log-determinant α-divergence function
% Zeineb Chebbi , Maher Moakher.


function A = logdet_mean(B)

K = size(B,3); % Nombre de matrices
A = mean(B,3);
tol = 10^-3;
imp = tol + 1;
while imp>tol
    fc = zeros(size(B,1));

    for i=1:K
        fc = fc + inv(0.5*B(:,:,i) + 0.5*A);
    end

    Anew = inv(fc/K);
    imp = distance_ld(Anew,A);
    A = Anew;
end

function B = mypseudoinv(A)
    [U S] = svd(A);
    k=1;
    K = size(A,1);
    B = zeros(K,K);
   while (k<=K) && (S(k,k)>10^-15)
        B = B+U(:,k)*(1/S(k,k))*U(:,k)';
        k=k+1;
    end
     

