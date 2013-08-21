function answer = get_orthonormal(m,n)
% Produces an m x n set of orthonormal vectors, 
% (thats n vectors, each of length m)
% 
% Inputs should be two scalars, m and n, where n is smaller than 
% or equal to m.
%
% Example: >> get_orthonormal(5,4)
%
% ans =
%     0.1503   -0.0884   -0.0530    0.8839
%    -0.4370   -0.7322   -0.1961   -0.2207
%    -0.3539    0.3098    0.7467   -0.0890
%     0.7890   -0.1023    0.0798   -0.3701
%    -0.1968    0.5913   -0.6283   -0.1585



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CHECK USER INPUT
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if ( (nargin==2) && (m>n) && (isnumeric(m)*isnumeric(n)) )
    
elseif ( nargin==1 && isnumeric(m) && length(m)==1 )
    
    n=m;
    
else
   error('Incorrect Inputs. Please read help text in m-file.')
end




% to get n orthogonal vectors (each of size m), we will first get a larger mxm 
% set of orthogonal vectors, and then just trim the set so it is 
% of size mxn,
%
% to get an mxm set of orthogonal vectors, 
% we can exploit the fact that the eigenvectors
% from distinct eigenspaces (corresponding to different eigenvalues) of a
% symmetric matrix are orthogonal



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Create the orthonormal vectors
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

count=0;
while (count==0)

    % generate an mxm matrix A, then make a symmetric mxm matrix B
    A=rand(m);
    B=A'*A ;

    % now take the eigenvectors of B, 
    % eigenvectors from different eigenspaces (corresponding to different
    % eigenvalues) will be orthogonal

    % there is a chance that there will be repeat eigenvalues, which would give
    % rise to non-orthogonal eigenvectors (though they will still be independent)
    % we will check for this below
    % if this does happen, we will just start the loop over again 
    % and hope that the next randomly created symmetric matrix will not
    % have repeat eigenvalues,
    % (another approach would be to put the non-orthogonal vectors
    % through the gram-schmidt process and get orthogonal vectors that way)

    % since matlab returns unit length eigenvectors, they will also be
    % orthonormal

    [P,D] = eig(B) ;

    % can double check the orthonormality, by taking the difference between 
    % P'P and I, if it is non-zero, then there is an error (caused by repeat
    % eigenvalues), repeat the loop over again

    if ((P'*P - eye(m))>eps) 
        % error, vectors not orthonormal, repeat the random matrix draw again
        count=0
    else
        % we want the first n of these orthonormal columns
        answer=P(:,1:n) ;
        count=1;
    end


end



