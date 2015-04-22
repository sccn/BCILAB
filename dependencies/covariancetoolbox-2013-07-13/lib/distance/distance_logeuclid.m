function a = distance_logeuclid(A,B)

a = norm((logm(B)-logm(A)),'fro')^2;