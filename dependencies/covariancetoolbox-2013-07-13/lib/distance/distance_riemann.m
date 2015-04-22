function a = distance_riemann(A,B)

a = sqrt(sum(log(eig(A,B)).^2));