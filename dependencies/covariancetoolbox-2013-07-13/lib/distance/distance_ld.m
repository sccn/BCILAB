function a = distance_ld(A,B)

a = sqrt(log(det((A+B)/2))-0.5*log(det(A*B)));