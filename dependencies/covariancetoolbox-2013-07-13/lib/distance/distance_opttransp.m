function a = distance_opttransp(A,B)

A12 = A^(0.5);
a = sqrt(trace(A) + trace(B) - 2 * trace( (A12*B*A12)^0.5));