function C = logeuclidean_geodesic(C1,C2,t)

    C = expm((1-t)*logm(C1) + t*logm(C2));
    