function C = opttransp_geodesic(C1,C2,t)

    I = eye(size(C1,1));
    C112 = C1^(0.5);
    
    Dx = C112*((C112*C2*C112)^(-0.5))*C112;
    
    D = euclidean_geodesic(I,Dx,t);
    
    C = D*C2*D;
    