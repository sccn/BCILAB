function C = riemann_geodesic(C1,C2,t)

    C1m12 = C1^(-0.5);
    C112 = C1^(0.5);
    
    C = C112*((C1m12*C2*C1m12)^t)*C112;
    
    