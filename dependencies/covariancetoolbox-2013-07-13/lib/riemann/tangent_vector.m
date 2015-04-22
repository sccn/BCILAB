function Tn = tangent_vector(P,Cov)
Tn = logm(P^-0.5*Cov*P^-0.5);
