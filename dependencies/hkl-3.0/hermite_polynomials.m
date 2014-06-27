function H = hermite_polynomials(k,x);
if k==1,
    H = x*0+1;
    return;
end

if k==2,
    H = 2*x;
    return;
end
    
H0 = x*0+1;
H1 = 2*x;

for j=3:k
   H = 2 * x .* H1 - 2 * ( j-2 ) * H0;
   H0 = H1;
   H1 = H;
    
end