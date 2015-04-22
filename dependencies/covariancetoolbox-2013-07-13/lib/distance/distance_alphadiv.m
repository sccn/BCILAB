function a = distance_alphadiv(A,B,alpha)

a = (4/(1-alpha^2))*log( det(((1-alpha)/2)*A + ((1+alpha)/2)*B)/((det(A)^((1-alpha)/2))*(det(A)^((1+alpha)/2))));