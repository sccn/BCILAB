function Ceq = vect_cons(X,F,tf,probinfo)

% Calculate h
numintervals = probinfo.numintervals;
h = tf/numintervals;

% Get Reference Indices
k    = probinfo.k;
kbp1 = probinfo.kbp1;
kp1  = probinfo.kp1;

Xk    = X(k,:);
Xkbp1 = X(kbp1,:);
Xkp1  = X(kp1,:);

Fk    = F(k,:);
Fkbp1 = F(kbp1,:);
Fkp1  = F(kp1,:);

% Dynamic Constraints
C1 = Xkbp1 - 1/2.*(Xkp1+Xk) - h/8.*(Fk-Fkp1);
C2 = Xkp1 - Xk - h/6.*(Fkp1 + 4.*Fkbp1 + Fk);

Ceq = [C1(:);C2(:)];