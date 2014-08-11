function L = lagrangian1(x,lambda)

f = objfun(x);
c = confun(x);

L = f + c*lambda;