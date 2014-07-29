function z = ldivide(x,y)
% CADA overloaded LDIVIDE function: calls cadabinaryarraymath
z = cadabinaryarraymath(y,x,0,1,'rdivide');