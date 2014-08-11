function [c,ceq] = confun(x)
% Constraints function

c(1) = 1.5 + x(1) * x(2) - x(1) - x(2); %Inequality constraints
c(2) = -x(1) * x(2)-10; 
% No nonlinear equality constraints
ceq=[];