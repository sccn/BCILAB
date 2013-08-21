function y = LSQRwrapper(x,transp,A,At)
% y = LSQRwrapper(x,trans,A,At)
% NESTA Version 1.1


% Stephen Becker, 11/12/09

if strcmpi(transp,'transp')
    y = At(x);
elseif strcmpi(transp,'notransp') 
    y = A(x);
end
