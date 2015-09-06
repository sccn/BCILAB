% loss_hsp - hyperbolic secant loss function
%
% Copyright(c) 2009-2011 Ryota Tomioka
%              2009      Stefan Haufe
% This software is distributed under the MIT license. See license.txt
function [floss, gloss]=loss_hsp(zz, bb)

% floss = -sum(log(sech(bb-zz)./pi));
% gloss = tanh(zz-bb);

zz = zz-bb;
mz = abs(zz);

floss = sum(mz + log(exp(zz-mz)+exp(-zz-mz))-log(2/pi));
gloss = (exp(zz-mz)-exp(-zz-mz))./(exp(zz-mz)+exp(-zz-mz));



