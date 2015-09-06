% loss_sqp - squared loss function
%
% Copyright(c) 2009 Ryota Tomioka
% This software is distributed under the MIT license. See license.txt
function [floss, gloss]=loss_sqp(zz, bb)

gloss = zz-bb;
floss = 0.5*sum(gloss.^2);