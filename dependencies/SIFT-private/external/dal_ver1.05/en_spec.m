% en_spec - spectrum function for the Elastic-net regularizer
% 
% Copyright(c) 2009 Ryota Tomioka
% This software is distributed under the MIT license. See license.txt

function nm=en_spec(ww,theta)

nm=theta*abs(ww)+0.5*(1-theta)*ww.^2;

