% transform a real array into a complex array; real and imaginary part are
% assumed to be stored in an interleaved way along the first dimension
%
% INPUT
%  a real    array of size [2*m,n,..]
%
% OUTPUT
%  a complex array of size [  m,n,..]
%
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2010 August 10

function z = re2cx(ab)

sz = size(ab); sz(1) = sz(1)/2;
z = reshape(ab(1:2:end) + 1i*ab(2:2:end), sz);  