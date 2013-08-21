% transform a complex array into a real array; real and imaginary part are
% stored in an interleaved way along the first dimension
%
% INPUT
%  a complex array of size [  m,n,..]
%
% OUTPUT
%  a real    array of size [2*m,n,..]
%
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2010 August 18

function ab = cx2re(z)

sz = size(z); sz(1) = 2*sz(1);
ab = reshape([real(z(:))'; imag(z(:))'], sz);                         % two rows