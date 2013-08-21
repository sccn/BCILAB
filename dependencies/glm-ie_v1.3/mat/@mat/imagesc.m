% base matrix: show as an image, relies on full
%
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2010 August 18

function imagesc(A)
  Af = full(A); 
  if ~isreal(Af), Af = [real(Af), -imag(Af); imag(Af), real(Af)]; end
  imagesc(Af)