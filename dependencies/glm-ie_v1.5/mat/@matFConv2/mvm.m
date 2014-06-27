% 2d convolution matrix: matrix vector multiplication
%
% Michael Hirsch * 20/07/2010
% wrapped by Hannes Nickisch, MPI for Biological Cybernetics, 2010 August 19

function y = mvm(A,x,ctransp)
%keyboard;
  sx = A.sx; sy = A.sy; sF = A.sF; sw = A.sw;                            % sizes
  if ctransp                                                        % transposed
    x = reshape(x,sy);
    y = x .* conj(A.F);
  else                                                          % not transposed
    y = reshape(x,sx) .* A.F;
  end