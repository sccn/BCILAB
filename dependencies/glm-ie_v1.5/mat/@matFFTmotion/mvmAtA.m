% Two and threedimensional (centered) Fourier transformation
% where the input can be subject to rigid body motion:
%       matrix vector multiplication with A'*A which is symmetric block Toeplitz
%
% (c) by Hannes Nickisch & Alexander Loktyushin, 
%                                  MPI for Biological Cybernetics, 2011 March 03

function y = mvmAtA(A,x)

if any([strfind(lower(A.type),'near'), strfind(lower(A.type),'lin'), ...
                                             strfind(lower(A.type),'cub')])
  y = mvmAtA(A.RF,x);
else
  y = []; error('Only implemented for types nearest, linear and cubic.');
end