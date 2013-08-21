% 2d convolution matrix: matrix vector multiplication
%
% Michael Hirsch * 20/07/2010
% wrapped by Hannes Nickisch, MPI for Biological Cybernetics, 2010 August 19

function y = mvm(A,x,ctransp)

  sx = A.sx; sy = A.sy; sF = A.sF; sw = A.sw;                            % sizes
  if ctransp                                                        % transposed
    x = reshape(x,sy);
    switch A.shape
      case 'circ'
        x = circshift(x, floor(sw/2));
      case 'same'
        sw2 = floor(sw/2);
        xx = zeros(sF);
        xx(sw2(1)+(1:sy(1)), sw2(2)+(1:sy(2))) = x;
        x = xx;
      case 'valid'
        x = [zeros(sw-1), zeros(sw(1)-1, sy(2));
             zeros(sy(1), sw(2)-1), x];                           % zero padding
    end
    y = ifftn(fftn(x) .* conj(A.F)); 
    y = y(1:sx(1), 1:sx(2));
  else                                                          % not transposed
    y = ifftn(fftn(reshape(x,sx),sF) .* A.F);
    switch A.shape
      case 'circ'
        y = circshift(y, -floor(sw/2));
      case 'same'
        sw2 = floor(sw/2);
        y = y((sw2(1)+1):end, (sw2(2)+1):end);
        y = y(1:sy(1), 1:sy(2));
      case 'valid'
        y = y(sw(1):end, sw(2):end);
    end
  end