% Two and threedimensional (centered) Fourier transformation
% where the input can be subject to rigid body motion:
%    matrix vector multiplication
%
% (c) by Alexander Loktyushin, MPI for Biological Cybernetics, 2011 February 19

function y = mvm(A,x,ctransp)

  t = A.t(:).*A.tc(:);                                             % translation
  if strcmp(A.dkind,'none')
    if ctransp                                % translation after after rotation
      y = A.F'*(conj(t).*x(:));
    else
      y = t.*(A.F*x(:));
    end
  elseif strcmp(A.dkind,'all')
      % NO TRANSPOSE
      
      if A.ndims == 3
        [y(:,1) y(:,2) y(:,3) y(:,4)] = trans_gpu(A.trans,x(:),A.imsz,1);
      else
        y = zeros(prod(A.imsz),3);
        y(:,1) = A.F*x(:);
        y(:,2:3) = A.dt.*repmat(A.tc(:).*y(:,1),1,size(A.dt,2));       % Translation
        y(:,1) = t.*y(:,1);                                            % No derivative
      end;
  else
    error('Transpose of rot derivative is not implemented.');
  end