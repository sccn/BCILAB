% Two and threedimensional (centered) Fourier transformation
% where the input can be subject to rigid body motion:
%    matrix vector multiplication
%
% (c) by Hannes Nickisch & Alexander Loktyushin, 
%                               MPI for Biological Cybernetics, 2011 February 19

function y = mvm(A,x,ctransp)
  
  t = A.t(:).*A.tc(:);                                             % translation
  if strcmp(A.dkind,'none')
    if ctransp                                % translation after after rotation
      y = A.RF'*( conj(t).*x(:) );
    else
      y = t.*(A.RF*x(:)); % Fourier transform followed by rotation & translation
    end
  else
    if strcmp(A.dkind,'rot')                             % rotational derivative
      dy = zeros(prod(A.imsz),A.ndims);                                % dy / dk
      if any(strfind(A.type,'lin')) || any(strfind(A.type,'cub')) ...
                                    || any(strfind(A.type,'fessler'))
        if ctransp
          y = zeros(prod(A.imsz),1);
          warning('No transpose of rotational derivatives implemented.')
        else
          for i=1:A.ndims, dy(:,i) = d(A.RF,i)*x(:); end
          if A.di==0                                           % all derivatives
            y = zeros(prod(A.imsz),size(A.dk,3));
            for dd=1:size(A.dk,3)
              y(:,dd) = t.*reshape(sum(dy.'.*A.dk(:,:,dd),1),[],1);      % trans
              if any(strfind(A.type,'fessler'))
                y(:,dd) = y(:,dd) + A.dtc(:,dd).*A.t(:).*(A.RF*x(:));
              end
            end
          else                                             % a single derivative
            y = t.*sum(dy.'.*A.dk(:,:,A.di),1).';                  % translation
            if any(strfind(A.type,'fessler'))
              y = y + A.dtc(:,A.di).*A.t(:).*(A.RF*x(:));
            end
          end
        end
      else
        y = zeros(prod(A.imsz),1);
        warning('Only linear/cubic/fessler rotational derivs implemented.')
      end
    else                                              % translational derivative
      if ctransp
        y = conj(A.tc(:)) .* x(:);
        if A.di==0                                             % all derivatives
          y = conj(A.dt) .* (y*ones(1,A.ndims));
          for dd=1:A.ndims, y(:,dd) = A.RF'*y(:,dd); end
        else
          y = A.RF'*( conj(A.dt(:,A.di)).*x(:) );
        end
      else
        y = A.RF*x(:);         % proceed as usual since rotation is not affected
        if A.di==0                                             % all derivatives
          y = A.dt.*repmat(A.tc(:).*y,1,size(A.dt,2));       % center correction
        else                                           % a particular derivative
          y = A.dt(:,A.di).*A.tc(:).*y;                      % center correction
        end
      end
    end
  end