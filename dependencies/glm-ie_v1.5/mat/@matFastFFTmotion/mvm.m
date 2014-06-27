% Two and threedimensional (centered) Fourier transformation
% where the input can be subject to rigid body motion:
%    matrix vector multiplication
%
% (c) by Alexander Loktyushin, MPI for Biological Cybernetics, 2011 February 19

function y = mvm(A,x,ctransp)

  t = A.t(:).*A.tc(:);                                             % translation
  if strcmp(A.dkind,'none')
    if ctransp                                % translation after after rotation
      y = A.f*(conj(t).*x(:));
      y = resample_mex(A.k(:),y(:),2*A.nuimsz,ctransp,0,2);  
      y = A.F'*y;
    else
      y = A.F*(A.f*x(:));
      if strcmp(A.mode,'mex')
          y = resample_mex(A.k(:),y(:),2*A.nuimsz,ctransp,0,2);  
      else
          y = resample_gpu(A.k(:),y(:),2*A.nuimsz,ctransp,0,2);
      end;
      y = t.*y; % Fourier transform followed by rotation & translation
    end
  elseif strcmp(A.dkind,'all')
      % NO TRANSPOSE
      
      if A.ndims == 3
        y = zeros(prod(A.imsz),7);
        dy = zeros(prod(A.imsz),A.ndims);
        if strcmp(A.mode,'mex')
          [y(:,1) dy(:,1) dy(:,2) dy(:,3)] = resample_mex(A.k(:),A.F*(A.f*x),2*A.nuimsz,ctransp,[0 1 2 3],2);
        else
          [y(:,1) dy(:,1) dy(:,2) dy(:,3)] = resample_gpu(A.k(:),A.F*(A.f*x),2*A.nuimsz,ctransp,[0 1 2 3],2);
        end;         
      
        for i=1:A.ndims
          dy(:,i) = 2*A.nuimsz(i)*dy(:,i);
        end  
          
        for dd=1:size(A.dk,3)
          y(:,dd+1) = t.*reshape(sum(dy.'.*A.dk(:,:,dd),1),[],1);        % Rotation
        end
   
        y(:,5:7) = A.dt.*repmat(A.tc(:).*y(:,1),1,size(A.dt,2));       % Translation
        y(:,1) = t.*y(:,1);                                            % No derivative
      else
        y = zeros(prod(A.imsz),4);
        dy = zeros(prod(A.imsz),A.ndims);
      
        if strcmp(A.mode,'mex')
          [y(:,1) dy(:,1) dy(:,2)] = resample_mex(A.k(:),A.F*(A.f*x),2*A.nuimsz,ctransp,[0 1 2],2);
        else
          [y(:,1) dy(:,1) dy(:,2)] = resample_gpu(A.k(:),A.F*(A.f*x),2*A.nuimsz,ctransp,[0 1 2],2);
        end;        
      
        for i=1:A.ndims
          dy(:,i) = 2*A.nuimsz(i)*dy(:,i);
        end  
          
        for dd=1:size(A.dk,3)
          y(:,dd+1) = t.*reshape(sum(dy.'.*A.dk(:,:,dd),1),[],1);        % Rotation
        end
   
        y(:,3:4) = A.dt.*repmat(A.tc(:).*y(:,1),1,size(A.dt,2));       % Translation
        y(:,1) = t.*y(:,1);                                            % No derivative
      end;
  else
    if strcmp(A.dkind,'rot')                             % rotational derivative
      if ctransp, error('Transpose of rot derivative is not implemented.'), end
      dy = zeros(prod(A.imsz),A.ndims);                                % dy / dk
      if any(strfind(A.type,'lin')) || any(strfind(A.type,'cub')) ...
                                    || any(strfind(A.type,'fessler'))
        if ctransp
          y = zeros(prod(A.imsz),1);
          warning('No transpose of rotational derivatives implemented.')
        else
            
         
          if A.ndims == 3
            if strcmp(A.mode,'mex')
              [dy(:,1) dy(:,2) dy(:,3)] = resample_mex(A.k(:),A.F*(A.f*x),2*A.nuimsz,ctransp,[1 2 3],2);
            else
              [dy(:,1) dy(:,2) dy(:,3)] = resample_gpu(A.k(:),A.F*(A.f*x),2*A.nuimsz,ctransp,[1 2 3],2);
            end;
          else
            if strcmp(A.mode,'mex')
              [dy(:,1) dy(:,2)] = resample_mex(A.k(:),A.F*(A.f*x),2*A.nuimsz,ctransp,[1 2],2);
            else
              [dy(:,1) dy(:,2)] = resample_gpu(A.k(:),A.F*(A.f*x),2*A.nuimsz,ctransp,[1 2],2);
            end;
          end;
          
          for i=1:A.ndims
            dy(:,i) = 2*A.nuimsz(i)*dy(:,i);
          end            

          %for i=1:A.ndims, dy(:,i) = d(A.RF,i)*x(:); end
          if A.di==0                                           % all derivatives
            y = zeros(prod(A.imsz),size(A.dk,3));
            for dd=1:size(A.dk,3)
              y(:,dd) = t.*reshape(sum(dy.'.*A.dk(:,:,dd),1),[],1);      % trans
            end
          else                                             % a single derivative
            y = t.*sum(dy.'.*A.dk(:,:,A.di),1).';                  % translation
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
          for dd=1:A.ndims
              y(:,dd) = A.f*(y(:,dd));
              y(:,dd) = resample_mex(A.k(:),y(:,dd),2*A.nuimsz,ctransp,0,2);
              y(:,dd) = A.F'*y(:,dd);
          end
        else
          y = A.f*( conj(A.dt(:,A.di)).*x(:));
          y = resample_mex(A.k(:),y(:),2*A.nuimsz,ctransp,0,2);
          y = A.F'*y;
        end
      else
        y = A.F*(A.f*x(:));         % proceed as usual since rotation is not affected
        
        if strcmp(A.mode,'mex')
            y = resample_mex(A.k(:),y(:),2*A.nuimsz,ctransp,0,2);
        else
            y = resample_gpu(A.k(:),y(:),2*A.nuimsz,ctransp,0,2);
        end;
        if A.di==0                                             % all derivatives
          y = A.dt.*repmat(A.tc(:).*y,1,size(A.dt,2));       % center correction
        else                                           % a particular derivative
          y = A.dt(:,A.di).*A.tc(:).*y;                      % center correction
        end
      end
    end
  end