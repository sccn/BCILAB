% Non-uniformly spaced 2d/3d FFT matrix: matrix vector multiplication
%
% (c) by Hannes Nickisch & Alexander Loktyushin, 
%                                  MPI for Biological Cybernetics, 2011 March 01

function y = mvm(A,x,ctransp)
  
  N = prod(A.imsz);
  if any(strfind(A.type,'fessler'))                            % Fessler's nufft
    if ctransp
      if A.di==0
        Ptx = A.P'*x(:);                         % sparse interpolation matrix P
      else
        Ptx = A.dP{A.di}'*x(:);    % derivative of sparse interpolation matrix P
      end
      Ptx = reshape(full(Ptx), A.fftsz);
      y = prod(A.fftsz)/sqrt(N) * ifftn(Ptx);            % apply oversampled fft
      if A.ndims==2                                               % 2 dimensions
        y = y(1:A.imsz(1),1:A.imsz(2),:);               % eliminate zero padding
      else
        y = y(1:A.imsz(1),1:A.imsz(2),1:A.imsz(3),:);     
      end
      y = conj(A.f) .* y;                         % apply diagonal filter matrix
    else
      Sx = A.f .* reshape(x,A.imsz);              % apply diagonal filter matrix
      FSx = fftn(Sx,A.fftsz)/sqrt(N);                    % apply oversampled fft
      if A.di==0
        y = A.P*FSx(:);                    % apply sparse interpolation matrix P
      else
        y = A.dP{A.di}*FSx(:);      % apply derivative of interpolation matrix P
      end
    end
  elseif any(strfind(A.type,'ferrara'))                        % Ferrara's nufft
    R=2;
    % M_sp is the length of the convolution kernel
    M_sp=A.accuracy; % This gives roughly 6 digits of accuracy
    % The variance, tau, of the Gaussian filter may be different in each dimension

    tau = (pi*M_sp./(A.imsz.*A.imsz*R*(R-.5))); % Suggested value  Greengard [1]
    
    %The length of the oversampled grid
    M_r = R*A.imsz;
    
    % Type-II NUFFT:
    % Precompute E_3, the constant component of the (truncated) Gaussian:
    E_3x(1,1:M_sp) = exp(-((pi*(1:M_sp)/M_r(1)).^2)/tau(1));
    % don't waste (slow) exponential calculations
    E_3x=[fliplr(E_3x(1:(M_sp-1))),1,E_3x];
    E_3y(1,1:M_sp) = exp(-((pi*(1:M_sp)/M_r(2)).^2)/tau(2));
    % don't waste (slow) exponential calculations
    E_3y=[fliplr(E_3y(1:(M_sp-1))),1,E_3y];
    
    % Precompute E_4 (for devonvolution after the FFT)
    kx_vec = (-floor(A.imsz(1)/2)):(ceil(A.imsz(1)/2)-1);
    
    % The Hadamard Inverse of the Fourier Transform of the truncated Gaussian
    E_4x(1:A.imsz(1),1)=sqrt(pi/tau(1))*exp(tau(1)*(kx_vec.^2));
    ky_vec = (-floor(A.imsz(2)/2)):(ceil(A.imsz(2)/2)-1);
    
    % The Hadamard Inverse of the Fourier Transform of the truncated Gaussian
    E_4y(1,1:A.imsz(2))=sqrt(pi/tau(2))*exp(tau(2)*(ky_vec.^2));%  
    
    k = A.k;
    
    %Ensure that the knots fit in a 2*pi interval
    M=size(k,1);
      
    kmin=min(k);
    kmax=max(k);
    bw=kmax-kmin;
    
    scale=(A.imsz-1)./bw;
    shift=-floor(A.imsz/2)-kmin.*scale;
    knots=repmat(scale,[M,1]).*k + repmat(shift,[M,1]);
    knots=mod(2*pi*knots./repmat(A.imsz,[M,1]),2*pi); %shift to [0,2*pi)  

    
    
    if A.ndims==2                                                 % 2 dimensions

      if ctransp                  %  Type I transform
        x = x(:);
        x = x(A.in);
        x = x(A.id);
        
        [f_taur,f_taui] = FGG_Convolution2D(real(x(:)),imag(x(:)),...
                         knots(:),E_3x,E_3y,[M_sp, tau(1), tau(2), M_r(1), M_r(2)]);
        f_tau = reshape(f_taur+sqrt(-1)*f_taui,[M_r(1), M_r(2)]);

        % Perform FFT and deconvolve the result
        F_tau=fftshift(fftn(fftshift(f_tau)))/sqrt(numel(f_tau));
        
        %Chop off excess pixels(the fine spacing in the frequency domain expanded
        %the image in the transform domain)
        F_tau(1:round(.5*(R-1)*A.imsz(1)),:)=[];
        F_tau(A.imsz(1)+1:end,:)=[];
        F_tau(:,1:round(.5*(R-1)*A.imsz(2)))=[];
        F_tau(:,A.imsz(2)+1:end)=[];
        
        %Deconvolve the FFT
    
        % AL: 27.01 Corrects for scaling
        %changed from:
        %y=F_tau.*(E_4x*E_4y)/(M*R*R); 
        % to:
        y=F_tau.*(E_4x*E_4y)/(prod(A.imsz)*R*R); 
    
        %Note that the scaling factor is 1/M, which is 
        %apparently not 1/M_r as shown in eqn (9) in [1]. This makes sense because
        %the scaling factor for the DFT we are attempting to approximate (in eqn 
        %(1) of [1]) is 1/M.      

        % Rescale by factor of 2
        y = y(:)*2;
      
      else                        %  Type II transform
      
        % Zero pad for convolution on finer mesh in frequency domain
        F = reshape(x, A.imsz(1), A.imsz(2));
    
        % AL: 27.01 Corrects for scaling
        %changed from:
        %F = F.*(E_4x*E_4y)/M;
        % to:
        F = F.*(E_4x*E_4y)/(prod(A.imsz));

        % AL: 04.04 padding for odd-sized images
        F = padarray(F,floor(A.imsz/2));
        F = padarray(F,[1 1],'pre');

        f_tau = fftshift(ifftn(ifftshift(F)))*sqrt(numel(F));          %   XX.A: sqrt(numel(F))   for orthonormality
        [f_r,f_i] = FGG_Convolution2D_type2(real(f_tau(:)),imag(f_tau(:)),...
                            knots,E_3x,E_3y,[M_sp, tau(1), tau(2), M_r(1), M_r(2)]);
    
        f = (f_r + 1i*f_i)/2;   % adjust the scaling
        
        yy = zeros(sum(A.in),1);
        yy(A.id) = f(:);
        y = zeros(numel(A.in),1);
        y(A.in) = yy;
      end

    else                                                          % 3 dimensions
        
      E_3z(1,1:M_sp) = exp(-((pi*(1:M_sp)/M_r(3)).^2)/tau(3));
      %don't waste (slow) exponential calculations
      E_3z=[fliplr(E_3z(1:(M_sp-1))),1,E_3z];

      kz_vec = (-floor(A.imsz(3)/2)):(ceil(A.imsz(3)/2)-1);
      %The Hadamard Inverse of the Fourier Transform of the truncated Gaussian
      E_4z(1:A.imsz(3),1)=sqrt(pi/tau(3))*exp(tau(3)*(kz_vec.^2));%
        
        
      if ctransp                  %  Type I transform
          
        x = x(:);
        x = x(A.in);
        x = x(A.id);          
     
        %Perform convolution onto finely-spaced grid
        [f_taur,f_taui]=...
                FGG_Convolution3D(double(real(x(:))),double(imag(x(:))),...
                double(knots(:)),double(E_3x), double(E_3y), double(E_3z), ...
                [double(M_sp), double(tau(1)), double(tau(2)), double(tau(3)),...
                double(M_r(1)), double(M_r(2)), double(M_r(3))]);
        f_tau = reshape(f_taur+sqrt(-1)*f_taui,[M_r(1), M_r(2), M_r(3)]);

        %Perform FFT
        F_tau=fftshift(fftn(ifftshift(f_tau)))/sqrt(numel(f_tau));   %  AL: Correct the scaling with sqrt(numel(f_tau));
        %Chop off excess pixels(the fine spacing in the frequency domain expanded
        %the image in the transform domain)
        F_tau(1:round(.5*(R-1)*A.imsz(1)),:,:)=[];
        F_tau(A.imsz(1)+1:end,:,:)=[];
        F_tau(:,1:round(.5*(R-1)*A.imsz(2)),:)=[];
        F_tau(:,A.imsz(2)+1:end,:)=[];
        F_tau(:,:,1:round(.5*(R-1)*A.imsz(3)))=[];
        F_tau(:,:,A.imsz(3)+1:end)=[];

        %Deconvolve the FFT by Hadamard multiplication with the Gaussian
        S=kron(E_4z,E_4x*E_4y);%Scalars for 3D kronecker product (matrix must be 
                                   %reshaped)
        D=permute(reshape(S,[A.imsz(1),A.imsz(3),A.imsz(2)]),[1,3,2]);
        
        % AL: 27.01.11 Changed from
        %y=F_tau.*D/(M*R*R*R); 
        % to:
        y=F_tau.*D/(prod(A.imsz)*R*R*R); 
        
        %Note that the scaling factor is 1/M, which is 
        %apparently not 1/M_r as shown in eqn (9) in [1]. This makes sense because
        %the scaling factor for the DFT we are attempting to approximate (in eqn 
        %(1) of [1]) is 1/M. 
    
        % AL: Scaling adjusted
        y = abs(y)*3;
        y = y/1.06;
    
        y = y(:);    
      
      else                        %  Type II transform
      
        %Deconvolve

        S=kron(E_4z,E_4x*E_4y);
        %Scalars for 3D kronecker product (matrix must be reshaped)
        D=permute(reshape(S,[A.imsz(1),A.imsz(3),A.imsz(2)]),[1,3,2]);
    
        % AL: 27.01.11 Changed from
        %F=reshape(x, A.imsz).*D/M;
        % to:
        F=reshape(x, A.imsz).*D/(prod(A.imsz));
        
        % AL: 04.04 padding for odd-sized images
        padF = padarray(F,floor(A.imsz/2));
        padF = padarray(padF,[1 1 1],'pre');
        
        f_tau = fftshift(ifftn(ifftshift(padF)))*sqrt(numel(padF));          %   XX.A: sqrt(numel(F))   for orthonormality

        [f_r,f_i]=...
                FGG_Convolution3D_type2(double(real(f_tau(:))),double(imag(f_tau(:))),...
                double(knots),E_3x,E_3y,E_3z,[M_sp, tau(1), tau(2), tau(3), M_r(1), M_r(2), M_r(3)]);

        f = (f_r + 1i*f_i)/3; % scaling adjusted
        f = f * 1.06;
        
        yy = zeros(sum(A.in),1);
        yy(A.id) = f(:);
        y = zeros(numel(A.in),1);
        y(A.in) = yy;
    
      end
    end
  else                                                 % nufft using matResample
    x = A.f*x; % rescaling needed due to zero padding
    if A.di==0
      if ctransp, y = A.F'*(  A.P'      *x); else y =   A.P       *(A.F *x); end
    else
      if ctransp, y = A.F'*(d(A.P,A.di)'*x); else y = d(A.P,A.di) *(A.F *x); end
      y = 2*y*A.imsz(A.di);               % account for pretransformation of A.k
    end
  end