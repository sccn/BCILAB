% base matrix: mtimes, relies on mvm, A.ctransp and A.complex
%
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2010 November 23

function C = mtimes(A,B)

  if isa(A,'mat')
    if isnumeric(B)
      if numel(B)==1 && size(A,2)>1               % scaled matrix: MATRIX*scalar
        C = mat(size(A),1,B,A,'scale');
      else                                                  % MVM: MATRIX*vector

        % here, we specify, how MVMs are evaluated; in general, our computation
        % tree has nodes {hcat,prod,scale,sub,sum,vcat} and leaves {full,zero}
        x = B(:);
        switch A.type
          case  {'hcat','vcat'}     % node: concatenation (horzcat.m, vertcat.m)
              if A.ctransp
                y = 0; j = 0;
                for i=1:length(A.args)
                  si1 = size(A.args{i},1);
                  id = j + (1:si1);
                  y  = y + [A.args{i}']*x(id);
                  j  = j + si1;
                end
              else
                y = []; for i=1:length(A.args), y = [y; A.args{i}*x]; end
              end
          case  'kron'                        % node: Kronecker product (kron.m)
              if A.ctransp
                B = A.args{1}; C = A.args{2};       % kron(B,C)'*vec(x) = C'*x*B
                sb = size(B); sc = size(C);
                x = reshape(x,sc(1),sb(1));
                Ctxt = zeros(sc(2),sb(1));
                for i=1:sb(1), Ctxt(:,i) = [C']*x(:,i); end, Ctxt = Ctxt';
                y = zeros(sb(2),sc(2));
                for i=1:sc(2), y(:,i) = [B']*Ctxt(:,i); end, y = y';
              else
                B = A.args{1}; C = A.args{2};        % kron(B,C)*vec(x) = C*x*B'
                sb = size(B); sc = size(C);
                x = reshape(x,sc(2),sb(2));
                Cxt = zeros(sc(1),sb(2));
                for i=1:sb(2), Cxt(:,i) = C*x(:,i); end, Cxt = Cxt';
                y = zeros(sb(1),sc(1));
                for i=1:sc(1), y(:,i) = B*Cxt(:,i); end, y = y';
              end    
          case  'prod'                                % node: product (mtimes.m)
              if A.ctransp
                y = [A.args{2}']*([A.args{1}']*x);
              else
                y = A.args{1}*(A.args{2}*x);
              end
          case  'rep'                             % node: replication (repmat.m)
              if A.ctransp
                C = A.args{1}; sb = A.args{2}; sc = size(C);  % kron(ones(r),C)'
                x = reshape(x,sc(1),sb(1));
                Ctx = zeros(sc(2),sb(1));
                for i=1:sb(1), Ctx(:,i) = [C']*x(:,i); end
                y = sum(Ctx,2)*ones(1,sb(2));
              else
                C = A.args{1}; sc = size(C); sb = A.args{2};   % kron(ones(r),C)
                x = reshape(x,sc(2),sb(2));
                Cx = zeros(sc(1),sb(2));
                for i=1:sb(2), Cx(:,i) = C*x(:,i); end
                y = sum(Cx,2)*ones(1,sb(1));
              end               
          case  'scale'       % node: scalar multiplication (mtimes.m, uminus.m)
              al = A.args{1}; Amat = A.args{2};
              if A.ctransp, y = al*([Amat']*x); else y = al*(Amat*x); end
          case  'sub'                               % node: indexing (subsref.m)
              Amat = A.args{1}; ii = A.args{2}; jj = A.args{3};
              if A.ctransp
                xf = zeros(size(Amat,1),1);
                xf(ii) = x; yf = [Amat']*xf; y = yf(jj);
              else
                xf = zeros(size(Amat,2),1);
                xf(jj) = x; yf = Amat*xf; y = yf(ii);
              end
          case  'sum'                              % node: sum (minus.m, plus.m)
              if A.ctransp
                y = [A.args{1}']*x + [A.args{2}']*x;
              else
                y =  A.args{1}  *x +  A.args{2}  *x;
              end
          case  'full'                                      % leaf: dense matrix
              if A.ctransp, y = [A.args']*x; else y = A.args*x; end
          case  'zero'              % leaf: matrix-vector-multiplication (mvm.m)
              if strcmp(class(A),'mat')
                y = zeros(size(A,1),1);
              else
                if numel(x)~=size(A,2), error('Wrong size'), end
                if A.ctransp, in = 1; out = 2; else in = 2; out = 1; end
                if A.complex(in)==3                         % input is real only
                  x = real(x);
                elseif A.complex(in)==2              % input represented by IR^2
                  x = re2cx(x);
                else                                       % input is complex IC
                  if A.complex(in)~=1, error('Unknown option for complex'), end
                end
                y = mvm(A,x,A.ctransp);
                if A.complex(out)==3                        % input is real only
                  y = real(y);
                elseif A.complex(out)==2            % output represented by IR^2
                  y = cx2re(y);
                else                                       % input is complex IC
                  if A.complex(out)~=1, error('Unknown option for complex'), end
                end
              end
          otherwise
            error('Unknown type.')
        end
        C = y(:);

      end
    elseif isa(B,'mat')
      sa = size(A); sb = size(B); sz = [size(A,1),size(B,2)];
      if sa(2)~=sb(1), error('Wrong matrix sizes'), end
      C = mat(sz,1,A,B,'prod');
    else
      error('Wrong argument type.')  
    end
  elseif isnumeric(A)                                         % B is of type mat
    if numel(A)==1 && size(B,1)>1                 % scaled matrix: scalar*MATRIX      
      C = mat(size(B),1,A,B,'scale');
    else                                                    % MVM: vector*MATRIX
      C = [mtimes([B'],A')'];
    end
  else
    error('Wrong argument type.')
  end