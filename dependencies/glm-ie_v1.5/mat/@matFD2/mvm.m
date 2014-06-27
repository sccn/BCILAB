% 2d finite difference matrix: matrix vector multiplication
%
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2011 October 20

function y = mvm(A,x,ctransp)
  if ctransp
    switch A.shape
      case 'zero'
        n = prod(A.imsz); 
        x1 = reshape(x(1:n),A.imsz);
        x2 = reshape(x(n+1:end),A.imsz);
        y = vertcat(zeros(1,size(x1,2)),x1(1:end-1,:))-x1 ...
          + horzcat(zeros(size(x1,2),1),x2(:,1:end-1))-x2;
      case 'circ'
        n = prod(A.imsz); 
        x1 = reshape(x(1:n),A.imsz);
        x2 = reshape(x(n+1:end),A.imsz);
        y = x1([end,1:end-1],:)-x1 ...
          + x2(:,[end,1:end-1])-x2;
      case 'valid'
        sz1 = A.imsz-[1,0]; n1 = prod(sz1); 
        x1 = reshape(x(1:n1),sz1);
        x2 = reshape(x(n1+1:end),A.imsz-[0,1]);
        y =   [-x1(1,:); x1(1:end-1,:)-x1(2:end,:); x1(end,:)] ...
            + [-x2(:,1), x2(:,1:end-1)-x2(:,2:end), x2(:,end)];
      case 'same'
        n = prod(A.imsz); 
        x1 = reshape(x(1:n),A.imsz);
        x2 = reshape(x(n+1:end),A.imsz);
        y = [-x1(1,:); x1(1:end-1,:)-x1(2:end,:)] ...
          + [-x2(:,1), x2(:,1:end-1)-x2(:,2:end)];
    end
  else
    x = reshape(x,A.imsz);
    switch A.shape
      case 'zero'
        y = [ vec( vertcat(x(2:end,:),zeros(1,size(x,2)))-x );
              vec( horzcat(x(:,2:end),zeros(size(x,1),1))-x ) ];
      case 'circ'
        y = [ vec( x([2:end,1],:)-x );
              vec( x(:,[2:end,1])-x ) ];
      case 'valid'
        y = [ vec( x(2:end,:)-x(1:end-1,:) );
              vec( x(:,2:end)-x(:,1:end-1) ) ];
      case 'same'
        y = [ vec( [x(2:end,:)-x(1:end-1,:); -x(end,:)] );
              vec( [x(:,2:end)-x(:,1:end-1), -x(:,end)] ) ];
    end
  end
  
function y = vec(x)
  y = x(:);