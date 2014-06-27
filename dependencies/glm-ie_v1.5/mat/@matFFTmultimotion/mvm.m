%    matrix vector multiplication

function y = mvm(A,x,ctransp)
  if A.di
      mask = A.masks(A.di,:,:);
      y = d(A.M{A.di},'all')*(mask(:).*x(:));
  else
    y = zeros(prod(A.imsz),1);
    if ctransp
      for i = 1:A.num_patches
        mask = A.masks(i,:,:);
        y = y + mask(:).*(A.M{i}'*x);
      end;
    else
      for i = 1:A.num_patches
        mask = A.masks(i,:,:);
        y = y + A.M{i}*(mask(:).*x(:));
      end;
    end;
  end;