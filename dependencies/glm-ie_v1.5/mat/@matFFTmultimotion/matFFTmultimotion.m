% Multi-patch non-rigid motion
%
% A = matFFTmultimotion(sz, type, trans, rot, masks)

% (c) Alexander Loktyushin,     MPI for Biological Cybernetics, 2012

function A = matFFTmultimotion(sz, type, trans, rot, masks)

  A.num_patches = size(trans,1);
  A.imsz = sz; sz = prod(A.imsz)*[1 1];
  A.masks = masks;
  A.M = cell(A.num_patches,1);
  
  for i = 1:A.num_patches
    if any(strfind(lower(type),'mex'))
      A.M{i} = matFastFFTmotion(A.imsz, trans{i}, rot{i}, [], type); 
    elseif any(strfind(lower(type),'gpu'))
      A.M{i} = matGPUFFTmotion(A.imsz, trans{i}, rot{i}, [], type); 
    elseif any(strfind(lower(type),'trans'))
      A.M{i} = matFFTtrans(A.imsz, trans{i}, []); 
    else
        % mat mode
    end
  end;
   
  %% derivatives
  A.di = 0;             % which patch derivative, 0 means all no derivative

  %% construct matrix class, inherit from mat
  A = class(A,mfilename,mat(sz,1));