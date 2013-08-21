% compilation script of the glm-ie package
%
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2011 February 17

function make(str)
  me = mfilename;                                          % what is my filename
  mydir = which(me); mydir = mydir(1:end-2-numel(me));      % where am I located
  cd(mydir)                                             % go to correct location

if nargin==0                                                           % compile
  cd @matNoise2/private
  mex -O realnoiselet.c
  cd(mydir)

  cd @matResample/private
  mex -O resample_mex.c CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
  cd(mydir)

  cd @matFFTnu/private
  mex -O FGG_Convolution2D.c
  mex -O FGG_Convolution3D.c
  mex -O FGG_Convolution2D_type2.c
  mex -O FGG_Convolution3D_type2.c
  cd(mydir)

  cd @matWav/private
  mex -O fwtn.c wavelet.c
  cd(mydir)
else                                                      % remove ALL mex files
  if strcmp('clean',str)
    !rm -f @matNoise2/private/realnoiselet.mex*
    !rm -f @matResample/private/resample_mex.mex*
    !rm -f @matFFTnu/private/FGG_Convolution*D*.mex*
    !rm -f @matWav/private/fwtn.mex*
  end
end