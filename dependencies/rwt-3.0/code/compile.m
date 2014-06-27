%    COMPILE compiles the c files and generates mex files.
%

if exist('OCTAVE_VERSION', 'builtin')
  mkoctfile --mex -v -DOCTAVE_MEX_FILE ../mex/mdwt.c   ../lib/src/dwt.c   ../lib/src/init.c ../lib/src/platform.c -I../lib/inc -o omdwt.mex
  mkoctfile --mex -v -DOCTAVE_MEX_FILE ../mex/midwt.c  ../lib/src/idwt.c  ../lib/src/init.c ../lib/src/platform.c -I../lib/inc -o omidwt.mex
  mkoctfile --mex -v -DOCTAVE_MEX_FILE ../mex/mrdwt.c  ../lib/src/rdwt.c  ../lib/src/init.c ../lib/src/platform.c -I../lib/inc -o omrdwt.mex
  mkoctfile --mex -v -DOCTAVE_MEX_FILE ../mex/mirdwt.c ../lib/src/irdwt.c ../lib/src/init.c ../lib/src/platform.c -I../lib/inc -o omirdwt.mex
else
  x = computer();
  if (x(length(x)-1:length(x)) == '64')
    mex -v -largeArrayDims ../mex/mdwt.c   ../lib/src/dwt.c   ../lib/src/init.c ../lib/src/platform.c -I../lib/inc -outdir ../bin
    mex -v -largeArrayDims ../mex/midwt.c  ../lib/src/idwt.c  ../lib/src/init.c ../lib/src/platform.c -I../lib/inc -outdir ../bin
    mex -v -largeArrayDims ../mex/mrdwt.c  ../lib/src/rdwt.c  ../lib/src/init.c ../lib/src/platform.c -I../lib/inc -outdir ../bin
    mex -v -largeArrayDims ../mex/mirdwt.c ../lib/src/irdwt.c ../lib/src/init.c ../lib/src/platform.c -I../lib/inc -outdir ../bin
  else
    mex -v -compatibleArrayDims ../mex/mdwt.c   ../lib/src/dwt.c   ../lib/src/init.c ../lib/src/platform.c -I../lib/inc -outdir ../bin
    mex -v -compatibleArrayDims ../mex/midwt.c  ../lib/src/idwt.c  ../lib/src/init.c ../lib/src/platform.c -I../lib/inc -outdir ../bin
    mex -v -compatibleArrayDims ../mex/mrdwt.c  ../lib/src/rdwt.c  ../lib/src/init.c ../lib/src/platform.c -I../lib/inc -outdir ../bin
    mex -v -compatibleArrayDims ../mex/mirdwt.c ../lib/src/irdwt.c ../lib/src/init.c ../lib/src/platform.c -I../lib/inc -outdir ../bin
  end
end
