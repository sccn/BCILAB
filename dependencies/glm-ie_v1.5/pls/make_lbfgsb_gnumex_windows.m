mex -v -g -c -f 'C:\gnumex\mexopts_f.bat' solver.f
mex -v -g -c -f 'C:\gnumex\mexopts_c.bat' arrayofmatrices.cpp
mex -v -g -c -f 'C:\gnumex\mexopts_c.bat' lbfgsb.cpp
mex -v -g -c -f 'C:\gnumex\mexopts_c.bat' matlabexception.cpp
mex -v -g -c -f 'C:\gnumex\mexopts_c.bat' matlabmatrix.cpp
mex -v -g -c -f 'C:\gnumex\mexopts_c.bat' matlabprogram.cpp
mex -v -g -c -f 'C:\gnumex\mexopts_c.bat' matlabscalar.cpp
mex -v -g -c -f 'C:\gnumex\mexopts_c.bat' matlabstring.cpp
mex -v -g -c -f 'C:\gnumex\mexopts_c.bat' program.cpp
mex -v -g -f    'C:\gnumex\mexopts_c.bat' -output lbfgsb *.obj
delete *.obj