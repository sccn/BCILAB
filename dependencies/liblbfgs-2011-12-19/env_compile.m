SupportFiles = {'liblbfgs-1.10/lib/lbfgs.c'}
IncludeDirectories = {'liblbfgs-1.10/include/'}

% note: if you have a really old platform that doesn't support these instrutions for some reason you can comment out the following line
Defines = {'HAVE_EMMINTRIN_H','HAVE_XMMINTRIN_H','USE_SSE'}

