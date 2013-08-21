
library_path = 'C:\Programmering\Optimization\liblbfgs-1.9';
cflags = [];

% You have (at least) three options. The third is perhaps simplest and
% therefor enabled by default.
%
% (You might need to run "mex -setup" in the Matlab prompt prior to running
%  this script)


% 1) Uncomment when liblfgs compiled using the GNU tools
%
% eval(['mex ' cflags ' -v lbfgs_.c -I' library_path '\include -L ' library_path '\lib\.libs\liblbfgs.a']);



% 2) Uncomment when liblbfgs compiled using MSVC (in Release mode here)
% 
% eval(['mex ' cflags '  -v lbfgs_.c -I' library_path '\include -L ' library_path '\Release\lbfgs.lib']);



% 3) Uncomment to compile the entire package from Matlab
%
%    Uncomment the next line to enable SSE or SSE2 optimizations available
%    in liblbfgs (I don't no whether they'll have any significant effect.
%                 Note also that the program may crash if these
%                 instructions are not available on your processor)                
%    
%cflags ='COMPFLAGS="/DUSE_SSE /D__SSE2__ $COMPFLAGS"';
eval(['mex ' cflags ' -v lbfgs_.c ' library_path '\lib\lbfgs.c -I' library_path '\include']); 