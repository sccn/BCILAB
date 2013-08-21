function COREmex(source, options)
%COREMEX           Compile CORE_ functions.
%   COREMEX(SOURCE), where SOURCE is the name of a CORE_ source
%   file, checks if the .mexw32 or .dll file of the same name is older 
%   than any of the following source files: {SOURCE, CORE_library.c, 
%   CORE_mextools.c}.  If any of these files have been modified since 
%   the mex library was compiled (or if no such compiled library exists), 
%   COREMEX calls:
%         mex SOURCE CORE_library.c CORE_mextools.c BLAS
%   where LAPACK includes the static LAPACK library definition files for
%   the currently chosen compiler (LCC or MSVC).  This compiles the MEX 
%   code in SOURCE and links to CORE_library functions and BLAS.
%   
%   COREMEX(SOURCE, OPTIONS) allows specification of command line options
%   to be passed to MEX.
%
%   Example:
%      If the current compiler is set to be LCC and CORE_testfile.c has
%      been modified since CORE_testfile.mexw32 was created (i.e.,
%      'Modified' in Windows),
%          COREmex('CORE_testfile.c', '-v -g')
%      calls
%          mex -v -g CORE_testfile.c ...
%                      CORE_library.c CORE_mextools.c lcc_libmwblas.lib
%
%      If not using the OPTIONS argument, the following syntax is valid:
%          COREMEX CORE_testfile.c

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parse Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (nargin < 2),  options = '';  end
[pat,nam,ext] = fileparts(source);
if (~exist(source, 'file')),
	error(['Source file ' source ' not found.']);
end

%%%%%%%%%%%%%%%%%%%%%% Check if Compile Needed %%%%%%%%%%%%%%%%%%%%%%%
src1 = dir(source);
src2 = dir('CORE_library.c');
src3 = dir('CORE_mextools.c');
try,  mexfile = dir([pat nam '.' mexext]);
catch mexfile = dir([pat nam '.dll']); % if MEXEXT fails, version is <7.1 ??
end

if (~isempty(mexfile))
	srclist = [datenum(src1.date) datenum(src2.date) datenum(src3.date)];
	if (all(datenum(mexfile.date) > srclist)),  return;   end
end

%%%%%%%%%%%%%%%%%%%%%% Choose BLAS definitions %%%%%%%%%%%%%%%%%%%%%%
% We need to figure out which compiler we're using so we can link the
% appropriate static BLAS library definition ...
[stat,rslt] = system('echo %USERPROFILE%');
mexopts = [rslt(1:end-1) '\Application Data\Mathworks\MATLAB\R2007b\mexopts.bat'];
origname = textread(mexopts, 'rem %s', 1, 'headerlines', 1);
switch (upper(origname{1}(1:3)))
	case 'MSV', blas = 'msvc_libmwblas.lib';
	%case 'BCC', blas = 'borland_libmwlapack.lib';
	case 'LCC', blas = 'lcc_libmwblas.lib';
	otherwise,  error('Unable to find the static BLAS definitions for that compiler.');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Compile %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eval(['mex ' options ' ' source ' CORE_library.c CORE_mextools.c ' blas]);
