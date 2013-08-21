function ok = hlp_trycompile(varargin)
% Try to auto-compile a set of binary files in a folder, and return the status.
% OK = hlp_trycompile(Options...)
%
% This function tries to ensure that a given set of functions or classes (specified by their
% MATLAB identifier), whose source files are assumed to be located in a given directory, are
% properly compiled.
%
% The Style parameter determines how the function proceeds: Either compilation is always done
% ('force'), or only if necessary ('eager', e.g. if invalid file or changed source code). The check
% for re-compilation may be done on every call ('eager'), or once per MATLAB session ('lazy').
%
% The most common use case is specifying a given directory (and omitting the identifiers). In this 
% case, all source files that have a mexFunction declaration are compiled, and all other source 
% files are also supplied to the compiler as additional files. In case of compilation errors,
% hlp_trycompile also tries to omit all additional files during compilation. Both the list of 
% additional files (or their file name patterns) or the considered file types can be specified. The 
% list of identifiers to consider can also be specified.
%
% Another possible use case is to omit both the identifiers and the directory. In this case, 
% hlp_trycompile assumes that a mex file with the same identifier (and path) as the calling function
% shall be compiled. This is would be used in .m files which directly implement some fallback code
% in case that the compilation fails (or which are just stubs to trigger the on-demand compilation).
%
%
% Since there can be many different versions of a mex binary under Linux (even with the same name),
% .mex files are by default moved into a sub-directory (named according to the hostname) after
% compilation. This does not apply to .class files, which are not platform-specific.
%
% The function supports nearly all features of the underlying MEX compiler, and can thus be used
% to compile a large variety of mex packages found in the wild (in some cases with custom defines,
% libraries, or include directories).
%
% If you are debugging code with it, it is best to set Verbose to true, so that you get compilation 
% output.
%
% Additional features of this function include:
%  * Figures out whether the -largeArrayDims switch should be used.
%  * Repeated calls of this function are very fast if used in 'lazy' mode (so that it can be used in 
%    an inner loop).
%  * Automatically rebuilds if the mex source has changed (does not apply to misc dependency files),
%    if used in 'eager' mode.
%  * By default uses the Mathworks versions of BLAS and LAPACK (if these libraries are pulled in).
%  * Supports both '/' and '\' in directory names.
%  * Behaves reasonably in deployed mode (i.e. gives warnings if files are missing).
%  * Also compiles .java files where appropriate.
%  * Supports test code.
%
% If this function produces errors for some mex package, the most common causes are:
%  * If the platform has never been used to compile code, an appropriate compiler may have to be 
%    selected using "mex -setup". If no supported compiler is installed (e.g. on Win64), it must 
%    first be doenloaded and installed (there is a free standard compiler for every platform).
%  * Some unused source files are in the directory which produce errors when they are automatically 
%    pulled in.
%    --> turn on verbose output and identify & remove these (or check the supplied make file for 
%        what files are actually needed)
%  * The functions require a custom define switch to work.
%    --> Check the make file, and add the switch(es) using the 'Defines' parameter.
%  * The functions use non-standard C code (e.g. // comments).
%    --> Tentatively rename the offending .c files into .cpp.
%  * The functions require a specific library to work.
%    --> Check the make file, and add the libraries using the 'Libaries' parameter.
%  * The functions require specific include directories to work.
%    --> Check the make file, and add the directories using the 'IncludeDirectories' parameter.
%  * The functions require additional files that are in a different directory.
%    --> Check the make file, and add these files using the 'SupportFiles' parameter. Wildcards are 
%        allowed (in particular the special '*' string, which translates into all source files in 
%        the Directory).
%  * The package assumes that mex is used with the -output option to use a custom identifier name
%    --> This type of make acrobatic is not supported by hlp_trycompile; instead, rename the source
%        file which has the mexFunction definition such that it matches the target identifier.
%  * The functions require specific library directories to work.
%    --> Check the make file, and add the directories using the 'LibraryDirectories' parameter.
%
%
% In:
%   Style : execution style, can be one of the following (default: 'lazy')
%           'force' : force compilation (regardless of whether the binaries are already there)
%           'eager' : compile only if necessary, check every time that this function is called
%           'lazy'  : compile only if necessary, and don't check again during this MATLAB session
%
%   --- target files ---
%
%   Directory : directory in which the source files are located
%               (default: directory of the calling function)
%
%   Identifiers : identifier of the target function/class, or cell array of identifiers that should
%                 be compiled (default: Calling function, if no directory given, or names of all
%                 compilable source files in the directory, if a directory is given.)
%
%   FileTypes : file type patterns to consider as sources files for the Identifiers
%               (default: {'*.f','*.c','*.cpp','*.java'})
%
%
%   --- testing conditions ---
%
%   TestCode : MATLAB code (string) which evaluates to true if the compiled code is behaving
%              correctly (and false otherwise), or alternatively a function handle which does the
%              same
%
%
%   --- additional compiler inputs ---
%
%   SupportFiles : cell array of additional / supporting source filenames to include in the compilation of all
%                  Identifiers (default: '*')
%                  Note: Any file listed here will not be considered part of the Identifiers, when
%                        all contents of a directory are to be compiled.
%                  Note: If there are support source files in sub-directories, include the full path
%                        to them.
%                  Note: If this is '*', all source files that are not mex files in the given 
%                        directory are used as support files.
%
%   Libraries : names of libraries to include in the compilation
%               (default: {})
%
%   IncludeDirectories : additional directories in which to search for included source files.
%                        (default: {})
%
%   LibraryDirectories : additional directories in which to search for referenced library files.
%                        (default: {})
%
%   Defines : list of defined symbols (either 'name' or 'name=value' strings)
%             (default: {})
%
%   Renaming : cell array of {sourcefile,identifier,sourcefile,identifier, ...} indicating that 
%              the MEX functions generated from the respective source files should be renamed to 
%              the given identifiers. Corresponds to MEX's -output option; does not apply to Java files.
%              (default: {})
%
%   Arguments : miscellaneous compiler arguments (default: {})
%               For possible arguments, type "help mex" in the command line
%
%   DebugBuild : whether to build binaries in debug mode (default: false)
%
%
%   --- user messages ---
%
%   ErrorMessage : the detail error message to display which describes what type of functionality
%                  will not be available (if any).
%                  Note: If you have a MATLAB fallback, mention this in the error message.
%
%   PreparationMessage : the message that will be displayed before compilation begins.
%                        (default: {})
%
%   Verbose : whether to display verbose compiler outputs (default: false)
%
%
%   --- misc options ---
%
%   MathworksLibs : whether to use the Mathworks versions of std. libraries instead of OS-supplied 
%                   ones, if present (applies to blas and lapack) (default: true)
%
%   DebugCompile : debug the compilation process; halts in the directory prior to invoking mex
%                  (default: false)
%
%
% Examples:
%   % try to compile all mex / Java files in a given directory:
%   hlp_trycompile('Directory','/Extern/MySources');
%
%   % as before, but restrict the set of identifiers to compile to a given set
%   hlp_trycompile('Directory','/Extern/MySources','Identifiers',{'svmtrain','svmpredict'});
%
%   % try to compile mex / Java files in a given directory, and include 2 libraries in the compilation
%   hlp_trycompile('Directory','/Extern/MySources','Libraries',{'blas','lapack'});
%
%   % like before, but this time include additional source files from two other directories
%   % (the single '*' means: include non-mex sources in the specified directory)
%   hlp_trycompile('Directory','/Extern/MySources','SupportFiles',{'*','../blas/*.c','../*.cpp'});
%
%   % like before, but this time add an include directory, a library directory, and some library
%   hlp_trycompile('Directory','/Extern/MySources', 'IncludeDirectories','/boost/include','LibraryDirectories','/boost/lib','Libraries','boost_date_time-1.44');
%
%   % like before, this time specifying some custom #define's
%   hlp_trycompile('Directory','/Extern/MySources','Defines',{'DEBUG','MAX=((a)>(b)?(a):(b))'});
%
%
% Use cases:
%   1) In addition to a source file mysvd.c (compiling into mysvd.mex*), a stub .m file of the same
%   name can be placed in the same directory, which contains code to compile the binary when needed.
%
%   function [U,S,V] = mysvd(X)
%   if hlp_trycompile
%       [U,S,V] = mysvd(X);
%   else
%       % either display an error message or implement some fallback code.
%   end
%
%   2) In a MATLAB function which makes use of a few mex files, ensure compilation of these files.
%   function myprocessing(X,y)
%   if ~hlp_trycompile('Identifiers',{'svmtrain.c','svmpredict.c'})
%       error('Your binary files could not be compiled.');
%   else
%       m = svmtrain(X,y);
%       l = svmpredict(X,m);
%       ...
%   end
%
%   3) In a startup script.
%   hlp_trycompile('Directory','/Extern/MySources');
%
% See also:
%   mex
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2011-03-09

% Copyright (C) Christian Kothe, SCCN, 2011, christian@sccn.ucsd.edu
%
% This program is free software; you can redistribute it and/or modify it under the terms of the GNU
% General Public License as published by the Free Software Foundation; either version 2 of the
% License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
% even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License along with this program; if not,
% write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
% USA


persistent results;           % a map of result-tag to OK value
persistent compiler_selected; % whether a compiler has been selected (true/false, or [] if uncertain)

% read options
o = hlp_varargin2struct(varargin, ...
    ... % overall behavior
    {'style','Style'}, 'lazy', ...
    ... % target files
    {'dir','Directory'}, [], ...
    {'idents','Identifiers'}, [], ...
    {'types','FileTypes'}, {'*.c','*.C','*.cpp','*.CPP','*.f','*.F','*.java','*.Java'}, ...
    ... % test condition
    {'test','TestCode'},'', ...
    ... % additional compiler inputs
    {'support','SupportFiles'}, {'*'}, ...
    {'libs','Libraries'}, {}, ...
    {'includedirs','IncludeDirectories'}, {}, ...
    {'libdirs','LibraryDirectories'}, {}, ...
    {'defines','Defines'}, {}, ...
    {'args','Arguments'}, '', ...
    {'renaming','Renaming'}, {}, ...
    {'debug','DebugBuild'}, false, ...
    ... % messages
    {'errmsg','ErrorMessage'}, {'Some BCILAB functionality will likely not be available.'}, ...
    {'prepmsg','PreparationMessage'}, {}, ...
    {'verbose','Verbose'}, false, ...
    ... % misc
    {'mwlibs','MathworksLibs'}, true, ...
    {'debugcompile','DebugCompile'}, false ...
    );

% support for parameterless calls
if isempty(o.dir)
    % if no dir given, use the calling function's directory
    [name,file] = hlp_getcaller();
    o.dir = fileparts(file);
    if isempty(file)
        error('If hlp_trycompile is called without a directory, it must be called from within a file.'); end
    % if neither idents nor dir given, use the calling function's identifier
    if isempty(o.idents)
        o.idents = name; end
end

% uniformize ident format
if isa(o.idents,'function_handle')
    o.idents = char(o.idents); end
if ischar(o.idents)
    o.idents = {o.idents}; end
if isempty(o.idents)
    o.idents = {}; end

% decide whether a re-check can be skipped based on identifiers and directory
if strcmp(o.style,'lazy') || isdeployed
    str = java.lang.String(sprintf('%s:',o.dir,o.idents{:}));
    tag = sprintf('x%.0f',str.hashCode()+3^31);
    if isfield(results,tag)
        ok = results.(tag);
        return;
    end
end

% uniformize directory format
o.dir = path_normalize(o.dir);

% verify style
if ~any(strcmp(o.style,{'force','eager','lazy'}))
    error('Unsupported style: %s',o.style); end

% uniformize test condition
if isa(o.test,'function_handle')
    o.test = char(o.test); end

% uniformize user messages
if ischar(o.errmsg)
    o.errmsg = {o.errmsg}; end
if ischar(o.prepmsg)
    o.prepmsg = {o.prepmsg}; end

% uniformize types
if ischar(o.types)
    o.types = {o.types}; end
for t=1:length(o.types)
    if o.types{t}(1) ~= '*'
        o.types{t} = ['*' o.types{t}]; end
end

% uniformize compiler inputs
if ischar(o.support)
    o.support = {o.support}; end
if ischar(o.includedirs)
    o.includedirs = {o.includedirs}; end
if ischar(o.libdirs)
    o.libdirs = {o.libdirs}; end
if ischar(o.libs)
    o.libs = {o.libs}; end
if ischar(o.defines)
    o.defines = {o.defines}; end
if ischar(o.args)
    o.args = {o.args}; end
for i=1:length(o.support)
    o.support{i} = path_normalize(o.support{i}); end
for i=1:length(o.includedirs)
    o.includedirs{i} = path_normalize(o.includedirs{i}); end
for i=1:length(o.libdirs)
    o.libdirs{i} = path_normalize(o.libdirs{i}); end

% if a support is given as '*'
starred = strcmp('*',o.support);
if any(starred)
    % list all in the given directory that are not mex files
    infos = [];
    for t = 1:length(o.types)
        if ~isempty(infos)
            infos = [infos; dir([o.dir filesep o.types{t}])]; 
        else
            infos = dir([o.dir filesep o.types{t}]); 
        end
    end
    fnames = {infos.name};
    supportfiles = ~cellfun(@(n)is_primary([o.dir filesep n]),fnames);
    o.support = [fnames(supportfiles) o.support(~starred)];
end

% infer directory, if not given (take it from the calling function)
if isempty(o.idents) && ~isempty(o.dir)
    % get all the source files in the given direcctory
    infos = [];
    for t = 1:length(o.types)
        if ~isempty(infos)
            infos = [infos; dir([o.dir filesep o.types{t}])]; 
        else
            infos = dir([o.dir filesep o.types{t}]); 
        end
    end
    fnames = {infos.name};
    % ... but exclude the support files
    fnames = setdiff(fnames,o.support);
    % and apply any renamings to get the corresponding identifiers
    if ~isempty(o.renaming)
        for n=1:length(fnames)
            fnames{n} = hlp_rewrite(fnames{n},o.renaming{:}); end
    end
    % ... and strip off the extensions
    for n=1:length(fnames)
        fnames{n} = hlp_getresult(2,@fileparts,fnames{n}); end
    o.idents = fnames;
end


ok = false;
missingid = []; % indices of missing identifiers (for dont-retry-next-time beacon files)
if isdeployed
    % --- deployed mode ---
    % Can not compile, but figure out whether everything needed is present. A special consideration
    % is that both the mex files calling functions are in a mounted .ctf archive.
    
    % check if all identifiers are present (either as mex or class)
    for i=1:length(o.idents)
        if ~any(exist(o.idents{i}) == [2 3 8])
            missingid(end+1) = i; end
    end
    ok = isempty(missingid);
        
    if ~isempty(missingid)
        % not all identifiers are compiled for this platform
        disp(['Note: The MEX functions/identifiers ' format_cellstr(o.idents(missingid)) ' are not included for your platform.']);
    elseif strcmp(o.style,'force')
        % in force mode, we remark that everything is already compiled
        disp_once(['The functions ' format_cellstr(o.idents) ' are properly compiled.']);
    end
    
else
    % --- regular mode ---
    % here, we *do* compile what needs to be compiled
    
    % find out a key configuration settings
    is64bit = ~isempty(strfind(computer,'64'));
    has_largearrays = is64bit && hlp_matlab_version >= 703;
    
    if ispc
        warning off MATLAB:FILEATTRIB:SyntaxWarning; end
    
    % add a few missing defines
    if hlp_matlab_version < 703
        o.defines = [o.defines {'mwIndex=int','mwSize=int','mwSignedIndex=int'}]; end
    
    % rewrite blas & lapack libs....
    if o.mwlibs
        % for each type of library
        for l={'blas','lapack'}
            lib = l{1};
            % note: this code is based on SeDuMi's compile script (by Michael C. Grant)
            in_use = strcmp(o.libs,lib);
            if any(in_use)
                if ispc
                    if is64bit
                        osdir = 'win64';
                    else
                        osdir = 'win32';
                    end
                    libpath = [matlabroot '\extern\lib\' osdir '\microsoft\libmw' lib '.lib'];
                    if ~exist(libpath,'file')
                        libpath = [matlabroot '\extern\lib\' osdir '\microsoft\msvc60\libmw' lib '.lib']; end
                    if exist(libpath,'file')
                        o.libs{in_use} = libpath; 
                    else
                        disp_once('Note: The Mathworks library %s was assumed to be in %s, but not found.',lib,libpath);
                    end
                else
                    o.libs{in_use} = ['mw' lib];
                end
            end
        end
    end
    
    try
        % remember the current directory & enter the target directory
        olddir = pwd;
        if hlp_matlab_version >= 706
            go_back = onCleanup(@()cd(olddir)); end
        if ~exist(o.dir,'dir')
            error(['The target directory ' o.dir ' does not exist.']); end
        cd(o.dir);
        
        % expand regex patterns in o.support
        for i=length(o.support):-1:1
            if any(o.support{i} == '*')
                found = dir(o.support{i});
                if ~isempty(found)
                    % ... and splice the results in
                    basepath = fileparts(o.support{i});
                    items = cellfun(@(x)[basepath filesep x],{found.name},'UniformOutput',false);
                    o.support = [o.support(1:i-1) items o.support(i+1:end)];
                end
            end
        end
        
        % find all source & target files for the respective identifiers...
        % (note that there might be multiple source files for each one)
        sources = cell(1,length(o.idents)); % list of all source file names for the corresponding identifiers (indexed like idents)
        targets = cell(1,length(o.idents)); % list of all target file names for the corresponding identifiers (indexed like idents)
        for i=1:length(o.idents)
            if ~isempty(o.renaming)
                % the renaming may yield additional source file names for the given identifiers
                idx = strcmp(o.idents{i},o.renaming(2:2:end));
                if any(idx)
                    % the identifier is a renaming target: add the corresponding source file name
                    filename = o.renaming{find(idx)*2-1};
                    % if a source file with this ident & type is present
                    if exist([o.dir filesep filename],'file')
                        % remember it & derive its respective target file name
                        sources{i}{end+1} = filename;
                        targets{i}{end+1} = [o.idents{i} '.' mexext];
                    end
                end
            end
            for t=1:length(o.types)
                filename = [o.idents{i} o.types{t}(2:end)];
                % if a source file with this ident & type is present
                if exist([o.dir filesep filename],'file')
                    % remember it
                    sources{i}{end+1} = filename;
                    % and also derive its respective target file name
                    if strcmp(o.types{t},'*.java')
                        targets{i}{end+1} = [o.idents{i} '.class'];
                    else
                        targets{i}{end+1} = [o.idents{i} '.' mexext];
                    end
                end
            end
            % check whether we have all necessary source files
            if isempty(sources{i})
                error('Did not find source file for %s',o.idents{i}); end
            if isempty(targets{i})
                error('Could not determine target file for %s',o.idents{i}); end
        end
        
        % check for existence (either .mex* or class) of all identifiers
        % and make a list of missing & present binary files; do this in a different directory,
        % to not shadow the mex files of interest with whatever .m files live in this directory
        cd ..
        binaries = {}; % table of existing binary file paths (indexed like idents)
        for i=1:length(o.idents)
            % get current file reference to this identifier
            binaries{i} = which(o.idents{i});
            % if it doesn't point to a .mex or .class file, ignore it
            if ~any(exist(o.idents{i}) == [3 8])
                binaries{i} = []; end
            if ~isempty(binaries{i})
                % check whether it is correct file path
                if ~any(binaries{i} == filesep)
                    error(['Could not determine the location of the mex file for: ' o.idents{i}]); end
                % check whether the referenced file actually exists in the file system
                if isempty(dir(binaries{i}))
                    binaries{i} = []; end
            end
            % if no binary found, record it as missing
            if isempty(binaries{i})
                missingid(end+1) = i; end
        end
        cd(o.dir);
        
        % check which of the existing binaries need to be re-compiled (if out of date)
        outdatedid = []; % indices of identifiers (in o.idents) that need to be recompiled
        for i=1:length(binaries)
            if ~isempty(binaries{i})
                % get the date of the binary file
                bininfo = dir(binaries{i});
                % find all corresponding source files
                srcinfo = [];
                for s=1:length(sources{i})
                    srcinfo = [srcinfo; dir([o.dir filesep sources{i}{s}])]; end
                if ~isfield(bininfo,'datenum')
                    [bininfo.datenum] = celldeal(cellfun(@datenum,{bininfo.date},'UniformOutput',false)); end
                if ~isfield(srcinfo,'datenum')
                    [srcinfo.datenum] = celldeal(cellfun(@datenum,{srcinfo.date},'UniformOutput',false)); end
                % if any of the source files has been changed
                if bininfo.datenum < max([srcinfo.datenum])
                    % check if their md5 hash is still the same...
                    if exist([binaries{i} '.md5'],'file')
                        try
                            contents = load([binaries{i} '.md5'],'-mat','srchash');
                            % need to do that over all source files...
                            srchash = [];
                            sorted_sources = sort(sources{i});
                            for s=1:length(sorted_sources)
                                srchash = [srchash hlp_cryptohash([o.dir filesep sorted_sources{s}],true)]; end
                            if ~isequal(srchash,contents.srchash)
                                % hash is different: mark binary as outdated
                                outdatedid(end+1) = i;
                            end
                        catch
                            % there was a probblem: mark as outdated
                            outdatedid(end+1) = i;
                        end
                    else
                        % no md5 present: mark as outdated
                        outdatedid(end+1) = i; 
                    end
                end
            end
        end
        
        % we try to recompile both what's missing and what's outdated
        recompileid = [missingid outdatedid];
        javainvolved = false; % for final error/help message generation
        mexinvolved = false;  % same
        
        if ~isempty(recompileid)
            % need to recompile something -- display a few preparatory messages...
            for l=1:length(o.prepmsg)
                disp(o.prepmsg{l}); end
            
            failedid = []; % list of indices of identifier that failed the build
            
            % for each identifier, try to compile it
            % and record whether it failed
            for i=recompileid
                success = false;
                fprintf(['Compiling the function/class ' o.idents{i} '...']);
                % for each source file mapping to that identifier
                for s=1:length(sources{i})
                    % check type of source
                    if ~isempty(strfind(sources{i}{s},'.java'))
                        % we have a Java source file: compile
                        [errcode,result] = system(['javac ' javac_options(o) sources{i}{s}]);
                        if errcode
                            javainvolved = true;
                            % problem: show display output
                            fprintf('\n');
                            disp(result);
                        else
                            success = true;
                            break;
                        end
                    else
                        
                        % generate MEX options
                        opts = mex_options(o);
                        supp = sprintf(' %s',o.support{:});
                        
                        if isempty(compiler_selected)
                            % not clear whether a compiler has been selected yet
                            if hlp_matlab_version >= 708
                                % we can find it out programmatically
                                try
                                    cconf = mex.getCompilerConfigurations; %#ok<NASGU>
                                    compiler_selected = true;
                                catch
                                    % no compiler has been selected yet...
                                    try
                                        % display a few useful hints to the user
                                        disp(' to compile this feature, you first need to select');
                                        disp('which compiler should be used on your platform.');
                                        if ispc
                                            if is64bit
                                                disp_once('As you are on 64-bit windows, you may find that no compiler is installed.');
                                            else
                                                disp_once('On 32-bit Windows, MATLAB supplies a built-in compiler (LLC), which should');
                                                disp_once('faithfully compile most C code. For broader support across C dialects (as well as C++), ');
                                                disp_once('you should make sure that a better compiler is installed on your system and selected in the following.');
                                            end
                                            disp_once('A good choice is the free Microsoft Visual Studio 2005/2008/2010 Express compiler suite');
                                            disp_once('together with the Microsoft Platform SDK (6.1 for 2008, 7.1 for 2010) for your Windows Version.');
                                            disp_once('See also: http://argus-home.coas.oregonstate.edu/forums/development/core-software-development/compiling-64-bit-mex-files');
                                            disp_once('          http://www.mathworks.com/support/compilers/R2010b/win64.html');
                                            disp_once('The installation is easier if a professional Intel or Microsoft compiler is used.');
                                        elseif isunix
                                            disp_once('On Linux/UNIX, the best choice is usually a supported version of the GCC compiler suite.');
                                        else
                                            disp_once('On Mac OS, you need to have a supported version of Xcode/GCC installed.');
                                        end
                                        % start the compiler selection tool
                                        mex -setup
                                        % verify that a compiler has been selected
                                        cconf = mex.getCompilerConfigurations; %#ok<NASGU>
                                        compiler_selected = true;
                                    catch
                                        compiler_selected = false;
                                    end
                                end
                            else
                                disp(' you may be prompted to select a compiler in the following');
                                disp('(as BCILAB cannot auto-determine whether one is selected on your platform).');
                            end
                        end
                        
                        if ~compiler_selected
                            fprintf('skipped (no compiler selected).\n');
                        else
                            
                            if o.verbose || isempty(compiler_selected)
                                % this variant will also be brought up if not sure whether a compiler
                                % has already been selected...
                                doeval = @eval;
                            else
                                doeval = @evalc;
                            end
                            
                            % try to build the file
                            try
                                mexinvolved = true;
                                % check if a renaming applies...
                                idx = strcmp(sources{i}{s},o.renaming);
                                if any(idx)
                                    rename = [' -output ' o.renaming{find(idx,1)+1} ' '];
                                else
                                    rename = '';
                                end
                                if o.debugcompile
                                    % display a debug console to allow the user to debug how their file compiles
                                    fprintf('\nExecution has been paused immediately before running mex.\n');
                                    disp(['You are in the directory "' pwd '".']);
                                    disp('The mex command that will be invoked in the following is:');
                                    if has_largearrays
                                        disp([' mex ' opts ' -largeArrayDims ' rename sources{i}{s} supp]);
                                    else
                                        disp([' mex ' opts rename sources{i}{s} supp]);
                                    end
                                    fprintf('\n\nTo proceed normally, type "dbcont".\n');
                                    keyboard;
                                end
                                if has_largearrays
                                    try
                                        % -largeArrayDims enabled
                                        try
                                            doeval(['mex ' opts ' -largeArrayDims ' rename sources{i}{s} supp]); % with supporting libaries
                                        catch
                                            doeval(['mex ' opts ' -largeArrayDims ' rename sources{i}{s}]); % without supporting libraries
                                        end
                                    catch
                                        % -largeArrayDims disabled
                                        try
                                            doeval(['mex ' opts rename sources{i}{s} supp]); % with supporting libaries
                                        catch
                                            doeval(['mex ' opts rename sources{i}{s}]); % without supporting libraries
                                        end
                                    end
                                else
                                    % -largeArrayDims disabled
                                    try
                                        doeval(['mex ' opts rename sources{i}{s} supp]); % with supporting libaries
                                    catch
                                        doeval(['mex ' opts rename sources{i}{s}]); % without supporting libraries
                                    end
                                end
                                
                                % compilation succeeded...
                                if any(i==outdatedid)
                                    % there is an outdated binary, which needs to be deleted
                                    try
                                        delete(binaries{i});
                                    catch
                                        disp(['Could not delete outdated binary ' binaries{i}]);
                                    end
                                end
                                % check whether the file is being found now
                                if exist(o.idents{i}) == 3
                                    success = true;
                                    compiler_selected = true;
                                    break;
                                end
                            catch
                                % build failed
                            end
                        end
                    end
                end
                
                
                % check if compilation of this identifier was successful
                if success
                    % if so, we sign off the binary with an md5 hash of the sources...
                    newbinary = which(o.idents{i});
                    try
                        srchash = [];
                        sorted_sources = sort(sources{i});
                        for s=1:length(sorted_sources)
                            srchash = [srchash hlp_cryptohash([o.dir filesep sorted_sources{s}],true)]; end
                        save([newbinary '.md5'],'srchash','-mat');
                        fprintf('success.\n');
                    catch
                        disp('could not create md5 hash for the source files; other than that, successful.');
                    end
                else
                    fprintf('failed.\n');
                    failedid(end+1) = i;
                end
            end
            
            embed_test = false;
            if isempty(failedid)
                % all worked: now run the test code - if any - to verify the correctness of the build
                if length(recompileid) > 1
                    if ~isempty(o.test)
                        fprintf('All files in %s compiled successfully; now testing the build outputs...',o.dir);
                    else
                        fprintf('All files in %s compiled successfully.\n',o.dir);
                    end
                elseif ~isempty(o.test)
                    fprintf('Now testing the build outputs...');
                end
                
                % test the output
                try
                    ans = true; %#ok<NOANS>
                    eval(o.test);
                catch
                    ans = false; %#ok<NOANS>
                end
                if ans %#ok<NOANS>
                    % the test was successful; now copy the files into a platform-specific directory                    
                    if ~isempty(o.test)
                        % only if we have a succeeding non-empty test
                        embed_test = true;
                        fprintf('success.\n'); 
                    end
                    retainid = recompileid;
                    eraseid = [];
                    ok = true;
                else
                    % the test was unsuccessful: remove all newly-compiled files...
                    if ~isempty(o.test)
                        fprintf('failed.\n'); end
                    disp('The code compiled correctly but failed the build tests. Reverting the build...');
                    disp('If this is unmodified BCILAB code, please consider reporting this issue.');
                    retainid = [];
                    eraseid = recompileid;
                end
            else
                if length(recompileid) > 1
                    if isempty(setdiff(recompileid,failedid))
                        fprintf('All files in %s failed to build; this indicates a problem in your build environment/settings.\n',o.dir);
                    else
                        fprintf('Some files in %s failed to build. Please make sure that you have a supported compiler; otherwise, please report this issue.\n',o.dir);
                    end
                else
                    disp('Please make sure that you have a supported compiler and that your build environment is set up correctly.');
                    disp('Also, please consider reporting this issue.');
                end
                % compilation failed; only a part of the binaries may be available...
                retainid = setdiff(recompileid,failedid);
                eraseid = [];
            end
            
            % move the mex files into their own directory
            moveid = retainid(cellfun(@exist,o.idents(retainid)) == 3);
            if ~isempty(moveid)
                % some files to be moved
                dest_path = [o.dir filesep 'build-' hlp_hostname filesep];
                % create a new directory
                if ~exist(dest_path,'dir')
                    if ~mkdir(o.dir,['build-' hlp_hostname])
                        error(['unable to create directory ' dest_path]); end
                    % set permissions
                    try
                        fileattrib(dest_path,'+w','a');
                    catch
                        disp(['Note: There are permission problems for the directory ' dest_path]);
                    end
                end
                % create a new env_add.m there
                try
                    filename = [dest_path 'env_add.m'];
                    fid = fopen(filename,'w+');
                    if embed_test
                        % if we had a successful test, we use this to control inclusion of the mex files
                        fprintf(fid,o.test);
                    else
                        % otherwise we check whether any one of the identifiers is recognized by 
                        % MATLAB as a mex function
                        fprintf(fid,'any(cellfun(@exist,%s)==3)',hlp_tostring(o.idents(moveid)));
                    end
                    fclose(fid);
                    fileattrib(filename,'+w','a');
                catch
                    disp(['Note: There were write permission problems for the file ' filename]);
                end                
                % move the targets over there...
                movefiles = unique(o.idents(moveid));
                for t = 1:length(movefiles)
                    [d,n,x] = fileparts(which(movefiles{t}));
                    movefile([d filesep n x],[dest_path n x]);
                    try
                        movefile([d filesep n x '.md5'],[dest_path n x '.md5']);
                    catch
                    end
                end
                % add the destination path
                addpath(dest_path);
                % and to be entirely sure, CD into that directory to verify that the files are being recognized...
                % (and don't get shadowed by whatever is in the directory below)
                cd(dest_path);
                all_ok = all(strncmp(dest_path,cellfun(@which,o.idents(moveid),'UniformOutput',false),length(dest_path)));
                cd(o.dir);
                % make sure that they are still being found...
                if ~all_ok
                	error('It could not be verified that the MEX file records in %s were successfully updated to their new sub-directories.',o.dir); end
            end

            % move the java class files into their own directory
            infos = dir([o.dir filesep '*.class']);
            movefiles = {infos.name};
            moveid = retainid(cellfun(@exist,o.idents(retainid)) ~= 3);
            if ~isempty(movefiles)
                % some files to be moved
                dest_path = [o.dir filesep 'build-javaclasses' filesep];
                % create a new directory
                if ~exist(dest_path,'dir')
                    if ~mkdir(o.dir,'build-javaclasses')
                        error(['unable to create directory ' dest_path]); end
                    % set permissions
                    try
                        fileattrib(dest_path,'+w','a');
                    catch
                        disp(['Note: There are permission problems for the directory ' dest_path]);
                    end
                end
                % create a new env_add.m there
                try
                    filename = [dest_path 'env_add.m'];
                    fid = fopen(filename,'w+'); fclose(fid);
                    fileattrib(filename,'+w','a');
                catch
                    disp(['Note: There were write permission problems for the file ' filename]);
                end                
                % move the targets over there...
                for t = 1:length(movefiles)
                    movefile([o.dir filesep movefiles{t}],[dest_path movefiles{t}]); 
                    try
                        movefile([o.dir filesep movefiles{t} '.md5'],[dest_path movefiles{t} '.md5']); 
                    catch
                    end                    
                end
                % add the destination path
                if isdeployed
                    warning off MATLAB:javaclasspath:jarAlreadySpecified; end
                javaaddpath(dest_path);
                
                % check whether the class is found
                if ~all(cellfun(@exist,o.idents(moveid)) == 8)
                    disp_once('Not all Java binaries in %s could be recognized by MATLAB.',dest_path); end
            end
            
            
            if ~isempty(eraseid)
                % some files need to be erased...
                for k=eraseid
                    for t=1:length(targets{k})
                        if exist([o.dir filesep targets{k}{t}])
                            try
                                delete(targets{k}{t});
                            catch
                                disp(['Could not delete broken binary ' binaries{i}{s}]);
                            end
                        end
                    end
                end
            end
        else
            % nothing to recompile
            ok = true;
            if strcmp(o.style,'eager') && ~isempty(o.idents) && o.verbose
                disp_once(['The functions ' format_cellstr(o.idents) ' are already compiled.']); end
        end
        
        % go back to the old directory
        if hlp_matlab_version < 706
            cd(olddir); end
    catch e
        ok = false; %#ok<NASGU>
        % go back to the old path in case of an error
        if hlp_matlab_version < 706
            cd(olddir); end
        rethrow(e);
    end
end

% store the OK flag in the results
if strcmp(o.style,'lazy') || isdeployed
    results.(tag) = ok; end

if ~ok
    if ~isdeployed
        % regular error summary
        if mexinvolved
            disp_once('\nIn case you need to use a better / fully supported compiler, please have a look at:');
            try
                v=version;
                releasename = v(find(v=='(')+1 : find(v==')')-1);
                if length(releasename) > 3 && releasename(1) == 'R'
                    releasename = releasename(2:end); end
                site = ['http://www.google.com/search?q=matlab+supported+compilers+' releasename];
            catch
                site = 'http://www.google.com/search?q=matlab+supported+compilers';
            end
            disp_once('  <a href="%s">%s</a>\n',site,site);
            if ispc
                if is64bit
                    disp_once('On 64-bit Windows, MATLAB comes with no built-in compiler, so you need to have one installed.');
                else
                    disp_once('On 32-bit Windows, MATLAB supplies a built-in compiler (LLC), which is, however, not very good.');
                end
                disp_once('A good choice is the free Microsoft Visual Studio 2005/2008/2010 Express compiler suite');
                disp_once('together with the Microsoft Platform SDK (6.1 for 2008, 7.1 for 2010) for your Windows Version.');
                disp_once('See also: http://argus-home.coas.oregonstate.edu/forums/development/core-software-development/compiling-64-bit-mex-files');
                disp_once('          http://www.mathworks.com/support/compilers/R2010b/win64.html');
                disp_once('The installation is easier if a professional Intel or Microsoft compiler is used.');
            elseif isunix
                disp_once('On Linux/UNIX, the best choice is usually a supported version of the GCC compiler suite.');
            else
                disp_once('On Mac OS, you need to have a supported version of Xcode/GCC installed.');
            end
        end
        if javainvolved
            disp_once('Please make sure that your system''s java configuration matches the one used by MATLAB (see "ver" command).');
        end
    end
end



% create javac options string from options struct
function opts = javac_options(o)
verbosity = hlp_rewrite(o.verbose,true,'-verbose',false,'');
if ~isempty(o.libdirs)
    cpath = ['-classpath ' sprintf('%s;',o.libdirs{:})];
    cpath(end) = [];
else
    cpath = '';
end
if ~isempty(o.includedirs)
    ipath = ['-sourcepath ' sprintf('%s;',o.includedirs{:})];
    ipath(end) = [];
else
    ipath = '';
end
debugness = hlp_rewrite(o.debug,true,'-g',false,'-g:none');
targetsource = '-target 1.6 -source 1.6';
opts = [sprintf(' %s',verbosity,cpath,ipath,debugness,targetsource) ' '];


% create mex options string from options struct
function opts = mex_options(o)
if ~isempty(o.defines)
    defs = sprintf(' -D%s',o.defines{:});
else
    defs = '';
end
debugness = hlp_rewrite(o.debug,true,'-g',false,'');
if ~isempty(o.includedirs)
    incdirs = sprintf(' -I"%s"',o.includedirs{:});
else
    incdirs = '';
end
if ~isempty(o.libs)
    if ispc
        libs = sprintf(' -l"%s"',o.libs{:});
    else
        libs = sprintf(' -l%s',o.libs{:});
    end
else
    libs = '';
end
if ~isempty(o.libdirs)
    libdirs = sprintf(' -L"%s"',o.libdirs{:});
else
    libdirs = '';
end
verbosity = hlp_rewrite(o.verbose,true,'-v',false,'');
opts = [sprintf(' %s',defs,debugness,incdirs,libs,libdirs,verbosity) ' '];


% format a non-empty cell-string array into a string
function x = format_cellstr(x)
if isempty(x)
    x = '';
else
    x = ['{' sprintf('%s, ',x{1:end-1}) x{end} '}'];
end


% check whether a given identifier is frozen in a ctf archive
function tf = in_ctf(ident) %#ok<DEFNU>
tf = isdeployed && strncmp(ctfroot,which(ident),length(ctfroot));


% normalize a directory path
function dir = path_normalize(dir)
if filesep == '\';
    dir(dir == '/') = filesep;
else
    dir(dir == '\') = filesep;
end
if dir(end) == filesep
    dir = dir(1:end-1); end


% determine if a given file is a mex source file or a java source file
% (and compiles into an identifier that is seen by MATLAB)
function tf = is_primary(filename)
if length(filename)>5 && strcmp(filename(end-4:end),'.java')
    tf = true;
    return; 
else
    tf = false;
end
fid = fopen(filename);
if fid ~= -1
    try
        contents = fread(fid);
        tf = ~isempty(strfind(char(contents)','mexFunction')); %#ok<FREAD>
        fclose(fid);
    catch
        fclose(fid);
    end
end


% act like deal, but with a single cell array as input
function varargout = celldeal(argin)
varargout = argin;


% for old MATLABs that can't properly move files...
function movefile(src,dst)
try
    builtin('movefile',src,dst);    
catch e
    if any([src dst]=='$') && hlp_matlab_version <= 705
        if ispc
            [errcode,text] = system(sprintf('move ''%s'' ''%s''',src,dst)); %#ok<NASGU>
        else
            [errcode,text] = system(sprintf('mv ''%s'' ''%s''',src,dst)); %#ok<NASGU>
        end
        if errcode
            error('Failed to move %s to %s.',src,dst); end
    else
        rethrow(e);
    end
end


% for old MATLABs that don't handle Java classes on the dynamic path...
function res = exist(obj,type)
if nargin > 1
    res = builtin('exist',obj,type);
    if  ~res && (hlp_matlab_version <= 704) && strcmp(type,'class') && builtin('exist',[obj '.class'],'file')
        res = 8; end
else
    res = builtin('exist',obj);
    if  ~res && (hlp_matlab_version <= 704) && builtin('exist',[obj '.class'])
        res = 8; end
end


% for old MATLABs that don't handle Java classes on the dynamic path...
function res = which(ident)
res = builtin('which',ident);
if ~any(res == filesep)
    if hlp_matlab_version <= 704
        if isempty(res) || ~isempty(strfind(res,'not found')) 
            res = builtin('which',[ident '.class']); end
    else
        if ~isempty(strfind(res,'Java')) 
            res = builtin('which',[ident '.class']); end
    end
end
