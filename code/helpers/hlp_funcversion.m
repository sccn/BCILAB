function v = hlp_funcversion(func,filename,versionformat)
% Get the version string of a MATLAB function, or an MD5 hash if unversioned.
% Version = hlp_funcversion(Function,Filename,VersionFormat)
%
% Returns just the string form of the input if a version cannot be determined
% or if the versionformat is passed in as false.
%
% In:
%   Function : name or handle to the function
%
%   Filename : optionally the filename of the function (can disambiguate cases where only the
%              identifier is known)
%
%   VersionFormat : regular expression of the version string inside the function's
%                   code. If containing the substring '$(funcname)', it will be 
%                   substituted by the name of the function. (default:'$(funcname)_version<\S+>')
%
% Out:
%   Version : Version string for the function, or MD5 hash of the code if no 
%             version string present. In deployed mode this function will return the MD5 hash
%             of the compiled code, unless the filename of the uncompiled code is provided.
%
%                                 Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                 2013-04-19

if nargin < 2
    if ischar(func)
        filename = which(func);
    else
        filename = getfield(functions(func),'file');
    end
end
if nargin < 3
    versionformat = '$(funcname)_version<\S+>'; end

try
    func = char(func);
    if ~isempty(filename)
        % open the source file
        f = fopen(filename,'r');
        try
            % read the code
            code = fread(f,Inf,'uint8=>char')';
            % check if it contains the version descriptor tag
            v = regexp(code,strrep(versionformat,'$(funcname)',func),'match');
            % otherwise we just hash the entire code
            if isempty(versionformat) || isempty(v)
                v = hlp_cryptohash(code); end
            if iscell(v)
                if length(v) > 1
                    warn_once('BCILAB:hlp_funcversion:multiple_tags','There is more than one version tag in file %s. Using the first one.',filename); end
                v = v{1}; 
            end
            fclose(f);
        catch %#ok<CTCH>
            try
                fclose(f);
            catch %#ok<CTCH>
            end
            v = func;
        end
    else
        % otherwise use the string representation as version
        v = func;
    end
catch
    v = func;    
end