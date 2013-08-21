function v = utl_funcversion(func,versionformat)
% Get the version string of a MATLAB function, or an MD5 hash if unversioned.
% Version = utl_funcversion(Function,VersionFormat)
%
% Returns just the string form of the input if a version cannot be determined
% or if the versionformat is passed in as false.
%
% In:
%   Function : name or handle to the function
%
%   VersionFormat : regular expression of the version string inside the function's
%                   code. If containing the substring '$(funcname)', it will be 
%                   substituted by the name of the function. (default:'$(funcname)_version<\S+>')
%
% Out:
%   Version : version string for the function, or MD5 hash of the code if no 
%             version string present.
%
%                                 Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                 2013-04-19

if ~exist('versionformat','var')
    versionformat = '$(funcname)_version<\S+>'; end

func = char(func);
filename = which(func);
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
        fclose(f);
    catch
        try
            fclose(f);
        catch,end
        v = func;
    end
else
    % otherwise use the string representation as version
    v = func;
end
