function hash = hlp_cryptohash(data,fromfile)
% Compute an MD5 hash of a file, string or generic data structure.
% Hash = hlp_cryptohash(Data,FromFile)
%
% In:
%   Data : data to be hashed; can be a filename, a string, or any other MATLAB data structure.
%
%   FromFile : if true, data is interpreted as a file name (default: false)
%
% Out:
%   Hash : MD5 hash (decimal, lowercase) of the Data
%
% Examples:
%   % calculate an md5 hash of a file
%   hlp_cryptohash('myscript.m',true);
%
%   % calculate an md5 hash of a data structure
%   hlp_cryptohash(lastmodel);
%
% See also:
%   hlp_fingerprint, hlp_serialize
%
%					      	     Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-10-10

% Copyright (C) Christian Kothe, SCCN, 2010, christian@sccn.ucsd.edu
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

% maximum amount of memory used inside the Java VM
max_java_memory = 2^26; % 64 MB

if nargin < 2
    fromfile = false; end    
if ~(isequal(fromfile,true) || isequal(fromfile,false))
    error('The given FromFile argument must be true or false.'); end

if fromfile
    % take data as a file name
    if ~ischar(data)
        error('To represent a file name, Data should be a string.'); end
    if ~exist(data,'file')
        error('The file %s does not exist.',data); end
    f = fopen(data,'r');
    try
        data = fread(f,Inf);
        fclose(f);
    catch e
        try 
            fclose(f);
        catch %#ok<CTCH>
        end
        rethrow(e);
    end
else
    % take data literally
    if ~ischar(data)
        data = hlp_serialize(data); end
end

% use Java to hash the data (idea from Michael Kleder)
hasher = java.security.MessageDigest.getInstance('MD5');
data = uint8(data);
if length(data) <= max_java_memory
    hasher.update(data);
else
    numsplits = ceil(length(data)/max_java_memory);
    for i=0:numsplits-1
        range = 1+floor(i*length(data)/numsplits) : min(length(data),floor((i+1)*length(data)/numsplits));
        hasher.update(data(range));
    end
end

hash = dec2hex(typecast(hasher.digest,'uint8'),2)';
hash = lower(hash(:)');
