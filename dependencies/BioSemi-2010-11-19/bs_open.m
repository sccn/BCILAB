function h = bs_open(basepath)
% Open a BioSemi connection
% [Handle] = bs_open(BasePath)
%
% In:
%   BasePath : Optionally the path to the folder in which the Win32/Win64/Linux32/... folders are
%              located (default: this path)
%
% Out:
%   Handle : handle to the connection
%            when the handle is deleted, the connection is automatically closed
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-12-09

persistent uid;
if isempty(uid)
    uid = 1;
else
    uid = uid+1;
end

% find out which DLL path we need
if ~exist('basepath','var') || isempty(basepath)
    basepath = [fileparts(mfilename('fullpath')) filesep]; end
if ispc
    if strfind(computer,'64')
        dllpath = [basepath 'Win64' filesep 'Labview_DLL.dll'];
    else
        dllpath = [basepath 'Win32' filesep 'Labview_DLL.dll'];
    end
elseif ismac
    dllpath = [basepath 'Mac' filesep 'liblabview_dll.0.0.1.dylib'];
elseif isunix
    if strfind(computer,'64')
        dllpath = [basepath 'Linux64' filesep 'liblabview_dll.so'];
    else
        dllpath = [basepath 'Linux32' filesep 'liblabview_dll.so'];
    end
else
    error('Operating System not supported by the BioSemi driver.');
end
        
% open the stream
[h.hDLL,h.hConn,h.hOPEN_DRIVER_ASYNC,h.hUSB_WRITE,h.hREAD_MULTIPLE_SWEEPS,h.hREAD_POINTER,h.hCLOSE_DRIVER_ASYNC,h.pBuffer,h.nbchan,h.srate,h.nbsync,h.nbtrig,h.nbeeg,h.nbexg,h.last_ptr] = bsb_open(dllpath);
% make sure that it gets auto-deleted when the handle is erased
h.cleanup = onCleanup(@()bsb_close(h.hDLL,h.hConn,h.hUSB_WRITE,h.hCLOSE_DRIVER_ASYNC,h.pBuffer));


% add misc fields
if h.srate == 0
    h.srate = 2048; end
h.nbextra = h.nbchan - (h.nbsync+h.nbtrig+h.nbeeg+h.nbexg);
h.uid = uid;

channels = {};
for k=1:double(h.nbsync)
    channels{end+1} = ['Sync' num2str(k)]; end
for k=1:double(h.nbtrig)
    channels{end+1} = ['Trig' num2str(k)]; end
for k=1:double(h.nbeeg)
    letters = 'A':'Z';
    channels{end+1} = [letters(ceil(k/32)) num2str(1+mod(k-1,32))];
end
for k=1:double(h.nbexg)
    channels{end+1} = ['EX' num2str(k)]; end
for k=1:double(h.nbextra)
    channels{end+1} = ['X' num2str(k)]; end
h.channels = channels;
