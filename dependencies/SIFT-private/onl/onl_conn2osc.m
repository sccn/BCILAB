function onl_conn2osc(varargin)
% write connectivity contents to a TCP port via Open Sound Control (OSC)
% based on run_writeosc.m from BCILAB (www.sccn.ucsd.edu/wiki/bcilab)
    
persistent conn;

% TODO: this function will be called by vis_causalBrainMovie3D. There will
% be an additional option in that vis_causalBrainMovie3D to pass data on via 
% OSC

% extract some stuff from inputs for arg defaults
% Conn = arg_extract(varargin,{'EEG','ALLEEG'},1);
% 
% if ~isempty(Conn)
%     Conn = EEG.CAT.Conn;
%     ConnNames   = hlp_getConnMethodNames(Conn);
%     conndef     = ConnNames{1};
%     freqrange   = [Conn.freqs(1) Conn.freqs(end)];
%     freqdef     = Conn.freqs; %['[' num2str(freqrange(1)) ':' num2str(Conn.freqs(2)-Conn.freqs(1)) ':' num2str(freqrange(end)) ']'];
%     timerange   = [Conn.erWinCenterTimes(1) Conn.erWinCenterTimes(end)];
%     timedef     = [timerange(1) timerange(end)];
%     clear Conn;
% else
%     ConnNames = {''};
%     conndef = '';
%     [freqrange, freqdef, timerange, timedef] = deal([]);
% end


% define arguments
arg_define(varargin, ...
    arg_norep({'EEG','ALLEEG'},mandatory,[],'EEGLAB dataset. Must contain Conn structure'), ...
    arg({'estimator','Estimator'},ConnNames{1},ConnNames,'Estimator to visualize','shape','row'), ...
    arg_subtoggle({'new_connection','NewConnection'},{},...
    {...
        arg({'out_hostname','OutputHost'}, '127.0.0.1',[],'Destination TCP hostname. Can be a computer name, URL, or IP address.'), ...
        arg({'out_port','OutputPort'}, 12346, [],'Destination TCP port. Depends on the destination application, usually > 1024.'), ...
        arg({'out_path','OutputPath'},'/sift',[],'Output OSC address. This is the local name of the OSC method that should be invoked (begins with a /)'), ...
        arg({'msg_format','MessageFormatter'},'@(D) num2cell(single(D(:))))',[],'Message Formatting function. Converts a BCI output (usually scalar or row vector, in some complex cases a matrix) into a cell array of OSC-compatible data types (float32,int32, ...).','type','expression'), ...
    },'create a new OSC connection'), ...
    arg_nogui({'osc_conn','OSC_Connection'},[],[],'Existing OSC connection. If specified, data will be written to this object.') ...
    );

% check the dataset
res = hlp_checkeegset(ALLEEG,{'conn'});
if ~isempty(res)
    error(['SIFT:' fcnName],res{1});
end

% convert format strings to formatting functions
if ischar(msg_format)
    msg_format = @(D) sprintf(msg_format,D); end

% connect to remote address via OSC
if ~exist('osc_new_address','file')
    try
        build_osc;
    catch
        error('The OSC library has not been built for your platform yet; see dependencies/OSC* for more info.'); 
    end
end

if new_connection.arg_selection && isempty(osc_conn)
    conn = osc_new_address(out_hostname,out_port);
elseif ~isempty(osc_conn)
    if ~isempty(conn)
        % free existing connection...
        osc_free_address(conn);
    end
    % ... and switch to new one
    conn = osc_conn;
elseif isempty(conn)
    error('A valid OSC connection object could not be found. Please create a new connection or provide an existing one');
end

% deleter = onCleanup(@()osc_free_address(conn)); % if the last reference to this is dropped, the connection is closed (on MATLAB 2008a+)

% prepare OSC object
connmat = EEG.CAT.Conn.(estimator);



% start background writer job
onl_write_background(@(y)send_message(y,conn,out_path,msg_format,deleter),in_stream,pred_model,out_form,update_freq,0,pred_name);

% background message sending function
function send_message(y,conn,opath,formatter,deleter)
msg = struct('path',opath,'data',{formatter(y)});
if osc_send(conn,msg) == 0
    error('OSC transmission failed.'); end
