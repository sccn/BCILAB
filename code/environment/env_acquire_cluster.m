function env_acquire_cluster
% Acquire a cluster using the currently configured acquire options.
%
% The behavior of this function is governed by the configuration variable acquire_options that can
% be set either in the config file or in the "Cluster Settings" GUI dialog. Generally, this function
% will (if supported by the OS), check which of the desired cluster resources are already up, and
% start what still needs to be started. There is no guarantee that the cluster resources actually
% start successfully (and don't crash, etc), the function just tries its best. After that, the
% function will start a "heartbeat" timer that periodically tells the cluster resources that it is
% still interested in keeping them alive. Depending on how the cluster was configured, it may shut
% down the resources after they are no longer needed (e.g., for cost reasons). 
%
% The heartbeat signal can be disabled by calling env_release_cluster or pressing the "request
% cluster availability" button in the GUI toolbar for a second time.
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2012-04-12

global tracking;
import java.io.*
import java.net.*
import java.lang.*

% interpret as cell array of arguments to par_getworkers_ssh
arguments = tracking.acquire_options;
if isempty(arguments)
    disp('No settings for acquiring the cluster have been specified, no startup command will be issued.'); 
elseif ~iscell(arguments)
    disp('The acquire_options parameter should be a cell array of name-value pairs which form arguments to par_getworkers_ssh.');
elseif isunix
    % invoke, but also impose default arguments (if unspecified)
    % by default, workers do not recruit (but only list) other workers, preventing a cascading effect
    try
        par_getworkers_ssh(arguments{:});
    catch e
        disp('Could not acquire worker machines; traceback: ');
        env_handleerror(e);
    end
else
    disp('Cannot automatically acquire hosts from a non-UNIX system.');
end

% start the heartbeat timer
if ~isempty(tracking.parallel.pool)
    fprintf('Initiating heartbeat signal... ');
    tracking.cluster_requested = {};
    % for each endpoint in the pool...
    for p=1:length(tracking.parallel.pool)
        % make a new socket
        sock = DatagramSocket();
        % and "connect" it to the worker endpoint (its heartbeat server)
        endpoint = hlp_split(tracking.parallel.pool{p},':');
        sock.connect(InetSocketAddress(endpoint{1}, str2num(endpoint{2})));
        tracking.cluster_requested{p} = sock;
    end
    % start a timer that sends the heartbeat (every 30 seconds)
    start(timer('ExecutionMode','fixedRate', 'Name','heartbeat_timer', 'Period',15, ...
        'TimerFcn',@(timer_handle,varargin)send_heartbeat(timer_handle)));
    disp('success.');
else
    tracking.cluster_requested = true;
end

% called periodically to send heartbeat messages over the network
function send_heartbeat(timer_handle)
import java.io.*
import java.net.*
import java.lang.*
global tracking;
try
    if ~isfield(tracking,'cluster_requested') || isempty(tracking.cluster_requested) || ~iscell(tracking.cluster_requested)
        error('ouch!'); end
    
    % for each socket in the list of heartbeat sockets...
    for p=1:length(tracking.cluster_requested)
        try
            % send the message
            tmp = DatagramPacket(uint8('dsfjk45djf'),10);
            tracking.cluster_requested{p}.send(tmp);
        catch e
            % some socket problem: ignored.
        end
    end
catch e
    % issue (most likely the request has been cancelled) stop & sdelete the heartbeat timer
    stop(timer_handle);
    delete(timer_handle);
    disp('Heartbeat signal stopped.');
end
