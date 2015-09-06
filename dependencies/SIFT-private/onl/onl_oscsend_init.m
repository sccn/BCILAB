function osc_conn = onl_oscsend_init(varargin)
% Intialize an Open Sound Control (OSC) connection.
% Returns an OSC connection object.
% 
% % Example: 
% % initialize local connection
% osc_conn = onl_oscsend_init('out_hostname','127.0.0.1','out_port',3334);
% % send struct as an OSC message bundle (one message per field)
% onl_oscsend('data',struct('msg1',[1 2 3],'msg2','test','msg3',logical(0)),'out_path','/mypath','osc_conn',osc_conn);
%
% Author: Tim Mullen, SCCN/INC/UCSD 2013
% 
% See Also: onl_oscsend()

arg_define(varargin, ...
    arg({'out_hostname','OutputHost'}, '127.0.0.1',[],'Destination TCP hostname. Can be a computer name, URL, or IP address.'), ...
    arg({'out_port','OutputPort'}, 12346, [],'Destination TCP port. Depends on the destination application, usually > 1024.') ...
    );

% connect to remote address via OSC
if ~exist('osc_new_address','file')
    try
        build_osc;
    catch
        error('The OSC library has not been built for your platform yet; see dependencies/OSC* for more info.'); 
    end
end

try
    osc_conn = osc_new_address(out_hostname,out_port);
catch err
    hlp_handleerror(err);
    disp('Unable to initialize OSC connection');
end
