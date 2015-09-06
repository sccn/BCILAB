function res = onl_oscsend(varargin)
% Write data to a UDP port via Open Sound Control (OSC)
%
% Example: 
% % initialize local connection
% osc_conn = onl_oscsend_init('out_hostname','127.0.0.1','out_port',3334);
% % send struct as an OSC message bundle (one message per field)
% onl_oscsend('data',struct('msg1',[1 2 3],'msg2','test','msg3',logical(0)),'out_path','/mypath','osc_conn',osc_conn);
%
% Author: Tim Mullen, SCCN/INC/UCSD 2013
%
% See Also: onl_oscsend_init()

% define arguments
arg_define(varargin, ...
    arg_norep({'data','Data'},mandatory,[],'Data to transmit. This is a structure where each field is a separate message to send in a bundle. The field contents can be scalar or 1-dimensional uniform datatype arrays (standard or cell) of the following types: int32,single,logical,char/string'), ...
    arg_subtoggle({'new_connection','NewConnection'},{},@onl_oscsend_init,'Create a new OSC connection'), ...
    arg({'out_path','OutputPath'},'/sift',[],'OSC output address. This is the base OSC routing address (begins with a /)'), ...
    arg({'routefields','RouteFieldNames'},false,[],'Route OSC messages based on field names. This will append to the OutputPath an additional routing suffix for each message based on its respective fieldname.'), ...
    arg_nogui({'osc_conn','OSC_Connection'},[],[],'Existing OSC connection. If specified, data will be written to this object.') ...
    );
% arg({'msg_format','MessageFormatter'},'@(D) num2cell(single(D(:))))',[],'Message Formatting function. Converts input (usually scalar or row vector, in some complex cases a matrix) into a cell array of OSC-compatible data types (float32,int32, ...). Can also be a sprintf-type format string (''%s'',''%d'',...)','type','expression'), ...
% convert format strings to formatting functions
% if ischar(msg_format)
%     msg_format = @(D) sprintf(msg_format,D); end
% deleter = onCleanup(@()osc_free_address(conn)); % if the last reference to this is dropped, the connection is closed (on MATLAB 2008a+)
% msg = structfun(@format_msg,data,'UniformOutput',false);

if new_connection.arg_selection && isempty(osc_conn)
    osc_conn = onl_oscsend_init(new_connection);
end

fn = fieldnames(data);
msg = cell(1,length(fn));
for fi=1:length(fn)
    msg{fi} = format_msg(data.(fn{fi}),fastif(routefields,['/' fn{fi}],''),out_path);
end

% send the message / bundle
res = osc_send(osc_conn,msg);
if res == 0
    error('OSC transmission failed.'); 
end


% Helper function to format each message
function msgo = format_msg(msg,routename,opath)

% convert datatype
[msg dtype] = convert_datatype(msg);
N = length(msg);
if isnumeric(msg) || islogical(msg)
    msg = {num2cell(msg)};
elseif ~iscell(msg)
    msg = {msg};
end

if length(dtype)~=N      
    dtype = repmat(dtype,[1 N]);
end
% convert to OSC struct format
msgo = struct('path',[opath routename],'typetags',dtype);
msgo.data = msg;

function [msg dtype] = convert_datatype(msg)

    switch class(msg)
        case 'cell'
            [msg dtype] = cellfun(@convert_datatype,msg,'UniformOutput',false);
            dtype = cellfun(@(x)x,dtype);
        case {'char','string'}
            dtype = 's';
        case {'double', 'single'}
            msg = single(msg);
            dtype = 'f';
        case 'logical'
            dtype = 'L';
        case {'int8','int16','int32','int64'}
            msg = int32(msg);
            dtype = 'i';
        otherwise 
            error('invalid osc data type ''%s''',cls);
    end
