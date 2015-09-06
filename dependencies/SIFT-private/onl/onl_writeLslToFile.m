function buffername = onl_writeLslToFile(varargin)
% Open a named LSL stream and write it to disk.
% If no arguments are supplied, a GUI is generated.
% This function requires the BCILAB toolbox (http://sccn.ucsd.edu/wiki/BCILAB)
%
% --------------------------------------------------------------------------------------------------
% Input                  Information                                                                    
% --------------------------------------------------------------------------------------------------
% StreamName:            LSL stream that should be written to disk                                      
%                        This is the name of the stream.                                                
%                        Possible values: ''                                                            
%                        Default value  : 'n/a'                                                         
%                        Input Data Type: string                                                        
%                                                                                                       
% OutputOptions:         Output options                                                                 
%                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
%     | FileName:        The file to write to. This may also contain a system or platform-independent path              
%                        Input Range  : Unrestricted                                                    
%                        Default value: 'data:/lastdata.set'                                                    
%                        Input Data Type: string                                                        
%                                                                                                       
%     | UpdateFrequency: Update frequency                                                               
%                        This is the rate at which data is written.                                     
%                        Input Range  : Unrestricted                                                    
%                        Default value: 1                                                               
%                        Input Data Type: real number (double)                                          
%                                                                                                       
%     | StartDelay:      Start-up delay                                                                 
%                        Delay before real-time processing begins; grace period until file is written.  
%                        Input Range  : Unrestricted                                                    
%                        Default value: 3                                                               
%                        Input Data Type: real number (double)  
%
% 
% 
% See Also: run_writedataset(), run_readlsl()
%
% Tim Mullen, SCCN/INC/UCSD 2013
% Inspired by BCILAB functions by Christian Kothe

% make sure that everything is on the path and LSL is loaded
if ~exist('arg_define','file')
    addpath(genpath(fileparts(mfilename('fullpath')))); end
% make sure that the library is also found if the toolbox is compiled
lib = lsl_loadlib(env_translatepath('dependencies:/liblsl-Matlab/bin'));


% handle default arguments
try
streamnames = find_streams;
catch err
    hlp_handleerror(err);
    streamnames = {''};
end
% handle inputs
opts = arg_define(varargin, ...
        arg({'streamname','StreamName'},streamnames{1},streamnames,'LSL stream that should be written to disk. This is the name of the stream.'), ...
        arg_sub({'outopts','OutputOptions'},{'out_filename','data:/lastdata.set'},@run_writedataset,'Output options','suppress',{'in_stream'}) ...
        );

if isempty(opts.streamname)
    return;
end
if isempty(varargin)
    % bring up GUI dialog if no arguments were passed
    arg_guidialog;
else
    % choose variable names to hold the stream's data (in the base workspace)
    taken = evalin('base','whos(''lsl*'')');
    buffername = genvarname(['lsl_' opts.streamname '_stream'],{taken.name});
    
    % create local stream to read data from LSL
    run_readlsl('MatlabStream',buffername,'SelectionProperty','name','SelectionValue',opts.streamname);
    
    % write data to file
    fprintf('Will write stream ''%s'' to %s\n',opts.streamname,env_translatepath(opts.outopts.out_filename));
    fprintf('To end recording, type:\n>> clear %s \n',buffername);
    run_writedataset(opts.outopts,'in_stream',buffername)
end

% find names of streams on the lab network...
function names = find_streams
    streams = lsl_resolve_all(lib,1);
    names = cellfun(@(s)s.name(),streams ,'UniformOutput',false);
    if isempty(names)
        error('There is no stream visible on the network.'); end

end

end
