function out_signal_dim = bci_Preflight(in_signal_dim)
% Check bci_Parameters variable (stored by BCI2000) for consistency.
% Out-Signal-Dim = bci_Preflight(In-Signal-Dim)
%
% In:
%  In-Signal-Dim  : a two-element vector of the form [#Channels, #Samples] which describes the input signal's dimensions
%
% Out:
%  Out-Signal-Dim : a two-element vector of the form [#Channels, #Samples] which describes the input signal's dimensions
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-09-08

global bci_Parameters;

% these variables are shared with Initialize and Process
global model meta_data outrate maxcpu outdims outform;

quicklog('bcilab-log.txt','called with in_signal_dim = [%i, %i]',in_signal_dim(1),in_signal_dim(2));

if in_signal_dim(1) <= 2
    quicklog('bcilab-log.txt','Haha, apparently bci_Preflight is being called with its own output by BCI2000... ignoring that for now.');
    out_signal_dim = in_signal_dim;
    return; 
end

% --- model checks ---

quicklog('bcilab-log.txt','checking if model file (%s) is present ...',bci_Parameters.Model{1});

% check if the bci model is present
if ~exist(bci_Parameters.Model{1},'file')
    error('BCILAB:BCI2000_plugin:model_not_found',['Could not find the bci model file specified as ' bci_Parameters.Model '.']); end

% --- stream meta-data checks ---

quicklog('bcilab-log.txt','checking if stream parameters are present...');
quicklog('bcilab-log.txt','declared parameters: %s',hlp_tostring(fieldnames(bci_Parameters)'));

% try to get the sampling rate
try
    srate = sscanf(bci_Parameters.SamplingRate{1},'%i');
    if ~isscalar(srate)
        quicklog('bcilab-log.txt','SamplingRate is ill-formatted: %s...',hlp_tostring(srate));
        error('Sampling rate has the wrong size...');
    end
    quicklog('bcilab-log.txt','SamplingRate is %.1f...',srate);
catch
    quicklog('bcilab-log.txt','SamplingRate is unspecified...');
    srate = [];
end

% try to get the channel names
try
    channels = bci_Parameters.ChannelNames;
    quicklog('bcilab-log.txt','ChannelNames are %s...',hlp_tostring(channels));
catch
    quicklog('bcilab-log.txt','ChannelNames are unspecified...');
    channels = [];
end

quicklog('bcilab-log.txt','checking if the meta-data file (%s) is present...',bci_Parameters.MetaFile{1});

% check if the meta-data file exists
if exist(bci_Parameters.MetaFile{1},'file')
    quicklog('bcilab-log.txt','found meata-data file, now loading...');
    % yes: use it
    meta_data = utl_parse_metadata('datasource',bci_Parameters.MetaFile{1},'srate',srate,'chanlocs',channels);
    quicklog('bcilab-log.txt','success.');
else
    % no: check if the BCI2000 parameters are sufficient to determine the stream format
    quicklog('bcilab-log.txt','Found no meata-data file.');
    if ~isempty(channels) && ~isempty(srate)
        quicklog('bcilab-log.txt','But ChannelNames and SamplingRate are fine, proceeding...');
        % we're fine...
        meta_data = utl_parse_metadata('srate',srate,'chanlocs',channels);
    else
        % we're NOT fine -- take a look at the model to find out the original (source) data set
        quicklog('bcilab-log.txt','ChannelNames or SamplingRate parameters are not specified correctly!');
        quicklog('bcilab-log.txt','Now checking if the original training set was specified in the model and can be retrieved for reference.');

        % try to load the model (jumping through a few hoops, if necessary....)
        try
            quicklog('bcilab-log.txt','Trying to load the model file...');
            model = io_load(bci_Parameters.Model{1});
            quicklog('bcilab-log.txt','success.');
            if ~isfield(model,'tracking') || ~isfield(model.tracking,'prediction_function')
                quicklog('bcilab-log.txt','But the model file is lacking the correct properties - apparently the file contains a variable that is the model.');
                % the loaded model is lacking the appropriate fields; check if there are variables
                % in the loaded data which are valid models
                candidates = {};
                for f = fieldnames(model)'
                    fname = f{1};
                    if isfield(model.(fname),'tracking') && isfield(model.(fname).tracking,'prediction_function')
                        candidates{end+1} = fname; end %#ok<AGROW>
                end
                if length(candidates) > 1
                    quicklog('bcilab-log.txt','Model file contains multiple candidate model variables; choice is ambiguous - aborting...');
                    error('BCILAB:bci_Preflight:ambiguous',['The file given as the model contains multiple candiate variables:\n' ...
                        hlp_tostring(candidates) '; please pass a file name which contains only one model variable.']);
                elseif isempty(candidates)
                    quicklog('bcilab-log.txt','Model file contains no candidate model variable - aborting...');
                    error('BCILAB:bci_Preflight:load_error','The given file contains no valid model.');
                else
                    quicklog('bcilab-log.txt','Model was successfully identified.');
                    model = model.(candidates{1});
                end
            end
        catch le
            quicklog('bcilab-log.txt','Error loading the model; traceback: %s',le.message);
            error('BCILAB:bci_Preflight:load_error',['The given model string could not be parsed; traceback: ' le.message]);
        end

        % now try to load the source data...
        if ~isfield(model,'source_data')
            quicklog('bcilab-log.txt','Model contains no source_data field, thus the source data set cannot be identified - aborting.');            
            error('BCILAB:bci_Preflight:load_error','The given model is lacking a source_data field (which should contain the name of the training set.'); 
        end
        if ~all(isfield(model.source_data,{'head','parts'}))
            quicklog('bcilab-log.txt','Model has a source_data field, but it does not refer to a single data set - aborting.');            
            error('BCILAB:bci_Preflight:load_error','The source_data field of the given model does not refer to a single data set (but some structure that is unrecognized by this plugin). You will have to specify the stream meta-data either via the MetaFile field or via SamplingRate and ChannelNames parameters.'); 
        end

        quicklog('bcilab-log.txt','Successfully extracted the meta-data from the source file.');            
        meta_data = exp_eval(model.source_data);
        
        if isfield(meta_data,'tracking')
            meta_data = rmfield(meta_data,'tracking'); 
            quicklog('bcilab-log.txt','Stripped ''tracking'' field from meta-data.');
        end
        if isfield(meta_data,'data')
            meta_data = rmfield(meta_data,'data'); 
            quicklog('bcilab-log.txt','Stripped ''data'' field from meta-data.');            
        end
        
        if ~isempty(srate)
            quicklog('bcilab-log.txt','Overriding the sampling rate with BCI2000''s %.1f',srate);
            meta_data.srate = srate; 
        end
    end
end

quicklog('bcilab-log.txt','Determined sampling rate is: %.1f',meta_data.srate);            
quicklog('bcilab-log.txt','Determined channel count: %i',length(meta_data.chanlocs));            
quicklog('bcilab-log.txt','Determined channel names are: %s',hlp_tostring({meta_data.chanlocs.labels}));            

% check that the input signal dimension matches the number of channel names
if length(meta_data.chanlocs) ~= in_signal_dim(1)
    quicklog('bcilab-log.txt','Number of channel names (%i) does not match the input signal channel # supplied by BCI2000 (%i) - aborting...',length(meta_data.chanlocs),in_signal_dim(1));
    error('BCILAB:BCI2000_plugin:metadata_inconsistent',['The number of channel names (' num2str(length(meta_data.chanlocs)) ', specified in the ChannelNames parameter or the MetaFile) does not match the signal''s channel count (' num2str(in_signal_dim(1)) ').']); 
end

quicklog('bcilab-log.txt','successfully matched specified and supplied channel numbers.');

% --- get remaining parameters ---
outforms = {'expectation','distribution','mode'};
outform = outforms{str2num(bci_Parameters.OutputForm{1})};
outrate = str2num(bci_Parameters.OutputRate{1});
maxcpu = str2num(bci_Parameters.MaxCPU{1});

quicklog('bcilab-log.txt','OutputForm is %s',outform);
quicklog('bcilab-log.txt','OutputRate is %.1f',outrate);
quicklog('bcilab-log.txt','MaxCPU is %.1f',maxcpu);

outdims = str2num(bci_Parameters.OutputChannels{1});
if outdims == 0
    quicklog('bcilab-log.txt','No number of output channels was specified; assuming univariate regression or classification...');
    if strcmp(outform,'distribution')
        quicklog('bcilab-log.txt','Since the OutputForm is ''distribution'', assuming that the output signal has 2 channels.');
        outdims = 2;
    else
        quicklog('bcilab-log.txt','Since the OutputForm is either ''mode'' or ''expectation'', assuming that the output signal has 1 channel.');
        outdims = 1;
    end
    % no output channel # specified; assume that we're doing univariate regression or classification
    disp(['Number of BCI output channels was not specified; assuming it is ' num2str(outdims)]);
end

out_signal_dim = [outdims, in_signal_dim(2)];

quicklog('bcilab-log.txt','Returing out_signal_dim = [%i, %i]',out_signal_dim(1),out_signal_dim(2));

if ~isempty(meta_data)
    quicklog('bcilab-log.txt','global meta_data variable is properly initialized.'); end
    
quicklog('bcilab-log.txt','finished successfully...');
