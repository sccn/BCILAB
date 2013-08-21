function res = run_pipedatariver(varargin)
% Run BCILAB's real-time engine as a module in the DataRiver real-time system.
% run_pipedatariver(Model,InputMetadata,InputStream,OutputStream,OutputForm,OutputSamplingRate,MaximumCPULoad)
%
% BCILAB reads raw data from some input stream and writes predictions to some output stream.
%
% In:
%	Model : a file or struct that contains a predictive model as previously computed by bci_train
%           (default: 'lastmodel')
%
%   InputMetadata : Meta-data of the input stream. This is a struct or cell array of name-value pairs 
%                   with meta-data fields to use. The mandatory fields are 'srate' and 'chanlocs',
%                   where chanlocs is either a channel locations struct array, or a cell array of
%                   channel names, or the number of channels (in which case a cell array of the form
%                   {'A1','A2', ..., 'A32','B1', ...} is created). Optionally, the field
%                   'datasource' can be set to point to a dataset on disk or in a MATLAB workspace
%                   variable.
%
%	InputStream  : name of the input stream to be used (e.g., from the acquisition system) (default:
%                  'C:/tmp/DataRiver')
%
%   OutputStream : name of the output stream (carrying the control signal); must be a valid file
%                  name (default: 'C:/tmp/BCILAB')
%
%   OutputForm : form of the output stream, can be one of the following:
%                'expectation' : the expected value (= posterior mean) of the outputs; can be multi-
%                                dimensional but is usually 1d (default) this mode is appropriate
%                                for simple applications that expect a smooth control signal, or for
%                                applications that expect a regression output
%                'distribution' : parameters of the output distribution; for discrete distributions,
%                                 this is one probability value for each target (adding up to 1)
%                                 this mode is appropriate for more advanced applications that use
%                                 the full output distribution (e.g., for decision-theoretical
%                                 processing)
%                'mode' : the most likely output value (currently only supported for discrete 
%                         distributions) this mode is appropriate for simple applications that take
%                         a non-probabilistic classifier decision (e.g., as from a Support Vector
%                         Machine)
% 
%   OutputSamplingRate :  rate, in Hz, at which the output stream shall be sampled (default: 10)
%                          note: if the sampling rate is higher than the rate at which BCILAB can
%                                process data, samples will be replicated (i.e., the same output
%                                comes multiple times)
%
%   MaximumCPULoad : the maximum fraction of CPU time allocated for this predictor (default: 0.6)
%                    note: this should leave some room for other important processes (e.g.,
%                          recording, data transfer, visualization, or other BCILAB predictors), but
%                          setting it too low (to sustain the given sampling rate) will lead to
%                          output samples being replicated
%
% Example:
%   % run BCILAB as a processing node in a datariver enironment; use the model stored in mypredictor.mat, 
%   % read from C:/tmp/DataRiver, with some given meta-data, and write to C:/tmp/BCILAB (default output stream)
%   run_pipedatariver('Model','mypredictor.mat', 'InputStream','C:/tmp/DataRiver', ...
%       'InputMetadata',{'srate',512,'channels',{'C3','C4','Cz','FP1', ...}});
%
% Notes:
%   The DataRiver real-time system is in active development, and therefore, only a beta version can
%   be obtained at this point.
%
% References:
%   [1] Vankov A., Bigdely-Shamlo N.,  Makeig S. "DataRiver ? A software platform for real time management of multiple data streams"
%       Fourth International BCI Meeting, Carmel, CA, June 2010  
%   [2] http://sccn.ucsd.edu/wiki/DataSuite
% 
%                               Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                               2010-09-08

% declare the name of this component (shown in the menu)
declare_properties('name','DataRiver');

default_prefix = fastif(ispc,'C:/tmp/','/tmp/');

% read options
opts = arg_define(varargin, ...
    arg({'model','Model','predictor','detector','classifier'}, 'lastmodel', [], 'Predictive model. As obtained via bci_train or the Model Calibration dialog.','type','expression'), ...
    arg_sub({'in_metadata','InputMetadata'},{},@utl_parse_metadata, 'Meta-data of the input stream. These are fields as they appear in EEGLAB data sets; only sampling rate and channel labels are mandatory.'), ...
    arg({'in_stream','InputStream'}, [default_prefix 'DataRiver'],[],'Input DataRiver stream. This is the stream that shall be analyzed and processed.'), ...
    arg({'out_stream','OutputStream'}, [default_prefix 'BCILAB'],[],'Output DataRiver stream. This is the stream that is produced by BCILAB.'), ...
    arg({'out_form','OutputForm'}, 'distribution',{'expectation','distribution','mode'},'Format of the produced output values. Can be the expected value (posterior mean) of the target variable, or the distribution over possible target values (probabilities for each outcome, or parametric distribution), or the mode (most likely value) of the target variable.'), ...
    arg({'out_srate','OutputSamplingRate'}, 10,[],'Output sampling rate. This is the rate at which estimate should be computed. If the sampling rate is higher than the CPU can handle, outputs will be repeated.'), ...
    arg({'max_cpu','MaximumCPULoad'}, 0.6,[],'Maximum CPU load. BCILAB will attempt to not exceed this value, which leaves room for other processes (e.g., recording, data transfer or other preditors).'));

% parse the meta-data spec
meta = utl_parse_metadata(opts.in_metadata);
   
% load DataSuite, if necessary
global ds_lib;
if isempty(ds_lib)
    disp('Loading DataSuite...');
    startup_ds; 
end

% set up a BCILAB online stream
onl_newstream('stream_datariver',meta);

% load the predictor
onl_newpredictor('predictor_datariver',opts.model,'stream_datariver');

% open a DataSuite input stream
in_stream = ds_OpenRead(opts.in_stream);
in_cleanup = onCleanup(@()ds_CloseRead(in_stream));

% open a DataSuite output stream
out_stream = ds_OpenWrite(opts.out_stream,opts.out_srate);
out_cleanup = onCleanup(@()ds_CloseWrite(out_stream));

% our temporary sample data structure
tmpsample = ds_array;
tmpsample.nItems = 1;
tmpsample.Data = 0;

% run
in_samples = 0;  % number of samples received so far
out_samples = 0; % number of samples produced so far
while 1
    t0 = tic;
    % while we have samples to read until the next output sample can be produced...
    while round(in_samples*(opts.out_srate/opts.in_stream.srate)) <= out_samples
        % read in at least one more sample
        received = false;
        while ~received
            % try to get a new sample
            [received,sample] = ds_Read(in_stream);
            if ~received
                pause(1/opts.in_metadata.srate); end
        end        
        % feed sample to BCI
        onl_append(double(sample.Data(2:sample.nItems)'));
        in_samples = in_samples+1;
    end
    
    % make a prediction to produce an output sample
    t1 = tic;
    prediction = onl_predict('predictor_datariver',opts.out_form);
    if ~(isnumeric(prediction) && isnan(prediction))
        % valid data: set up the sample
        tmpsample.nItems = length(prediction);
        tmpsample.Data = prediction;
    else
        % keep the dimensionality of the last output (or dummy output), but set data to NaN's
        tmpsample.Data = tmpsample.Data.*NaN;
    end
    % write out the data
    ds_Write(out_stream,tmpsample);
    out_samples = out_samples+1;
    
    % make sure that we are not exceeding our CPU budget
    computation_time = toc(t1);
    if opts.max_cpu < 1
        pause(computation_time * (1 - opts.max_cpu) / opts.max_cpu); end
    
    % check how much processing time is left for this output sample
    time_spent = toc(t0);
    time_left = 1/opts.out_srate - time_spent;
    if time_left < 0
        % we exceeded our time budget: replicate samples
        for k=1:ceil(-time_left*opts.out_srate)
            ds_Write(out_stream,tmpsample);
            out_samples = out_samples+1;
        end
    end
end
