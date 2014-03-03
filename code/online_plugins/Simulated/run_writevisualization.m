function run_writevisualization(varargin)
% Produce a real-time visualization within MATLAB.
% run_writevisualization(Model,SourceStream,VisFunction,UpdateFrequency,OutputForm,CreateFigure,StartDelay,PredictorName)
%
% In:
%   Model : A model data structure (as obtained from bci_train) based on which the predictor shall be 
%           created; typically this is a model struct, but for convenience it can be a file name, 
%           variable name in the base workspace, or a cell array of {file name, variable name} to 
%           refer to a variable inside a .mat file. The model is not modified by this function.
%           (default: 'lastmodel')
%
%   SourceStreamNames : Optional names of stream data structures in the MATLAB base workspace to
%                       consider as possible data sources (previously created with onl_newstream); 
%                       if a stream contains all channels that are needed by the predictor, or 
%                       alternatively has the right number and type of channels it will be considered 
%                       as a potential source stream unless ambiguous. (default: 'laststream')
%
%   VisFunction : visualization function (default: 'bar(y);ylim([0 1]);')
%                 y is the current signal value, f is the figure handle
%                 can also be given as a unary or binary function handle
%
%   UpdateFrequency : The rate at which new outputs will be computed. (default: 10)
%
%   OutputForm : output data form, see onl_predict (default: 'distribution')
%
%   CreateFigure : Whether a figure should be created before the visualization function is invoked.
%                  (default: true)
%
%   StartDelay : Delay before real-time processing begins; grace period until figure is created 
%                (default: 1)
%
%   PredictorName : Name of the predictor to be created; a variable of this name will be created in 
%                   the MATLAB base workspace to hold the predictor's state. If a variable with this
%                   name already exists it will be overridden. (default: 'lastpredictor')
%
% Examples:
%   % open a new processing stream, reading from an input stream called 'mystream', and using a 
%   % predictive model called 'mymodel', and display the outputs in a simple bar plot
%   run_writevisualization('mymodel','mystream')
%
%   % as before, but use a custom visualization function
%   run_writevisualization('mymodel','mystream',@myvisualizer)
%
%   % as before, but pass the arguments by name
%   run_writevisualization('Model','mymodel','SourceStream','mystream','VisFunction',@myvisualizer)
%
%   % as before, but update at 25Hz (instead of 10Hz)
%   run_writevisualization('Model','mymodel','SourceStream','mystream','VisFunction',@myvisualizer,'UpdateFrequency',25)
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-11-19

% declare the name of this component (shown in the menu)
declare_properties('name','MATLAB visualization');

% define arguments
opts = arg_define(varargin, ...
    arg({'pred_model','Model'}, 'lastmodel', [], 'Predictive model. As obtained via bci_train or the Model Calibration dialog.','type','expression'), ...
    arg({'in_stream','SourceStreamNames','SourceStream'}, 'laststream',[],'Input Matlab stream name(s). Optional names of stream data structures in the MATLAB base workspace to consider as possible data sources (previously created with onl_newstream); if a stream contains all channels that are needed by the predictor, or alternatively has the right number and type of channels it will be considered as a potential source stream unless ambiguous.'), ...
    arg({'vis_func','VisFunction'},'bar(y);ylim([0 1]);',[],'Visualization function. Function of y (the current prediction) and possibly f (the target figure); can be an expression or a function handle.'), ...
    arg({'update_freq','UpdateFrequency'},10,[0 Inf],'Update frequency. The rate at which new outputs will be computed.'), ...
    arg({'out_form','OutputForm'},'distribution',{'expectation','distribution','mode'},'Form of the produced output values. Can be the expected value (posterior mean) of the target variable, or the distribution over possible target values (probabilities for each outcome, or parametric distribution), or the mode (most likely value) of the target variable.'), ...
    arg({'create_fig','CreateFigure'}, true, [],'Create a figure. If the VisFunction assumes that a figure is present, this should be set to true.'), ...
    arg({'start_delay','StartDelay'}, 1, [],'Start-up delay. Delay before real-time processing begins; grace period until figure is created.'), ...
    arg({'predict_at','PredictAt'}, {},[],'Predict at markers. If nonempty, this is a cell array of online target markers relative to which predictions shall be made. If empty, predictions are always made on the most recently added sample.','type','expression'), ...
    arg({'pred_name','PredictorName'}, 'lastpredictor',[],'Name of new predictor. A variable of this name will be created in the MATLAB base workspace to hold the predictor''s state. If a variable with this name already exists it will be overridden.'), ...
    arg({'verbose','Verbose'}, false,[],'Verbose output. If false, the console output of the online pipeline will be suppressed.'));

% visualization function given as function handle?
if opts.vis_func(1) == '@'
    opts.vis_func = eval(opts.vis_func); end

% optionally create figure
if opts.create_fig
    fig = figure;
    newplot(fig);
    drawnow;
end

% start background writer job
onl_write_background( ...
    'ResultWriter',@(y)visualize(y,opts.vis_func,fig),...
    'MatlabStream',opts.in_stream, ...
    'Model',opts.pred_model, ...
    'OutputFormat',opts.out_form, ...
    'UpdateFrequency',opts.update_freq, ...
    'PredictorName',opts.pred_name, ...
    'PredictAt',opts.predict_at, ...
    'Verbose',opts.verbose, ...
    'StartDelay',opts.start_delay,...
    'EmptyResultValue',[]);

% background visualization function
function visualize(yy,visfunc,fig)
% for each prediction...
for k=1:size(yy,1)
    y = yy(k,:);
    % forward it to the visualization
    try
        if ischar(visfunc)
            % given as a string (expression)
            assignin('base','y',y);
            assignin('base','f',fig);
            evalin('base',visfunc);
        else
            % given as a function handle
            if nargin(visfunc) > 1
                visfunc(y,fig);
            else
                visfunc(y);
            end
        end
    catch e
        disp('Error in visualization function: ');
        hlp_handleerror(e);
    end
end
drawnow;
