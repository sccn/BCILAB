function bci_visualize(model,varargin)
% Visualize a model
% bci_visualize(Model)
%
% In:
%   Model : a predictive model as computed via bci_train
%
%   Options... : options for the respective paradigm's visualization
%                (if empty, a GUI will be opened)
%
% Examples:
%   % assuming that a model has been learned previously, visualize it
%   bci_visualize(lastmodel)
%
% See also:
%   bci_train
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-11-09

if ~isfield(model,'paradigm')
    error('The given data structure is not a valid model (lacking the required field .paradigm)'); end

% instantiate model class
try
    instance = eval(model.paradigm);
catch e
    error('Failed to instantiate the given BCI paradigm %s with error: %s.',char(model.paradigm),e.message);
end

if isempty(varargin) && arg_supported(@instance.visualize)
    arg_guidialog(@instance.visualize,'Parameters',{'Model',model});
else
    instance.visualize(model,varargin);
end
