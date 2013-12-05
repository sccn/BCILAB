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

% instantiate model class...
instance = eval(model.paradigm);

if isempty(varargin)
    arg_guidialog(@instance.visualize,'Parameters',{'Model',model});
else
    instance.visualize(model,varargin);
end
