function bci_visualize(model,varargin)
% Visualize a model
% bci_visualize(Model)
%
% In:
%   Model : a predictive model as computed via bci_train
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
instance.visualize(model,varargin{:});