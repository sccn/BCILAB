function [parameters,states] = bci_ConstructReal
% Handle a BCI2000 Construct() event: declare parameters for BCI2000 management.
% [Parameters,States] = bci_Construct()
%
% Out:
%   Parameters : cell array of strings, which are BCI2000 parameter declaration lines
%                (see http://www.bci2000.org/wiki/index.php/Technical_Reference:Parameter_Definition)
%
%   States : cell array of strings, which are BCI2000 state declaration lines
%            (see http://www.bci2000.org/wiki/index.php/Technical_Reference:State_Definition)
%
%                               Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                               2010-09-08

quicklog('bcilab-log.txt','============================================================');
quicklog('bcilab-log.txt','======= called, checking whether toolbox is loaded... ======');

% make sure that the toolbox is properly loaded
global tracking;
if isempty(tracking)
    error('BCILAB toolbox has not yet been loaded... please make sure that this is being called through the bci_Construct.m in the path/to/BCI2000/progs directory'); end

quicklog('bcilab-log.txt','determining default model to load...');

% determine the default model file to be loaded
models = dir(env_translatepath('resources:/models/*.mat'));
if isempty(models)
    modelfile = env_translatepath('resources:/models/lastmodel.mat');
else
    [dummy,order] = sort([models.datenum]);
    modelfile = env_translatepath(['resources:/models/' models(order(end)).name]);
end

quicklog('bcilab-log.txt','declaring parameters...');

% declare parameters
states = {};
parameters = {...
    ['BCILAB string Model= ' modelfile ' // The BCI model to be used for making predictions, previously generated via BCILAB, usually from some calibration data. (inputfile)'],...
    'BCILAB string MetaFile= (as-in-training-data) // An optional EEGLAB-supported file which specifies the meta-data of the stream, i.e. channel labels, sampling rate, etc. (inputfile)', ...
    'BCILAB float OutputRate= 10 10 0 1000 // Rate at which the predictor''s outputs will be sampled, in Hz', ...
    'BCILAB int OutputForm= 1 1 1 3 // Form/type of the BCI outputs. See help of utl_formatprediction.m: 1 expectation, 2 distribution, 3 mode (enumeration)', ...
    'BCILAB int OutputChannels= 0 0 0 1000 // Number of BCI output channels, if known. Note that this depends on the output form and the type of machine learning module that was used in BCILAB.', ...
    'BCILAB float MaxCPU= 1 1 0 1 // Maximum fraction of CPU time allocated for this predictor', ...
    };


quicklog('bcilab-log.txt','finished successfully...');
