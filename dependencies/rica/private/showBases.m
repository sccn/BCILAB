function [] = showBases(theta, infoStruct, stuff, varargin)
global params
if mod(infoStruct.iteration,50) == 0
    W = reshape(theta, params.numFeatures, params.n);
    display_network(W');
end
