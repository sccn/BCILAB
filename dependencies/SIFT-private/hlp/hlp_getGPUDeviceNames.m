function [GPUDeviceNames, GPUDeviceIdx, GPUDeviceCC] = hlp_getGPUDeviceNames(minComputeCapability)
% Return the names of all GPU devices on this machine
% This requires the Parallel Computing Toolbox
% Devices with Compute Capability less than minComputeCapability are
% ignored
% The Index of each device is returned in GPUDeviceIdx
%
% Author: Tim Mullen, SCCN/INC/UCSD 2014

if nargin<1
    minComputeCapability = 1.3; end
ws = warning('off','all');
dcount = gpuDeviceCount;
GPUDeviceNames = [];
GPUDeviceIdx   = [];
GPUDeviceCC    = [];

% save selected device
gs = gpuDevice();
for k=1:dcount
    % select the device
    if gs.Index==k
        g=gs;
    else
        g=gpuDevice(k); 
    end
    % skip devices with compute capability less than minComputeCapability
    GPUDeviceCC(k) = str2double(g.ComputeCapability);
    if GPUDeviceCC(k) < minComputeCapability
        continue; end
    % get name and index of the device
    GPUDeviceNames{end+1} = g.Name;
    GPUDeviceIdx(end+1)   = g.Index;
end
% restore gpu device
gpuDevice(gs.Index);
warning(ws);