function v = env_version
% Get the current version of BCILAB
v = '1.3-testing';
if isdeployed
    v = [v ' compiled']; end
