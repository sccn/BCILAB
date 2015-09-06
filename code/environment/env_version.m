function v = env_version
% Get the current version of BCILAB
v = '1.4-devel';
if isdeployed
    v = [v ' compiled']; end
