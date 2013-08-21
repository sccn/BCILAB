function v = env_version
% Get the current version of BCILAB
v = '1.1-beta';
if isdeployed
    v = [v ' compiled']; end
