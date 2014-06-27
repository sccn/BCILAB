function v = env_version
% Get the current version of BCILAB
v = '1.2-alpha';
if isdeployed
    v = [v ' compiled']; end
