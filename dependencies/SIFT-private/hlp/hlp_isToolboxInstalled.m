function res = hlp_isToolboxInstalled(toolboxName)
% check if a specific toolbox is installed (ignoring case)
% returns logical true if the toolbox is installed, else false

persistent allToolboxNames;

if isempty(allToolboxNames)
    v = ver;
    allToolboxNames = {v.Name};
end

res = any(strcmpi(toolboxName, allToolboxNames));