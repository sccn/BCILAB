function helpText = hlp_readHelpInDeployedMode(functionName)
try
    filename = utl_whichfile(functionName);
catch
    error('SIFT:hlp_getFcnPreambleText:badFcnName','Unable to find function %s.m',functionName);
end
helpText = '';
fid = fopen(filename,'r');
finishup = onCleanup(@() fclose(fid));
if ~feof(fid), line = fgetl(fid);end
if ~isempty(line) && line(1) == '%', helpText = sprintf('%s\n%s',helpText,line(2:end)); end
if ~feof(fid), line = fgetl(fid);end
if ~isempty(line) && line(1) == '%', helpText = sprintf('%s\n%s',helpText,line(2:end)); end
if ~feof(fid), line = fgetl(fid);end
if ~isempty(line) && line(1) == '%', helpText = sprintf('%s\n%s',helpText,line(2:end)); end
while ~feof(fid)
    line = fgetl(fid);
    if isempty(line), break;end
    if line(1) ~= '%'
        break
    else
        helpText = sprintf('%s\n%s',helpText,line(2:end));
    end
end