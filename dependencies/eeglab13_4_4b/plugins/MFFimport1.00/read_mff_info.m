function info = read_mff_info(filePath)

mffjarpath = which('MFF-1.0.d0004.jar');
if isempty(mffjarpath)
    error('MFF-1.0.d0004.jar not found in path! Make sure the mffreader directory is in your MATLAB path');
end

javaaddpath(mffjarpath);

if ~exist(filePath,'dir')
    error('Unable to open %s.',filePath);
end

info = struct();

infoObj = mff_getObject(com.egi.services.mff.api.MFFResourceType.kMFF_RT_Info, 'info.xml', filePath);

datetimestr = char(infoObj.getRecordTime);
[datestr,timestr] = strtok(datetimestr,'T');
timestr = strtok(timestr(2:end),'.');
info.date = datestr;
info.time = timestr;
info.ampserial = char(infoObj.getAmpSerialNumber);
info.ampfirmware = char(infoObj.getAmpFirmwareVersion);
info.moviedelta = char(infoObj.getMovieDelta);

% fprintf('Information about %s:\n',filePath);
% fprintf('Date and time: %s %s\n',datestr,timestr);
% fprintf('Amp serial number: %s\n',char(infoObj.getAmpSerialNumber));
% fprintf('Amp firmware version: %s\n',char(infoObj.getAmpFirmwareVersion));
% fprintf('Movie delta: %s\n',char(infoObj.getMovieDelta));

