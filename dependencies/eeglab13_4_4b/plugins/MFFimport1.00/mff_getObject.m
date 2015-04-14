function theObject = mff_getObject(objType, filename, path)
URI = [path filesep filename];

oldSchool = true;
switch objType
    case com.egi.services.mff.api.MFFResourceType.kMFF_RT_Signal
        objStr = 'com.egi.services.mff.api.Signal';
    case com.egi.services.mff.api.MFFResourceType.kMFF_RT_EventTrack
        objStr = 'com.egi.services.mff.api.EventTrack';
    case com.egi.services.mff.api.MFFResourceType.kMFF_RT_Info
        objStr = 'com.egi.services.mff.api.Info';
    case com.egi.services.mff.api.MFFResourceType.kMFF_RT_SensorLayout
        objStr = 'com.egi.services.mff.api.SensorLayout';
    case com.egi.services.mff.api.MFFResourceType.kMFF_RT_Categories
        objStr = 'com.egi.services.mff.api.Categories';
    case com.egi.services.mff.api.MFFResourceType.kMFF_RT_Subject
        objStr = 'com.egi.services.mff.api.Subject';
    otherwise
        oldSchool = false;
end
if oldSchool
    switch objStr
        case {'com.egi.services.mff.api.Info',...
         'com.egi.services.mff.api.EventTrack',...
         'com.egi.services.mff.api.SensorLayout',...
         'com.egi.services.mff.api.Categories',...
         'com.egi.services.mff.api.Subject'}
            theObject = javaObject(objStr, true);
            theObject = theObject.unmarshal(URI, true);
        case 'com.egi.services.mff.api.Signal'
            theObject = javaObject(objStr);
            theObject = theObject.unmarshal(URI);
%         case 'com.egi.services.mff.api.Epochs'
%             theObject = javaObject(objStr);
%             theObject = theObject.unmarshal(URI,true);
    end
else
    delegate = javaObject('com.egi.services.mff.api.LocalMFFFactoryDelegate');
    factory = javaObject('com.egi.services.mff.api.MFFFactory', delegate);
    resourceVal = objType;
    resourceType = javaObject('com.egi.services.mff.api.MFFResourceType', resourceVal);
%    fprintf('%s %s\n', char(URI), char(resourceType));
    theObject = factory.openResourceAtURI(URI, resourceType);
    theObject.loadResource();
end

% com.egi.services.mff.api.MFFResourceType.kMFF_RT_Info = 7
% com.egi.services.mff.api.MFFResourceType.kMFF_RT_EventTrack = 3
% com.egi.services.mff.api.MFFResourceType.kMFF_RT_Categories = 9
% com.egi.services.mff.api.MFFResourceType.kMFF_RT_Signal = 2
% com.egi.services.mff.api.MFFResourceType.kMFF_RT_Epochs = 4
% com.egi.services.mff.api.MFFResourceType.kMFF_RT_InfoN = 8
% com.egi.services.mff.api.MFFResourceType.kMFF_RT_SensorLayout = 10
% com.egi.services.mff.api.MFFResourceType.kMFF_RT_Coordinates = 11
% com.egi.services.mff.api.MFFResourceType.kMFF_RT_History = 6
% com.egi.services.mff.api.MFFResourceType.kMFF_RT_MFFFile = 1
% com.egi.services.mff.api.MFFResourceType.kMFF_RT_Photogrammetry = 12
% com.egi.services.mff.api.MFFResourceType.kMFF_RT_Subject = 5
% com.egi.services.mff.api.MFFResourceType.kMFF_RT_Any = 0
% com.egi.services.mff.api.MFFResourceType.kMFF_RT_Unknown = -1
