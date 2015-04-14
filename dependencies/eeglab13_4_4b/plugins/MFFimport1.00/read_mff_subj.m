function subj = read_mff_subj(filePath)

if ~exist(filePath,'dir')
    error('Unable to open %s.',filePath);
end

subj = struct();

subjObj = mff_getObject(com.egi.services.mff.api.MFFResourceType.kMFF_RT_Subject, 'subject.xml', filePath);
fieldVec = subjObj.getFields;

fprintf('Subject information in %s:\n',filePath);

for f = 0:fieldVec.size-1
    subjField = fieldVec.elementAt(f);
    fprintf('%s (%s): %s\n',char(subjField.getName),char(subjField.getDataType),char(subjField.getData));
    subj.(char(subjField.getName)) = char(subjField.getData);
end