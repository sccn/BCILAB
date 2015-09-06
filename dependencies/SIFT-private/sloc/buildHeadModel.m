
function hmObj = buildHeadModel(eeg,surfFile,templateFile,plotModel,scalingFactor)
% surfFile is the path to the output surface file that will be generated
% templateFile is the path to an existing template file

% electrode locations [X Y Z]
elocs   = [[eeg.chanlocs.X]' [eeg.chanlocs.Y]' [eeg.chanlocs.Z]']*scalingFactor;
labels  = {eeg.chanlocs.labels};

% remove any electrodes which lack locations
unknownLocs = find(arrayfun(@(x) isempty(x.X), eeg.chanlocs));
if ~isempty(unknownLocs)
    lbl = cellfun(@(x) [x ', '],labels(unknownLocs),'UniformOutput',false);
    lbl{end} = lbl{end}(1:end-2);
    fprintf('Electrodes [%s] do not have locations. Removing...\n',char(lbl{:})');
    labels(unknownLocs) = [];
end

% 
% % build the head model object
% hmObj = headModel('surfaces',surfFile,'atlas',template.atlas,'fiducials',template.fiducials,'channelSpace',elocs,'label',labels);
% if plotModel
%     plotHeadModel(hmObj);
% end

hmObj = headModel('channelSpace',elocs,'label',labels);
hmObj.warpTemplate2channelSpace(templateFile,surfFile);
% hmObj.warpChannelSpace2Template(templateFile,surfFile);

if plotModel
    plotHeadModel(hmObj);
end

%% solving the forward problem with OpenMEEG
conductivity = [0.33 0.022 0.33]; % brain and scalp = 0.33 S/m, skull = 0.022 S/m; these conductivies were taken from
                                  % Valdes-Hernandez et al., 2009, Oostendrop TF, 2000; Wendel and Malmivuo, 2006      
normal2surface = true;
hmObj.computeLeadFieldBEM(conductivity,normal2surface);

