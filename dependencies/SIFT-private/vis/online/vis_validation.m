function vis_validation(whitestats,PCstats,stability,whitenessCriteria)
% display results of model validation

% fprint validation results
fprintf('\n\n------------------------------\n');
fprintf('Model Validation:\n');

if ~isempty(whitestats)
    variableNames     = hlp_variableize(lower(whitenessCriteria));
    if ~iscell(variableNames)
        variableNames = {variableNames};
    end
    fprintf('\tResidual Whiteness:\n');
    for i=1:length(whitenessCriteria)
        fprintf('\t\t%s:\tp=%0.10g\t(%s)\n', ...
            whitenessCriteria{i}, ...
            whitestats.(variableNames{i}).pval, ...
            fastif(whitestats.(variableNames{i}).w,'white','not white'));
    end
end
if ~isempty(stability)
    fprintf('\tModel Stability:\n')
    fprintf('\t\tln|Lambda|=%0.10g\t(%s)\n',max(stability.lambda),fastif(stability.stability,'stable','not stable'));
end
if ~isempty(PCstats)
    fprintf('\tModel Consistency:\n')
    fprintf('\t\tPC=%0.10g%%\n',PCstats.PC);
end
fprintf('------------------------------\n\n');