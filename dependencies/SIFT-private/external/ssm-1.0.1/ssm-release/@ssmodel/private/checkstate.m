function checkstate(model, alpha, name)
if ~isnumeric(alpha), error('ssm:ssmodel_checkstate:InputError', [name ' must be numeric.']); end
if size(model.T, 1) ~= size(alpha, 1), error('ssm:ssmodel_checkstate:IncompatibleModel', 'state dimension of model does not match state vector sequence.'); end
if model.n ~= 1 && size(alpha, 2) > model.n, error('ssm:ssmodel_checkstate:IncompatibleModel', 'model time duration is less than state vector sequence.'); end
