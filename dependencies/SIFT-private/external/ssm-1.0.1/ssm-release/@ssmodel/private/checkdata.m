function checkdata(model, y, name)
if ~isnumeric(y), error('ssm:ssmodel_checkdata:InputError', [name ' must be numeric.']); end
if model.p ~= size(y, 1), error('ssm:ssmodel_checkdata:IncompatibleModel', 'data dimension of model does not match data.'); end
if model.n ~= 1 && size(y, 2) > model.n, error('ssm:ssmodel_checkdata:IncompatibleModel', 'model time duration is less than data.'); end
