function model = ssm_genair(nparam, nfreq, freq)

%SSM_GENAIR Create SSMODEL object for generalized airline model.
%   model = SSM_GENAIR(nparam, nfreq, freq)

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 3, error('ssm:ssm_genair:NotEnoughInputs', 'Insufficient input arguments.'); end
if ~isscalar(nparam) || ~isnumeric(nparam) || all(nparam ~= 3:4), error('ssm:ssm_genair:InputError', 'nparam must be 3 or 4.'); end
if ~isscalar(nfreq) || ~isnumeric(nfreq) || all(nfreq ~= 3:5), error('ssm:ssm_genair:InputError', 'nfreq must be 3 to 5.'); end
if ~isnumeric(freq) || (length(freq) ~= 6-nfreq), error('ssm:ssm_genair:InputError', ['freq must specify ' int2str(6-nfreq) ' frequencies.']); end
model   = ssm_freqspec(0, 1, 1, 0, 1, nparam, nfreq, freq);

