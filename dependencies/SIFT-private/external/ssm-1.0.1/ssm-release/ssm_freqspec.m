function model = ssm_freqspec(p, d, q, P, D, nparam, nfreq, freq)

%SSM_FREQSPEC Create SSMODEL object for frequency specific SARIMA model.
%   model = SSM_FREQSPEC(p, d, q, P, D, nparam, nfreq, freq)

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 8, error('ssm:ssm_freqspec:NotEnoughInputs', 'Insufficient input arguments.'); end
if ~isscalar(p) || ~isnumeric(p), error('ssm:ssm_freqspec:InputError', 'p must be a scalar.'); end
if ~isscalar(d) || ~isnumeric(d), error('ssm:ssm_freqspec:InputError', 'd must be a scalar.'); end
if ~isscalar(q) || ~isnumeric(q), error('ssm:ssm_freqspec:InputError', 'q must be a scalar.'); end
if ~isscalar(P) || ~isnumeric(P), error('ssm:ssm_freqspec:InputError', 'P must be a scalar.'); end
if ~isscalar(D) || ~isnumeric(D), error('ssm:ssm_freqspec:InputError', 'D must be a scalar.'); end
if ~isscalar(nparam) || ~isnumeric(nparam) || all(nparam ~= 3:4), error('ssm:ssm_freqspec:InputError', 'nparam must be 3 or 4.'); end
if ~isscalar(nfreq) || ~isnumeric(nfreq) || all(nfreq ~= 3:5), error('ssm:ssm_freqspec:InputError', 'nfreq must be 3 to 5.'); end
if ~isnumeric(freq) || (length(freq) ~= 6-nfreq), error('ssm:ssm_freqspec:InputError', ['freq must specify ' int2str(6-nfreq) ' frequencies.']); end
[Z T R P1 Tmask Rmask P1mask]   = mat_sarima(p, d, q, P, D, 1, 12, false);
[fun gra psi]                   = fun_freqspec(p, q, P, nparam, nfreq, freq);
model                           = [ssm_null ssmodel(struct('type', 'freqspec', 'p', p, 'd', d, 'q', q, 'P', P, 'D', D, 'nparam', nparam, 'nfreq', nfreq, 'freq', freq), 0, Z, ssmat(T, Tmask), ssmat(R, Rmask), ssmat(0, true), [repmat('T', 1, p+P>0) 'RQP1'], fun, gra, psi, [], ssmat(P1, P1mask))];

