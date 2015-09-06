function cancel = hlp_confirmWaitbarCancel(waitbarTitle,message)
% Confirm cancellation of a waitbar with specified title
% If cancellation is confirmed, waitbar will be closed. otherwise, cancel
% status is reset
% Author: Tim Mullen, 2012, SCCN/INC, UCSD

if nargin < 2 || isempty(message)
    message = 'Are you sure you want to cancel?';
end
if strcmpi('yes',questdlg2( ...
        message, ...
        'Confirmation','Yes','No','No'));
    multiWaitbar(waitbarTitle,'Close');
    cancel = true;
else
    multiWaitbar(waitbarTitle,'ResetCancel',true);
    cancel = false;
end

