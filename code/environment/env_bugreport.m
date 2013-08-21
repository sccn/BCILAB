function env_bugreport(varargin)
% File a bug report for BCILAB.
% env_bugreport(Summary,Description,Severity,ReproductingScript,DataAvailable)
% 
% This function uses MATLAB's built-in email capability, which relies on the SMTP server (and user
% credentials) being specified in MATLAB's preferences. See the documentation of sendmail() to see 
% how to set these up. If you are working in a research lab, these settings are probably already pre-
% configured for your account.
% 
% The data that is implicitly associated with the bug report is the MATLAB and toolbox versions, the 
% OS and compiler version, BCILAB version, as well as a snapshot of the MATLAB console output for 
% the current session.
%
% In:
%   Summary : executive summary of the bug (default: '')
%
%   Description : description of the bug (default: '')
%
%   Severity : severity level of the bug ('enhancement','minimal','normal','major','critical','blocker')
%
%   ReproducingScript : script file to reproduce the data (default: '')
%
%   DataAvailable : is data available upon request for fixing the bug? (default: false)
%
% Examples:
%   % file a very quick bug report (assuming that the console output is informative enough)
%   env_bugreport('tutorial_erd1 crashes on our Solaris mainframe');
%   
%   % file a report for a bug that just happened, and include with it a script that reproduces it
%   env_bugreport('Review/Edit GUI crashes when clicking Cancel 3x in a row', ...
%       'This happens every time when the attached script has previously been executed until line 9.', ...
%       'normal', ...
%       'bcilab:/userscripts/myscript.m',true);
%
%   % file a bug report for a set of problematic scripts script
%   env_bugreport('Summary','Scripts producing tons of annoying warnings','ReproducingScript',{'myscript1.m','myscript2.m'},'Severity','enhancement');
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2011-01-20

% define arguments
arg_define(varargin,...
    arg({'bug_summary','Summary'},'(summary)',[],'Executive summary. Should be shorter than 80 characters'), ...
    arg({'bug_description','Description'},'(description)',[],'Problem description. Extended description of the bug (excluding command line output).'), ...
    arg({'bug_severity','Severity'},'normal',{'enhancement','minimal','normal','major','critical','blocker'},'Problem severity level.'), ...
    arg({'bug_script','ReproducingScript'},'',[],'Script file that reproduces the error. A BCILAB script to reproduce the bug. If additional data / toolboxes required, please point out.','shape','row','type','char'), ...
    arg({'bug_data_avail','DataAvailable'},false,[],'Data is available. Whether data to reproduce the bug is available upon request.'));

% for now, done via sendmail
try
    % create attachments
    io_saveworkspace('bugreport.mat',true);    
    attachments = {'bugreport.mat'};
    if ~isempty(bug_script)
        if ~iscell(bug_script)
            bug_script = {bug_script}; end
        attachments = [attachments cellfun(@env_translatepath,bug_script,'UniformOutput',false)];
    end
    
    if strcmp('Yes', questdlg2('This will send an email with your command line output (for this session) and system specs to the BCILAB developer. Do you want to proceed?','Confirmation'))
        % send the email...
        sendmail('christiankothe@gmail.com',sprintf('[BCILAB %s] bug: %s / %s',env_version,bug_severity,bug_summary),[bug_summary 10 13 bug_description ...
            10 13 'Data available: ' hlp_rewrite(bug_data_avail,false,'no',true,'yes')],attachments);
    else
        disp('You can review the bug report contents in the file bugreport.mat and submit it manually.');
    end
catch e
    if strcmp(e.message, 'MATLAB:sendmail:SMTPServerIndeterminate')
        disp('You need to configure your email preferences in matlab; see documentation of sendmail.');
        if ~isdeployed
            doc sendmail; end
    else
        disp(['Error while trying to send email: ' e.message]);
    end
end
