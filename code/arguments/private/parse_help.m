function help = parse_help(help,arg_name,summary_len)
% helper function for the arg* specifiers, to parse the help into a first and second part.
% Help = parse_help(Help,SummaryLength)
%
% In:
%   Help: some help specification (as it appears in the arg* functions
%
%   ArgumentName : name of the argument, for diagnostics
%   
%   SummaryLength : the maximum length for the executive summary portion
%
% Out:
%   Help: a cell array of {executive summary, description}
%
% See also:
%   arg_define
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-09-24

% Copyright (C) Christian Kothe, SCCN, 2010, christian@sccn.ucsd.edu
%
% This program is free software; you can redistribute it and/or modify it under the terms of the GNU
% General Public License as published by the Free Software Foundation; either version 2 of the
% License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
% even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License along with this program; if not,
% write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
% USA

if nargin < 3
    summary_len = 60; end

try
    initial_help = help;
    if ischar(help)
        % string: split at the end of the first sentence and put into a cell array
        [a,b] = regexp(help,'\.\s+[A-Z]','once');
        if ~isempty(a)
            help = {help(1:a-1), help(b-1:end)};
        else
            help = {help};
        end
    elseif ~iscellstr(help)
        error('The help text must be a string.');
    end

    % remove trailing dot
    if length(help{1}) > 1 && help{1}(end) == '.'
        help{1} = help{1}(1:end-1); end

    % check for length limit
    if length(help{1}) > summary_len
        % Note: The first sentence in the description is used in some GUIs which have a size limit;
        %       to prevent text from being cut off, please use a shorter wording in the first sentence.
        %
        %       Also note that if this sentence is not followed by a capital letter, the remaining part 
        %       is not considered separate.
        error('The executive summary for the given argument is too long.'); 
    end
catch e
    fprintf('Problem with the help text for argument %s: %s\n(text: %s)\n',arg_name,e.message,hlp_tostring(initial_help));
    help = {};
end
