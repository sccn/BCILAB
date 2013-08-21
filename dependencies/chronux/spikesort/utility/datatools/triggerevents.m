function [events,trigger] = triggerevents(trace, trigger, pre, post, gap)
%TRIGGEREVENTS     Extracts data samples before/after a trigger signal.
%   E = TRIGGEREVENTS(X, TRIGGER, PRE, POST), where X is an M x 1 vector
%   and TRIGGER is a length M vector with P 0->~0 transitions, returns a
%   P x (PRE+1+POST) matrix E.  The j-th row of E contains values of X
%   before and after the index of the j-th 0->~0 transition in TRIGGER;
%   PRE and POST specify the number pre- and post-transition samples,
%   respectively.  While X and TRIGGER can be of any numeric type, E will
%   always be of type double.
%
%   E = TRIGGEREVENTS(X, TRIGGER, PRE, POST, GAP) ignores 0->~0
%   transitions in TRIGGER which occur <= GAP samples after a previous  
%   unignored transition.
%
%   [E,TRIGLIST] = TRIGGEREVENTS(X, TRIGGER, ...) also returns indices
%   of TRIGGER at which 0->~0 transitions occur.
%
%   E = TRIGGEREVENTS(X, TRIGLIST, PRE, POST) directly specifies trigger
%   crossing indices rather than inferring them from 0->~0 transitions.
%   Here, TRIGLIST is a length P vector (P ~= M) containing only values
%   between 1 and M that are interpreted as row indices into X.  The j-th
%   row of the resulting E matrix contains values of X indexed relative to
%   the j-th index in TRIGLIST.  The GAP syntax described above is not
%   allowed in this case.
%
%   In all of the above cases, if X is an M x N matrix, E will be of size
%   P x (PRE+1+POST) x N and E(:,:,k) will be equal to the result of
%   TRIGGEREVENTS(X(:,k), ...).
%
%   When a trigger is fewer than PRE samples after the start of the X or
%   fewer than POST samples before the end, NaN values are returned for
%   the invalid X samples.
%
%   Examples:
%       triggerevents([1 2 3 4 5 6 7 8 9]', [1 1 0 0 1 1 1 1 1]', 2, 3)
%   and triggerevents([1 2 3 4 5 6 7 8 9]', [1 5], 2, 3)
%              both return   [NaN NaN 1 2 3 4;  3 4 5 6 7 8].
%
%       triggerevents([1 2 3 4 5 6 7 8 9]', [1 1 0 0 1 0 1 1 1]', 2, 3, 4)
%                   returns  [NaN NaN 1 2 3 4;  5 6 7 8 9 NaN].
%
%       triggerevents([1 2 3 4 ; 5 6 7 8]', [2], 1, 1)
%                   returns  cat(3, [1 2 3], [5 6 7]).
%
%   See also LEADINGEDGES.

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parse Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (ndims(trace) > 2),    error('TRACE can not have more than 2 dimensions.');  end;
if ((length(trigger) > 1) && ~isvectord(trigger)),  error('TRIGGER can not be a matrix.');  end;
if (pre < 0 || post < 0), error('PRE and POST must be non-negative.');  end;

[M,N] = size(trace);
trigger = trigger(:);   % force column vector

if (length(trigger) == M)  % if we haven't been given a trigger list ...
	trigger = find(leadingedges(trigger));    % ... find them
else
	if (~all((trigger >= 1) & (trigger <= M)))
		error('All entries in a TRIGLIST must be between 1 and the # of rows in TRACE.');
	end
    if (nargin > 4), error('GAP can not be used with a TRIGLIST.');  end;
end
P = length(trigger);
Q = pre + post + 1;

if ((nargin > 4) && (gap < 0)),  error('GAP must be non-negative.');  end;

%%%%%%%%%%%%%%%%%%%%%% Filter out events too close %%%%%%%%%%%%%%%%%%%
if (nargin > 4)
	A = 1;  B = 2;   skip = repmat(false, size(trigger));
	while (B <= P)   % skip trigs too close to last unskipped trig
		if ((trigger(B)-trigger(A)) <= gap),  skip(B) = 1;
        else                                  A = B;
		end
		B = B + 1;
	end
	trigger(skip) = [];
	P = length(trigger);  % this might have changed
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%% Extract Events %%%%%%%%%%%%%%%%%%%%%%%%%%
events = repmat(0, [P, Q, N]);    % Force events to be type double

% Do each of these cases separately for efficiency; wasteful to
% check for overruns when we don't need to.
left = trigger <= pre;       right = trigger > M-post;
leftevents  = find(left)';   rightevents = find(right)';
safeevents  = find(~(left | right))';

window = (-pre:post);
for k = leftevents
	inds = window + trigger(k);    mask = (inds < 1);   inds(mask) = 1;
	events(k,:,:) = trace(inds,:);
	events(k,mask,:) = NaN;
end
for k = safeevents
	inds = window + trigger(k);
	events(k,:,:) = trace(inds,:);
end
for k = rightevents
	inds = window + trigger(k);    mask = (inds > M);   inds(mask) = M;
	events(k,:,:) = trace(inds,:);
	events(k,mask,:) = NaN;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Clean Up %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Screwy Matlab indexing: if inds is N x M, trace(inds) will be N x M ...
%   ... unless inds is 1 x M, in which case the dimensions of trace(inds) 
%   are M x 1 if trace is a column vector.  So ...
if (size(events,2) == 1),  events = events';  end