function n = progress(rate,title)
%PROGRESS   Text progress bar
%   Similar to waitbar but without the figure display.
%
%   Start:
%   PROGRESS('init'); initializes the progress bar with title 'Please wait...'
%   PROGRESS('init',TITLE); initializes the progress bar with title TITLE
%   PROGRESS(RATE); sets the length of the bar to RATE (between 0 and 1)
%   PROGRESS(RATE,TITLE); sets the RATE and the TITLE
%   PROGRESS('close'); (optionnal) closes the bar
%
%   Faster version for high number of loops:
%   The function returns a integer indicating the length of the bar.
%   This can be use to speed up the computation by avoiding unnecessary
%   refresh of the display
%   N = PROGRESS('init'); or N = PROGRESS('init',TITLE);
%   N = PROGRESS(RATE,N); changes the length of the bar only if different
%   from the previous one
%   N = PROGRESS(RATE,TITLE); changes the RATE and the TITLE
%   PROGRESS('close'); (optionnal) closes the bar
%
%   The previous state could be kept in a global variable, but it is a bit
%   slower and doesn't allows nested waitbars (see examples)
%
%   Known bug: Calling progress('close') shortly afer another call of the
%   function may cause strange errors. I guess it is because of the
%   backspace char. You can add a pause(0.01) before to avoid this.
%
%   Examples:
%       progress('init');
%       for i=1:100
%           progress(i/100, sprintf('loop %d/100',i));
%
%           % computing something ...
%           pause(.1)
%       end
%       progress('close'); % optionnal
%
%
%       % Inside a script you may use:
%       n = progress('init','wait for ... whatever');
%       for i=1:100
%           n = progress(i/100,n);
%
%           % computing something ...
%           pause(.1)
%       end
%       progress('close');
%
% 
%       % Add a time estimation:
%       progress('init','Processing...');
% 		tic       % only if not already called
% 		t0 = toc; % or toc if tic has already been called
% 		tm = t0;
% 		L  = 100;
% 		for i=1:L
% 			tt = ceil((toc-t0)*(L-i)/i);
% 			progress(i/L,sprintf('Processing... (estimated time: %ds)',tt));
% 
% 			% computing something ...
% 			pause(.1)
% 		end
% 		progress('close');
% 
% 
%       % Add a faster time estimation:
% 		n  = progress('init','Processing...');
% 		tic       % only if not already called
% 		t0 = toc; % or toc if tic has already been called
% 		tm = t0;
% 		L  = 100;
% 		for i=1:L
% 			if tm+1 < toc % refresh time every 1s only
% 				tm = toc;
% 				tt = ceil((toc-t0)*(L-i)/i);
% 				n  = progress(i/L,sprintf('Processing... (estimated time: %ds)',tt));
% 			else
% 				n  = progress(i/L,n);
% 			end
% 
% 			% computing something ...
% 			pause(.1)
% 		end
% 		progress('close');
%
%       % Nested loops:
%       % One loop...
% 		n1 = progress('init','Main loop');
% 		for i=0:7
% 			n1 = progress(i/7,n1);
% 
% 			% ... and another, inside the first one.
% 			n2 = progress('init','Inside loop');
% 			for j=0:50
% 				n2 = progress(j/50,n2);
% 
% 				% computing something ...
% 				pause(.01)
% 			end
% 			progress('close');
% 		end
% 		pause(.01)
% 		progress('close');

%   31-08-2007
%   By Joseph martinot-Lagarde
%   joseph.martinot-lagarde@m4x.org

%   Adapted from:
%   MMA 31-8-2005, martinho@fis.ua.pt
%   Department of Physics
%   University of Aveiro, Portugal

%% The simplest way to bypass it...
% n = 0; return

%% Width of the bar
%If changes are made here, change also the default title
lmax=70;  % TM: changed from lmax=50;

%% Erasing the bar if necessary
% not needed, but one could find it prettier
if isequal(rate,'close')
	% there were 3 '\n' added plus the title and the bar itself
	fprintf(rep('\b',2*lmax+3))
	return
end

%% The init
if isequal(rate,'init') % If in init stage
	cont = 0;           % we don't continue a previous bar
	back = '\n';        % add a blank line at the begining
	rate = 0;           % start from 0
else
	cont = 1;           % we continue a previous bar
end

%% No need to update the view if not necessary
% optional, but saves a LOT of time

% length of the displayer bar in number of char
% double([0,1]) to int([0,lmax-1])
n = min(max( ceil(rate*(lmax-2)) ,0),lmax-2);

% If the 2nd arg is numeric, assumed to be the previous bar length
if nargin >=2 && isnumeric(title)
	if n == title % If no change, do nothing
		return
	else          % otherwise continue
		n_ = title;
		clear title
	end
else % draw the whole bar
	n_ = -1;
end

%% The title
% If a new title is given, display it
if exist('title','var')
	Ltitle = length(title);
	if Ltitle > lmax % If too long, cut it
		title = [title(1:lmax) '\n']
	else             % otherwise center it
		title = [rep(' ',floor((lmax-Ltitle)/2)) title rep(' ',ceil((lmax-Ltitle)/2)) '\n'];
	end
	if cont % If not in init stage, erase the '\n' and the previous title
		back = rep('\b',lmax+1);
	end
else
	if cont % If not in init stage, give a void title
		title = '';
		back  = ''; % has to be set
	else    % else set a default title
		title = '                  Please wait...                  \n';
	end
end

%% The bar
% '\f' should draw a small square (at least in Windows XP, Matlab 7.3.0 R2006b)
% If not, change to any desired single character, like '*' or '#'
if ~cont || n_ == -1 % at the begining disp the whole bar
	str = ['[' rep('*',n) rep(' ',lmax-n-2) ']\n'];
	if cont % If not in init stage, erase the previous bar
		back = [back, rep('\b',lmax+1)];
	end
else % draw only the part that changed
	str  = [rep('*',n-n_) rep(' ',lmax-n-2) ']\n'];
	back = [back, rep('\b',lmax-n_)];
end

%% The print
% Actually make the change
fprintf([back title str]);
return

%% Function to repeat a char n times
function cout = rep(cin,n)
if n==0
	cout = [];
	return
elseif length(cin)==1
	cout = cin(ones(1,n));
	return
else
	d    = [1; 2];
	d    = d(:,ones(1,n));
	cout = cin(reshape(d,1,2*n));
	return
end