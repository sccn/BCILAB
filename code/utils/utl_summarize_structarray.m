function [text,results] = utl_summarize_structarray(data,varargin)
% Summarize a struct array's fields in text form.
% Text = utl_summarize_structarray(Data)
%
% In:
%   Data : a struct array to be displayed
%   Options... : additional name-value pairs; possible names are:
%
%                 --- content control ---
%                 'retain_items' : cell array of checker functions (as in cellfun) that all retained
%                                  data items must satisfy (default: {@isvector,@isnumeric,@isfinite})
%
%                 'retain_fields' : cell array of checker functions (as in cellfun) that all 
%                                   retained fields must satisfy (default: {@isvector})
%
%                 'rewrite_fields' : cell array of rewrite rules for field names
%                                   (e.g.{'mcr','Mis-classification rate','auc','Area under Curve'})
%
%                 --- formatting control ---
%
%                 'lhs' : format of the left-hand side of each line (default: '%s')
%
%                 'sep' : separator string between left-hand side and right-hand side (default: ' : ')
%
%                 'rhs' : expression for the right-hand side of each line (in terms of X, the data)
%                         (default: 'sprintf(''%.3f +/- $.3f (N=%d)'',mean(X),std(X),length(X))')
%
%                 'leftspacing' : number of characters to reserve on the left-hand-side; 
%                                minimum marging to first character, and minimum spacing to separator 
%                                (default: [5 30])
%
%                 'rightspacing': number of characters to reserve on the right-hand-side; 
%                                minimum marging from last character, and minimum spacing from 
%                                separator (default: [5 15])
%
% Out:
%   Text : a textual summary, akin to what disp(data(1)) produces, but with means and std. devs
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-11-08

% get the options
hlp_varargin2struct(varargin, 'retain_items', {@(x)isscalar(x),@(x)~isempty(x),@isnumeric,@isfinite}, 'retain_fields',{@(x)~isempty(x)}, ...
    'rewrite_fields',{}, 'lhs','%s', 'sep',' : ', 'rhs','sprintf(''%.2f +/- %.2f (N=%d)'',mean(X),std(X),length(X))', ...
    'leftspacing',[5 0], 'rightspacing', [5 0]);

% get the field names
fnames = fieldnames(data)';
% get the raw contents
raw = struct2cell(data(:));
% reformat the contents
for r=1:size(raw,1)
   conts{r} = {raw{r,:}}; end

% reduce the data items
for r=1:length(conts)
    retain{r} = 1:length(conts{r}); end
for checker = retain_items
    for r=1:length(conts)
        retain{r} = retain{r}(cellfun(checker{1},conts{r}(retain{r}))); end
end
% reduce conts
for r=1:length(conts)
    fullconts{r} = cell(1,length(conts{r}));
    conts{r} = conts{r}(retain{r}); 
    fullconts{r}(retain{r}) = conts{r};
end


% reduce the fields
mask = 1:length(fnames);
for checker = retain_fields
    mask = mask(cellfun(checker{1},conts(mask))); end
fnames = fnames(mask);
fullconts = fullconts(mask);
conts = conts(mask);

% rewrite field names
for f=1:length(fnames)
    fnames{f} = hlp_rewrite(fnames{f},rewrite_fields{:}); end %#ok<USENS>

% create the text
text = cell(1,length(conts));
for r=1:length(conts)
    X = [conts{r}{:}]; %#ok<NASGU>
    left{r} = sprintf(lhs,fnames{r});
    right{r} = eval(rhs);
end

leftspacing = max(leftspacing(2),leftspacing(1)+max(cellfun('length',left))); %#ok<NODEF>
rightspacing = max(rightspacing(2),rightspacing(1)+max(cellfun('length',right))); %#ok<NODEF>
for r=1:length(conts)
    left{r} = [repmat(' ',1,leftspacing - length(left{r})) left{r}];
    right{r} = [right{r} repmat(' ',1,rightspacing - length(right{r})) ];
    text{r} = [left{r} sep right{r}];
end

if nargout == 0
    fprintf('\n');
    % print the text to the console
    for r=1:length(text)
        fprintf([text{r} '\n']); end
    fprintf('\n');
else
    % create results structure
    results.text = text;
    results.fnames = fnames;
    results.conts = conts;
    results.fullconts = fullconts;
end