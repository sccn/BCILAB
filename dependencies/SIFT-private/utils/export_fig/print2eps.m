%PRINT2EPS  Prints figures to eps with improved line styles
%
% Examples:
%   print2eps filename
%   print2eps(filename, fig_handle)
%   print2eps(filename, fig_handle, options)
%
% This function saves a figure as an eps file, with two improvements over
% MATLAB's print command. First, it improves the line style, making dashed
% lines more like those on screen and giving grid lines their own dotted
% style. Secondly, it substitutes original font names back into the eps
% file, where these have been changed by MATLAB, for up to 11 different
% fonts.
%
%IN:
%   filename - string containing the name (optionally including full or
%              relative path) of the file the figure is to be saved as. A
%              ".eps" extension is added if not there already. If a path is
%              not specified, the figure is saved in the current directory.
%   fig_handle - The handle of the figure to be saved. Default: gcf.
%   options - Additional parameter strings to be passed to print.

% Copyright (C) Oliver Woodford 2008-2012

% The idea of editing the EPS file to change line styles comes from Jiro
% Doke's FIXPSLINESTYLE (fex id: 17928)
% The idea of changing dash length with line width came from comments on
% fex id: 5743, but the implementation is mine :)

% 14/11/2011 Fix a MATLAB bug rendering black or white text incorrectly.
% Thanks to Mathieu Morlighem for reporting the issue and obtaining a fix
% from TMW.

% 8/12/2011 Added ability to correct fonts. Several people have requested
% this at one time or another, and also pointed me to printeps (fex id:
% 7501), so thank you to them. My implementation (which was not inspired by
% printeps - I'd already had the idea for my approach) goes
% slightly further in that it allows multiple fonts to be swapped.

% 14/12/2011 Fix bug affecting font names containing spaces. Many thanks to
% David Szwer for reporting the issue.

% 25/1/2012 Add a font not to be swapped. Thanks to Anna Rafferty and Adam
% Jackson for reporting the issue. Also fix a bug whereby using a font
% alias can lead to another font being swapped in.

% 10/4/2012 Make the font swapping case insensitive.

function print2eps(name, fig, varargin)
options = {'-depsc2'};
if nargin < 2
    fig = gcf;
elseif nargin > 2
    options = [options varargin];
end
% Construct the filename
if numel(name) < 5 || ~strcmpi(name(end-3:end), '.eps')
    name = [name '.eps']; % Add the missing extension
end
% Find all the used fonts in the figure
font_handles = findall(fig, '-property', 'FontName');
fonts = get(font_handles, 'FontName');
if ~iscell(fonts)
    fonts = {fonts};
end
% Map supported font aliases onto the correct name
fontsl = lower(fonts);
for a = 1:numel(fonts)
    f = fontsl{a};
    f(f==' ') = [];
    switch f
        case {'times', 'timesnewroman', 'times-roman'}
            fontsl{a} = 'times-roman';
        case {'arial', 'helvetica'}
            fontsl{a} = 'helvetica';
        case {'newcenturyschoolbook', 'newcenturyschlbk'}
            fontsl{a} = 'newcenturyschlbk';
        otherwise
    end
end
fontslu = unique_bc(fontsl);
% Determine the font swap table
matlab_fonts = {'Helvetica', 'Times-Roman', 'Palatino', 'Bookman', 'Helvetica-Narrow', 'Symbol', ...
                'AvantGarde', 'NewCenturySchlbk', 'Courier', 'ZapfChancery', 'ZapfDingbats'};
matlab_fontsl = lower(matlab_fonts);
require_swap = find(~ismember_bc(fontslu, matlab_fontsl));
unused_fonts = find(~ismember_bc(matlab_fontsl, fontslu));
font_swap = cell(3, 0);
for a = 1:min(numel(require_swap), numel(unused_fonts))
    ind = find(strcmp(fontslu{require_swap(a)}, fontsl));
    n = numel(ind);
    font_swap(1,end+1:end+n) = reshape(mat2cell(font_handles(ind), ones(n, 1)), 1, []);
    font_swap(2,end-n+1:end) = matlab_fonts(unused_fonts(a));
    font_swap(3,end-n+1:end) = reshape(fonts(ind), 1, []);
end
% Swap the fonts
for a = 1:size(font_swap, 2)
    set(font_swap{1,a}, 'FontName', font_swap{2,a});
end
% Set paper size
old_mode = get(fig, 'PaperPositionMode');
set(fig, 'PaperPositionMode', 'auto');
% MATLAB bug fix - black and white text can come out inverted sometimes
% Find the white and black text
white_text_handles = findobj(fig, 'Type', 'text');
M = get(white_text_handles, 'Color');
if iscell(M)
    M = cell2mat(M);
end
M = sum(M, 2);
black_text_handles = white_text_handles(M == 0);
white_text_handles = white_text_handles(M == 3);
% Set the font colors slightly off their correct values
set(black_text_handles, 'Color', [0 0 0] + eps);
set(white_text_handles, 'Color', [1 1 1] - eps);
% Print to eps file
print(fig, options{:}, name);
% Reset the font colors
set(black_text_handles, 'Color', [0 0 0]);
set(white_text_handles, 'Color', [1 1 1]);
% Reset paper size
set(fig, 'PaperPositionMode', old_mode);
% Correct the fonts
if ~isempty(font_swap)
    % Reset the font names in the figure
    for a = 1:size(font_swap, 2)
        set(font_swap{1,a}, 'FontName', font_swap{3,a});
    end
    % Replace the font names in the eps file
    font_swap = font_swap(2:3,:);
    try
        swap_fonts(name, font_swap{:});
    catch
        warning('swap_fonts() failed. This is usually because the figure contains a large number of patch objects. Consider exporting to a bitmap format in this case.');
        return
    end
end
% Fix the line styles
try
    fix_lines(name);
catch
    warning('fix_lines() failed. This is usually because the figure contains a large number of patch objects. Consider exporting to a bitmap format in this case.');
end
return

function swap_fonts(fname, varargin)
% Read in the file
fh = fopen(fname, 'r');
if fh == -1
    error('File %s not found.', fname);
end
try
    fstrm = fread(fh, '*char')';
catch ex
    fclose(fh);
    rethrow(ex);
end
fclose(fh);

% Replace the font names
for a = 1:2:numel(varargin)
    fstrm = regexprep(fstrm, [varargin{a} '-?[a-zA-Z]*\>'], varargin{a+1}(~isspace(varargin{a+1})));
end

% Write out the updated file
fh = fopen(fname, 'w');
if fh == -1
    error('Unable to open %s for writing.', fname2);
end
try
    fwrite(fh, fstrm, 'char*1');
catch ex
    fclose(fh);
    rethrow(ex);
end
fclose(fh);
return
