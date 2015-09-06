function printmat(mat, mmask, dmmask, ngmmask, fmmask, indent)

%PRINTMAT Print SSMAT/SSDIST/SSFUNC content with optional indentation.
%   PRINTMAT(mat, mmask, dmmask, ngmmask, fmmask, indent)

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

MAXCOL      = 90;

if nargin < 6, indent = 1; end
if nargin < 5 || isempty(fmmask), fmmask = zeros(size(mat)); end
if nargin < 4 || isempty(ngmmask), ngmmask = zeros(size(mat)); end
if nargin < 3 || isempty(dmmask), dmmask = false(size(mat)); end
if nargin < 2 || isempty(mmask), mmask = false(size(mat)); end

AllInteger  = isequal(round(mat), mat);
if AllInteger, COL = 5; else COL = 9; end
MAXEL       = floor((MAXCOL-4*indent)/COL);

k   = 0;
while size(mat, 2) > MAXEL*k
    columns = MAXEL*k+1 : min([MAXEL*(k+1) size(mat, 2)]);
    if size(mat, 2) > MAXEL
        if isscalar(columns), fprintf(1, '\tColumn %d\n', columns);
        else fprintf(1, '\tColumns %d through %d\n', columns(1), columns(end)); end
        fprintf(1, '\n');
    end
    for i = 1 : size(mat, 1)
        fprintf(1, '\t');
        for j = columns
            if fmmask(i, j) > 0
                fline   = sprintf('FN%d', fmmask(i, j));
                fprintf(1, '%s%s', repmat(' ', 1, COL-length(fline)), fline);
            elseif ngmmask(i, j) > 0
                fline   = sprintf('NG%d', ngmmask(i, j));
                fprintf(1, '%s%s', repmat(' ', 1, COL-length(fline)), fline);
            elseif dmmask(i, j), fprintf(1, '%sDYN', repmat(' ', 1, COL-3));
            elseif mmask(i, j), fprintf(1, '%sVAR', repmat(' ', 1, COL-3));
            else
                if round(mat(i, j)) == mat(i, j), fline = sprintf(' %%%dd', COL-1);
                else fline = sprintf(' %%%d.4f', COL-1); end
                fprintf(1, fline, mat(i, j));
            end
        end
        fprintf(1, '\n');
    end
    k   = k + 1;
    fprintf(1, '\n');
end
fprintf(1, '\n');

