function displaytable(data, colheadings, wid, fms, rowheadings, fid, colsep, rowending)
% Prints formatted matrix of numerical data with headings
% 
% Syntax
% 
% displaytable(data, colheadings, wid, fms, rowheadings, fid, colsep, rowending)
% 
% Input
% 
% data: a matrix or cell array, containing the data to be put in the table.
% If a matrix the numbers in the table will be printed using the default
% format specifier (f), or the specifier(s) supplied in fms. data can also
% be a cell array containing a mixture of strings and numbers. each cell in
% this case can contain only a single scalar value. In this case numbers
% are printed as in the matrix case, while strings are printed up to the
% maximum width of the column.
% 
% colheadings: a cell array of strings for the headings of each column. Can
% be an empty cell array if no headings are required. If no column widths
% are supplied (see below), the columns will be made wide enough to
% accomdate the headings.
% 
% wid: (optional) scalar or vector of column widths to use for the table.
% If scalar, every column will have the same width. If a vector it must be
% of the same length as the number of columns of data. If not supplied, and
% column headers are supplied, the width of the widest column header will
% be used. If not supplied and column headers are not supplied, a default
% with of 16 characters is used.
% 
% fms: (optional) a string, or cell array of strings containing format
% specifiers for formatting the numerical output in each column. If a
% single string, the same specifier is used for every column. If not
% supplied, the 'g' specifier is used for every column.
% 
% rowheadings: (optional) a cell array of strings for the start of each
% row. Can be an empty cell array if no row headings are required. If row
% headings are supplied, the first column will be made wide enough to
% accomodate all the headings.
% 
% fid: (optional) the file id to print to. Use 1 for stdout (to print to
% the command line). 
% 
% colsep: (optional) A string or character to insert between every column.
% The default separation string is ' | ', i.e. a space followed by a
% vertical bar, followed by a space. A table suitible for inclusion in a
% LaTeX document can be created using the ' & ' string, for example.
% 
% rowending: (optional) An optional string or character to be appended at
% the end of every row. Default is an empty string.
% 
% Example 1 - Basic useage
% 
% colheadings = {'number of projects','sales','profit'};
% rowheadings = {'Jimmy Slick', 'Norman Noob'}
% data = [3 rand(1) rand(1); 1 rand(1) rand(1)];
% 
% To format the first number in each row as a decimal (%d), the second
% number %16.4f, and the third as %16.5E do the following:
% 
% wid = 16;
% fms = {'d','.4f','.5E'};
% 
% In this case 16 will be the field width, and '.5E' is what to use for the
% fms argument
% 
% fileID = 1;
% 
% >> displaytable(data,colheadings,wid,fms,rowheadings,fileID);
%             |number of projec |           sales |          profit 
% Jimmy Slick |               3 |          0.4502 |    5.22908E-001 
% Norman Noob |               1 |          0.9972 |    2.78606E-002 
%
% Example 2 - Produce a latex table
% 
% colheadings = {'number of projects','sales','profit'};
% rowheadings = {'Jimmy Slick', 'Norman Noob'};
% data = [3 rand(1) rand(1); 1 rand(1) rand(1)];
% wid = 16;
% fms = {'d'};
% 
% colsep = ' & ';
% rowending = ' \\';
% 
% fileID = 1;
% 
% >> displaytable(data,colheadings,wid,fms,rowheadings,fileID,colsep,rowending);
%
%             & number of projec &            sales &           profit \\
% Jimmy Slick &                3 &    6.948286e-001 &    3.170995e-001 \\
% Norman Noob &                1 &    9.502220e-001 &    3.444608e-002 \\
% 
% Example 3 - Mixed numeric and strings
%
% colheadings = {'number of projects','sales','profit'};
% rowheadings = {'Jimmy Slick', 'Norman Noob'};
% 
% data = {3, rand(1), rand(1); 
%         1, 'none', 0};
%     
% wid = 16;
% fms = {'d','.4f','.5E'};
% 
% fileID = 1;
% 
% >> displaytable(data,colheadings,wid,fms,rowheadings,fileID);
%             | number of projec |            sales |           profit
% Jimmy Slick |                3 |           0.4854 |     8.00280E-001
% Norman Noob |                1 |             none |     0.00000E+000
%
% See Also: DISP, FPRINTF

% Created by Richard Crozier 2012

    if nargin < 8
        rowending = '';
    end

    if nargin < 7
        colsep = ' | ';
    end

    % do some basic checking of input
    if nargin < 6
        % print to the command line
        fid = 1;
    end
    
    if nargin < 5
        % no row headings supplied, use empty cell array
        rowheadings = {};
    end
    
    if nargin < 4
        % no format specifiers supplied, use 'g' for all columns
        fms = 'g';
    end
    
    if nargin < 3
        % default width is 10, this will be modified if column headers are
        % supplied
        wid = 10;
    end
    
    if nargin < 2
        colheadings = {};
    end

    if nargin >= 6 && ~iscellstr(rowheadings)
        error ('row headings must be cell array of strings');
    end
    
    % get the numbers of rows and columns of data
    [nRowD,nColD] = size(data);
    
    % check that sensible format specifiers have been supplied
    if ~iscellstr(fms)
        
        if ischar(fms)
            fms = repmat({fms},1,nColD);
        else
            error('fms must be a string or cell array of strings.');
        end
        
    elseif isempty(fms)
        
        fms = repmat({'f'},1,nColD);
        
    elseif numel(fms) ~= nColD
        
        if numel(fms) == 1
            fms = repmat(fms,1,nColD);
        else
           error('fms does not have the same number of format specifiers as there are columns.');
        end
        
    end
    
    % replace empty format specifiers with 'f'
    for i = 1:numel(fms)
        if isempty(fms{i})
            fms{i} = 'f';
        end
    end
    
    [nRowFms,nColFms] = size(fms);
    if(nRowFms>1)
        error ('fms can not have more than one row');
    end
    

    if ~isempty(rowheadings)
        
        if ~iscellstr(rowheadings)
            error('rowheadings must be a cell array of strings');
        end
        
        if numel(rowheadings) ~= size(data, 1)
            error('Rowheadings must be a cell array of strings of the same size as the number of rows in data.')
        end
    
        rhwid = 0;

        for r = 1:numel(rowheadings)
            % find the maximum width the first column must be to accomodate
            % the row headings
            rhwid = max(rhwid, length(rowheadings{r}));
        end
        
    end
    
    if isscalar(wid)
        
        tempwid = zeros(1, numel(fms));

        for i = 1:numel(fms)

            % get the number of decimal places specified
            [start_idx, end_idx, extents, matches] = regexp(fms{i}, '(?<=\.)\d+');

            if isempty(start_idx)

                % no numbers were found in the format spec, just use the
                % supplied width
                tempwid(i) = wid;

            else

                % some numbers were supplied, use the larger of the width
                % or the number of desired decimal places plus two (to
                % allow for leading number plus decimal point)
                tempwid(i) = max(wid, round(str2double(matches{1})) + 2);

            end

        end

        % replace scalar width with array of width values
        wid = tempwid;
        
    end
        
    if ~isempty(colheadings)
        
        [nRowCH,nColCH] = size(colheadings);
        if(nRowCH>1)
            error ('column headings can not have more than one row');
        end
        
        if ~iscellstr(colheadings)
            error('If not empty, colheadings must be a cell array of strings');
        end
        
        if(nColCH ~= nColD)
            error ('data must have same number of columns as headings');
        end
    
%         fmt = arrayfun(@(x) ['%',num2str(wid(x)),'s |'], 1:nColD, 'UniformOutput', false);

        if ~isempty(rowheadings)
            fprintf(fid, ['%s',colsep], repmat(' ', 1,rhwid));
        end
        
        if nargin < 3
            
            % only data and column headings have been supplied, so
            % determine a sensible value for the width of each column
            
            tempwid = zeros(size(wid));
            for i = 1:numel(colheadings)

                % get a column width which is the minimum to accept the
                % column heading length or the default width specification
                tempwid(i) = max(length(colheadings{i}), wid(i));
                
                if tempwid < 1
                    error('Column width is less than 1, and the column header is empty.')
                end
                
            end
            
            wid = tempwid;
            
        end
        
        % Now loop through the column headings printing them out with the
        % desired column separator
        for i = 1:numel(colheadings)
            
            str = sprintf(['%',num2str(wid(i)),'s'],colheadings{i});

            if i == numel(colheadings)
                % If we are at the end of a row, don't print the column
                % separator
                fprintf(fid, '%s', str(1:wid(i)) );
            else
                % print the column header and column separator
                fprintf(fid, ['%s',colsep], str(1:wid(i)) );
            end
            
        end
        
        fprintf(fid, '%s\n', rowending);

    end
    
    fmt = arrayfun(@(x) ['%',num2str(wid(x)),fms{x}],1:nColD,'UniformOutput',false);
    
    for i = 1:size(data,1)
        
        % first print a row header if necessary
        if ~isempty(rowheadings)
            fprintf(fid, ['%',num2str(rhwid),'s',colsep], rowheadings{i});
        end
            
        % now loop through the data formatting and printing as appropriate
        for j = 1:size(data,2)

            if iscell(data)
                % data is a cell array, possibly containing mixed data
                % types

                if ischar(data{i,j})
                    
                    str = sprintf(['%',num2str(wid(j)),'s'],data{i,j});
                    
                    % we print only what fits in the column width of the
                    % string
                    fprintf(fid, ['%s',colsep], str(1:wid(j)));
                    
                elseif isscalar(data{i,j})
                    
                    % write the number 
                    str = sprintf(fmt{j},data{i,j});
                    
                    if length(str) > wid(j)

                        % convert to scientific notation as it doesn't fit
                        str =  sprintf(['%',num2str(wid(j)),'g'],data(i,j));

                        if length(str) > wid(j)
                            % fill space with #s as the number stil doesn't
                            % fit in
                            str = repmat('#', 1, wid(j));
                        end

                    end
                    
                    if j == size(data,2)
                        % do not print the last column separator at the end
                        % of a row
                        fprintf(fid, '%s', str(1:wid(j)));
                    else
                        fprintf(fid, ['%s',colsep], str(1:wid(j)));
                    end
                    
                else
                    % we can only tabulate strings and scalars, so throw an
                    % error
                    error('CROZIER:displaytable:badcellcontents', ...
                          'each cell in cell array data must contain only individual strings or scalar values.')
                
                end

            else
                
                % data is a numerical matrix

                str = sprintf(fmt{j},data(i,j));

                if length(str) > wid(j)
                    
                    % convert to scientific notation as it doesn't fit
                    str =  sprintf(['%',num2str(wid(j)),'g'],data(i,j));
                    
                    if length(str) > wid(j)
                        % fill space with #s as the number doesn't fit in
                        str = repmat('#', 1, wid(j));
                    end
                    
                end
                
                if j == size(data,2)
                    % do not print the last column separator at the end of
                    % a row
                    fprintf(fid, '%s', str(1:wid(j)));
                else
                    fprintf(fid, ['%s',colsep], str(1:wid(j)));
                end
            end
            
        end
        
        % end of line so put in a new row
        fprintf(fid, '%s\n', rowending);
        
    end

end