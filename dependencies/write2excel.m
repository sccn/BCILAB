function write2excel(fileloc,promptforsave,varargin)
% write2excel(fileloc,promtforsave,range1,data1,range2,data2,...)
%
%
%             Uses ActiveX commands to write data_n into range_n in an existing Excel
%             spreadsheet. Inputs (excluding fileloc and promptforsave) must be paired.
%             As of 10/04 update, you may provide the target range (upper left
%             cell to lower right cell) OR just the upper right cell. If
%             the range is specified, the function will verify that the
%             corresponding data block is the correct size, and give an
%             error if not. (This may be useful for error checking, for instance.)
%             If only the upper left cell is provided, write2excel will
%             compute the target range.
%
%             (Please use caution, as you can now overwrite data pretty easily.)
%
%             Additionally, you may now specify cells by address (eg., 'H3') OR row, column
%             (eg, '[3,8]').
%
% FILELOC:            Enter a string representing the location of an Excel file.
%                     Example: 'c:\brett\my archives\test1.xls'
% PROMPTFORSAVE: binary variable. 1 (DEFAULT) = Prompt before saving
%                                 0           = No prompt required
% RANGE SPECIFIER(S): Enter the range(s) to read. You may use Excel
%             cell-references, as in:
%                              'B1:P5'
%                              'B1:B1' (or simply 'B1')
%             OR, alternatively, you may specify the row, column values, as
%             in:
%                              '[8,3]:[12,4]' (to write from row 8, column
%                                     3 to row 12, column 4);
%                              '[8,3]' to write from row 8, column 3
%                                     TO WHATEVER RANGE IS REQUIRED FOR THE DATA BLOCK.
%
% DATA:       NOTE: To enter multiple strings, use cell arrays. Size compatibility is verified
%             if a cell range is given, and is not if only the starting
%             position is specified.
%
% EXAMPLES:   write2excel('c:\brett\my archives\test1.xls', 1, 'C1:E3',magic(3));
%             write2excel('c:\brett\my archives\test1.xls', 0, 'C1:E3',{'string1','string 2', 'string3'});
%             write2excel('c:\brett\my archives\test1.xls', 0, '[3,8]',magic(4));
%
% Written by Brett Shoelson, Ph.D.
% Last update: 1/04.
%              10/04: Allow "dynamic specification" of cell ranges, and
%                     allow specification of cells by row, column format.
%
% SEE ALSO: readfromexcel

if nargin < 4
	msgstr = sprintf('At a minimum, you must specify three input arguments.\nThe first is a string indicating the location of the excel file,\nthe second is a range to be written, and the third contains the data to write.');
	error(msgstr);
elseif ~iseven(nargin-2)
	msgstr = sprintf('Please enter input variables in pairs...\n''write range'',data,''write range'',data')
	error(msgstr)
end
tmp = varargin;
sheetchanges = [];counter = 1;
for ii = 1:length(tmp)
	if ischar(tmp{ii}) & (strcmp(tmp{ii},'sheet') | strcmp(tmp{ii},'sheetname'))
		sheetchanges(counter) = ii;
		counter = counter + 1;
	end
end
if ~isempty(sheetchanges)
	[sheetnames{1:length(sheetchanges)}] = deal(varargin{sheetchanges+1});
end

[pathstr,name,ext] = fileparts(fileloc);
if isempty(ext)
	fileloc = [fileloc,'.xls'];
end
if isempty(pathstr)
	fileloc = which(fileloc,'-all');
	if size(fileloc,1) ~= 1
		error('File was either not located, or multiple locations were found. Please reissue readfromexcel command, providing absolute path to the file of interest.');
	end
end

% Ensure that range sizes and data are size-matched
for ii = 1:2:nargin-2
	if ismember(ii,sheetchanges) | ismember(ii,sheetchanges + 1)
		continue
	end
	% How are cells specified?
	if any(ismember(double(varargin{ii}),[65:90,97:122]))
		addrtype = 'letternumber';
	else
		addrtype = 'rowcol';
	end
	% Is range provided, or should it be auto-calculated?
	autorange = isempty(findstr(varargin{ii},':'));
	switch addrtype
		case 'letternumber'
			if autorange
				r1{ii} = varargin{ii};
				[rx1,cx1] = an2nn(r1{ii});
				rx2 = rx1 + size(varargin{ii+1},1)-1;
				cx2 = cx1 + size(varargin{ii+1},2)-1;
				r2{ii} = nn2an(rx2,cx2);
			else
				tmp = findstr(varargin{ii},':');
				r1{ii} = varargin{ii}(1:tmp-1);
				r2{ii} = varargin{ii}(tmp+1:end);
				[rx1,cx1] = an2nn(r1{ii});
				[rx2,cx2] = an2nn(r2{ii});
			end
		case 'rowcol'
			if autorange
				r1{ii} = varargin{ii};
				[t,r]=strtok(r1{ii},',');
				rx1 = str2num(t(2:end));
				cx1 = str2num(r(2:end-1));
				r1{ii} = nn2an(rx1,cx1);
				rx2 = rx1 + size(varargin{ii+1},1)-1;
				cx2 = cx1 + size(varargin{ii+1},2)-1;
				r2{ii} = nn2an(rx2,cx2);
			else
				tmp = findstr(varargin{ii},':');
				r1{ii} = varargin{ii}(1:tmp-1);
				[t,r]=strtok(r1{ii},',');
				rx1 = str2num(t(2:end));
				cx1 = str2num(r(2:end-1));
				r2{ii} = varargin{ii}(tmp+1:end);
				[t,r]=strtok(r2{ii},',');
				rx2 = str2num(t(2:end));
				cx2 = str2num(r(2:end-1));
				r1{ii} = nn2an(rx1,cx1);
				r2{ii} = nn2an(rx2,cx2);
			end
	end
	if ~autorange % Validate size match for target range, data block
		sz = [rx2 - rx1 + 1, cx2 - cx1 + 1];
		switch class(varargin{ii+1})
			case {'double','cell'}
				sz2 = size(varargin{ii+1});
			case 'char'
				sz2 = [size(varargin{ii+1},1),1];
		end
		if ~isequal(sz,sz2)
			error(sprintf('Mismatched range/data size for input pair %d. Specified range is %d x %d, data block is %d x %d.',(ii+1)/2,sz(1),sz(2),sz2(1),sz2(2)));
		end
	end
end

Excel = actxserver('Excel.Application'); 
Excel.Visible = 0; 
w = Excel.Workbooks; 


try
	excelarchive = invoke(w, 'open', fileloc);
catch
	invoke(Excel, 'quit'); 
	release(w); 
	delete(Excel);
	error(sprintf('Sorry...unable to open file %s',fileloc));
end
Sheets = Excel.ActiveWorkBook.Sheets;

archive = Excel.Activesheet; 
initval = get(archive,'Index');
archive.Unprotect; 

% Read appropriate ranges into output variables
chgcount = 1;
for ii = 1:2:nargin-2
	if ismember(ii,sheetchanges)
		try
			sheet = get(Sheets,'Item',sheetnames{chgcount});
			invoke(sheet,'Activate');
			archive = Excel.Activesheet;
			chgcount = chgcount + 1;
			continue
		catch
			invoke(Excel, 'quit'); 
			release(w); 
			delete(Excel);
			error(sprintf('\nUnable to find/open sheet %s.',sheetnames{chgcount}));
		end
	elseif ismember(ii,sheetchanges + 1)
		continue
	end
	
	archiverange = get(archive, 'Range', r1{ii}, r2{ii}); 
	set(archiverange, 'value', varargin{ii+1}); 
	release(archiverange);
end

sheet = get(Sheets,'Item',initval);
invoke(sheet,'Activate');

if ~promptforsave
	invoke(excelarchive,'save');
end
invoke(Excel, 'quit'); 
release(excelarchive); 
release(w); 
delete(Excel);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SUBFUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function k=iseven(x)
k = x/2==floor(x/2);
return

function [r, c] = an2nn(cr)
% convert alpha, number format to number, number format
t = find(isletter(cr)); 
t2 = abs(upper(cr(t))) - 64; 
if(length(t2) == 2), t2(1) = t2(1) * 26; end
c = sum(t2); r = str2num(cr(max(t) + 1:length(cr)));
return

function cr = nn2an(r, c)
% convert number, number format to alpha, number format
%t = [floor(c/27) + 64 floor((c - 1)/26) - 2 + rem(c - 1, 26) + 65]; 
t = [floor((c - 1)/26) + 64 rem(c - 1, 26) + 65]; 
if(t(1)<65), t(1) = []; end
cr = [char(t) num2str(r)]; 