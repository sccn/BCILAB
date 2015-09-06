function result = wildcardsearch(rootdir, searchstr, casesensitive, strict)
% WILDCARDSEARCH Searches the file system for files matching a wildcard
% pattern.
%    WILDCARDSEARCH(ROOTDIR, SEARCHSTRING) Searches the file system from
%    the starting directory ROOTDIR for files and folders matching the
%    pattern SEARCHSTRING. The function WILDCARDSEARCH is similar to the
%    search file option in Explorer. 
%
%    WILDCARDSEARCH(ROOTDIR, SEARCHSTRING, CASESENSITIVE) The extra boolean 
%    CASESENSITIVE allows for case senstive searching of files and folders. 
%
%    WILDCARDSEARCH(ROOTDIR, SEARCHSTRING, CASESENSITIVE, STRICT) The extra
%    boolean STRICT allows for only precise matches to be returned. The
%    default behaviour is loose. 
%
%    Example:
%     rootdir = 'C:\'
%     searchstr = '*.jpg;*.tif';
%     files = wildcardsearch(rootdir, searchstr);
% 
%    This wil return all JPG and TIF files on the C:\ drive and maybe some
%    folders or files also matching this pattern e.g. 'somefile.jpg.old'
%
%    Example:
%     rootdir = 'C:\'
%     searchstr = 'MatLab';
%     files = wildcardsearch(rootdir, searchstr, true);
% 
%    This wil return only return files or folders matching the pattern 
%    'MatLab'. This could also match a file named 'MatLabRC'.
%
%    Example:
%     rootdir = 'C:\'
%     searchstr = '*\matlab';
%     files = wildcardsearch(rootdir, searchstr, true, true);
% 
%    This wil return only return files named exactly 'matlab'.
%
%    Example:
%     rootdir = 'C:\'
%     searchstr = '*\';
%     files = wildcardsearch(rootdir, searchstr, true, true);
%
%    To find only folders do a strict search and end your search string
%    with the file separator (\ on Windows, / on UNIX). The above example
%    will return a complete directory tree of the C-drive.
%
%    Have a look at the source code for more information. For more 
%    info on this function and how to use it, bug reports, feature
%    requests, etc. feel free to contact the author.
%
%    See also REGEXPDIR

%==========================================================================
% B.C. Hamans (b.c.hamans@rad.umcn.nl)                                 2007
%==========================================================================

% Check input
error(nargchk(2, 4, nargin));

% Create the regular expression
beginstr='('; endstr=')';
if ~exist('strict','var'); strict = false; end
if strict; beginstr=['^' beginstr]; endstr=[endstr '$']; end
if ~exist('casesensitive','var'); casesensitive = false; end
if casesensitive; beginstr = ['(?-i)' beginstr]; end
regexpstr=[beginstr strrep(regexptranslate('wildcard', searchstr), pathsep, [endstr '|' beginstr]) endstr];

% Search
result = regexpdir(rootdir, regexpstr, true);

%==========================================================================
% Changelog:
% 03-09-2007 v1.00 (BCH)  Initial release
%==========================================================================