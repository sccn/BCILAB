%% introduction

% CPRINTF	us 11-Jun-2009 08:58:55 / us@neurol.unizh.ch
%
%		- converts an ND array of any MATLAB data type
%		  to a 2D character array
%		- the input may be a cell array formatted as a table with
%		  row/column labels and row/column separators
%		- any input can be formatted as a table using
%		  any combination of
%		  row/column labels and row/column separators
%
%		  use ML's web browser with html proportional text
%		  set to a fixed font to see the HTML file with
%		  proper formatting
%% --- syntax

%		[T,TC,AC,P] = CPRINTF(A,OPT1,...,OPTn)
%				converts A
%		 P          = CPRINTF
%				returns  the engine parameters
%		              CPRINTF
%				displays the help
%% --- input

% A	:	an ND array of
%  		- real and/or complex full   numeric data
%  		- real and/or complex sparse numeric data
%  		- logical data
%  		- char strings
%		- structures
%  		- other objects
%  		an ND cell array of any combination of the above
%% --- output

% T	:	character array with same number of rows and
%		delimiter separated columns as A
% TC	:	cell array of T
% AC	:	cell array of A (useful if A is displayed as a table)
% P	:	structure with engine parameters
%
% NOTE	:	unless an output argument is present, CPRINTF
%		DOES print the result even if an [;] is added
%		at the end to suppress printing
%% --- summary of options

% OPTION	argument	description		default
% ----------------------------------------------------------------
%		SC		a single CHAR
%		CS		a CHAR string
%		FS		a format spec
% ----------------------------------------------------------------
% CELL				data type
% ----------------------------------------------------------------
%    -c	:	FS		character string	'%s'
%    -n	:	FS		numeric real		'%g'
%   -cr	:	FS		numeric complex real	'%g'
%   -ci	:	FS		numeric complex imag	'%+gi'
%   -cd	:	FS		numeric complex delim	' '
%    -s	:	FS		numeric sparse indices	'(%g %g)'
%    -f	:	FS		false			'F'
%    -t	:	FS		true			'T'
%    -E	:	FS		empty CELL		'E(class)'
%    -I	:	FS		±Inf			'±INF'
%    -N	:	FS		NaN			'NAN'
%   -hs	:	T|F		convert to single hex	[F]
%   -hd	:	T|F		convert to double hex	[F]
%  -nex	:	T|F		no char CELL expansion	[F]
%    -O	:	FS		other objects		[built-in]
%   -Or	:	T|F		other objects raw mode	[built-in]
%    -C	:	FS		text surrounding CELLs	'%s'
%   -la	:	T|F		cell content alignment	[F]
%  -cla	:	T|F		complex alignment	[F]
%				F = right align
%				T = left  align
%
% ROW				content
% ----------------------------------------------------------------
%    -L	:	CS		leading   row text	''
%    -T	:	CS		trailing  row text	''
%    -d	:	SC		separator between CELLs	' '
%   -dt	:	SC		separator table columns	' '
%   -nd	:	T|F		show ND page indices	[F]
%
% TABLE				content / processing
% ----------------------------------------------------------------
%   -Ct	:	FS		text surrounding body	'%s'
%				but not label CELLs
%   -nc	:	FS		numeric real col	'%g'
%   -nr	:	FS		numeric real row
%   -Lh	:	{tn}		table  name		{' '}
%   -Lc	:	{c1...cn}	column labels
%   -Lr	:	{r1...rn}	row    labels
%  -Lcs	:	SC		column separator	''
%  -Lrs	:	SC		row    separator	''
%   -it	:	T|F		input is a       table	[F]
%   -mt	:	T|F		convert input to table	[F]
%   -ic	:	T|F		column width		[F]
%				F = max of all  cols
%				T = max of each col
%
% OUTPUT FILES
% ----------------------------------------------------------------
%   -fa	:	name		append to output file	[]
%   -fc	:	name		create    output file	[]
%   -fi	:	name		input  file		[]
%   -fm	:	marker		insert result at marker	[]
%				in file [-fin]
%   -fr	:	{t1,r1,...}	replace token tx with	[]
%				value rx
%
% PROCESSING
% ----------------------------------------------------------------
%    -p	:	T|F		do NOT use preferences	[F]
%  -opt	:	struct		use struct.option	[]
%  -ini	:	name		read options from file	[]
%  -sav	:	name		save options to   file	[]
%  -tab	:	n		use n SPACES/TAB	[8]
% -ntab	:	T|F		keep TABs in CELLs	[F]
%    -q	:	n		do not display result	[F]
%   -so	:	T|F		collect all output in	[F]
%				a structure
%   -db	:	T|F		show processing stages	[F]
%% arrays and tables
%% --- definitions

% each element of an array of any data type including a cell array
%    is considered a CELL
% CELLs can be formatted in various ways according to the data type
% a separator is printed between two CELLs
%
% CPRINTF optionally prints the input as a table consisting of
%    default or user definable
%    - table name
%    - row and/or column labels
%    - row and/or column separators
%    in any combination
% any input of min size 3x3 may be	printed   as a table	[-it]
% any input of any size may be		converted to a table	[-mt]
%
% CPRINTF optionally prints LEADING and/or TRAILING text
%    before and after each row of an array or a table
%% --- formatting

% sequence of formatting
%
% 1: CELL			each CELL's content is displayed according to
%				   the various CELL formatting options
%				by default, CELLs of Mx1 or ND character strings
%				   are expanded unless the [-nex] option is used
%				by default, all column widths of the array are
%				   the maximun width of all columns after conversion
%				   unless the [-ic] option is used, in which case
%				   each column's max width is used
%				a CELL's text may be left aligned [-la] with
%				   respect to the column width
% 2: CELL surroundings		each CELL's output may be embedded in additional
%				   text options [-C]/[-Ct]
% 3: CELL separators		a separator is inserted between two CELLs according
%				   to the [-d] option
% 4: ROW/COL labels		each label's numeric content is displayed according
%				   to table [-nc]/[-nr] options
% 5: ROW/COL separators		separators are displayed according to
%				   various options and depend 
%				   on the content of the separator as well as
%				   the content of the individual column separator
%				   of a predefined table [-it]
% 6: table name			is added at position (1,1) and formatted according
%				   to the common [-n] option
% 7: L/T text			row text delimiters are added
%				   before (L) and/or after (T) the array or table
%
% options to create a format string
%
% 1:				an ANSI C format string, eg,
%				      'a number %g'
% 2:				an anonymous function returning a string, eg,
%				      @(x) sprintf('a number %10d',x)
% 3:				a handle to a user defined function
%				   returning a string, eg,
%				      @mynumber
%				   where mynumber.m is a function file
% 4:				a handle to a ML stock function
%				   returning a string, eg,
%				      @dec2bin
%% --- arrays

% NOTE:	all data is collected in structure S to not clutter the workspace
%%

% simple real numeric matrix

	clear s;
	s.r=reshape(1:4*6,[4,6]);

	disp(s.r);
%%

% CPRINTF default output
% - all CELL widths have the same max CELL width
% - a CELL separator [-d] is inserted between CELLs 2:end-1 including
%      row labels and separators if the output is a table

% NOTE:	the [-tab] option (see section -engine options-) has NO effect
%	on the column separator [-d] if it contains TABs

	cprintf(s.r);
%% --- creating tables using row/column options

% label and separator options
%-------------------------------------------------------------------------------
% -Lh		name of the table at position (1,1)	{name}
%		- numeric data is converted
%		  according to the [-n] option
%		- names are always left aligned
% -Lc		col labels				{cl1,...,cln}
% -Lr		row labels				{rl1,...,rln}
%		- label(s) will be repeated or cut
%		  automatically to fit the table size
% -nc		format of numeric col labels		FS		['%g']
% -nr		format of numeric row labels		FS		['%g']
% -Lcs		col  separator				'csep'
% -Lrs		row  separator				'rsep'
% -d		CELL separator				'csep'
%		- also apply to labels and
%		  row/col separators
%		- a '\t' is printed as TAB
% -dt		CELL separator				'csep'
%		- only apply to col labels
%		  and separators in a table
%		- a '\t' is printed as TAB

% tables are constructed according to these input options

	s.clab={'apple','orange','banana','kiwi','peach','pear','NOTUSED'};
	s.rlab={'sweet','sour','salty','bitter'};

	cprintf(s.r,...
		'-Lc',	s.clab);
%%
	cprintf(s.r,...
		'-Lr',	s.rlab);
%%
	cprintf(s.r,...
		'-Lr',	s.rlab,...
		'-Lc',	s.clab,...
		'-Lcs',	'=');
%%
	cprintf(s.r,...
		'-Lr',	s.rlab,...
		'-Lc',	s.clab,...
		'-Lrs',	'|',...
		'-Lcs',	'=',...
		'-Lh',	{'A table'});
%%

% alignment of CELL content
% by default, each CELL's content is right aligned within the cols max width
%-------------------------------------------------------------------------------
% -la:		left align				T|F	[F]

	cprintf(s.r,...
		'-Lr',	s.rlab,...
		'-Lc',	s.clab,...
		'-Lrs',	'|',...
		'-Lcs',	'=',...
		'-la',	1);
%%

% widths of CELL content
% by default, the width of all CELLs is the max of all CELLs' width
%-------------------------------------------------------------------------------
% -ic		adjust each column width according	T|F	[F]
%		to its max CELL width

	cprintf(s.r,...
		'-Lr',	s.rlab,...
		'-Lc',	s.clab,...
		'-Lrs',	'|',...
		'-Lcs',	'=',...
		'-ic',	1);
%%

% CELL separators
%-------------------------------------------------------------------------------
% -d		CELL separator				'csep'
% -dt		CELL separator of table body CELLs only	'csep'

	cprintf(s.r,...
		'-Lr',	s.rlab,...
		'-Lc',	s.clab,...
		'-Lrs',	'|',...
		'-Lcs',	'=',...
		'-d',	'<.>');				% global CELL separator
%%

% using different separators in a table

	cprintf(s.r,...
		'-Lr',	s.rlab,...
		'-Lc',	s.clab,...
		'-Lrs',	'#',...
		'-Lcs',	'=',...
		'-dt',	'|',...				% table body CELLs' separator
		'-d',	'.');				% global CELL separator
%%

% if the separators contain TAB(s), CELL separators will not be continuous
%    for reasons of consistency

	cprintf(s.r,...
		'-Lr',	s.rlab,...
		'-Lc',	s.clab,...
		'-Lrs',	'|',...
		'-Lcs',	'=',...
		'-dt',	'|',...				% table body CELLs' separator
		'-d',	'\t');				% global CELL separator
%%

% formatting of numeric row/column labels (see also [-n] below)

	s.cclab=s.clab;
	s.rrlab=s.rlab;
	s.cclab([1,2])={1,20};
	s.rrlab([1,4])={pi,-4*pi};
	s.cfmt='col:%d';				% a char spec
	s.rfmt=@(x) sprintf('row:%d/%8.4f',sign(x),x);	% a function handle

	cprintf(s.r,...
		'-Lr',	s.rrlab,...
		'-Lc',	s.cclab,...
		'-Lrs',	'|',...
		'-Lcs',	'=',...
		'-nc',	s.cfmt,...
		'-nr',	s.rfmt);
%% --- creating tables using the default format

% an array can be displayed as a table using default options
%-------------------------------------------------------------------------------
% -mt		create a table				T|F	[F]

	cprintf(s.r,...
		'-mt',	1);
%%

% any row [-Lr]/[-Lrs] and/or column option [-Lc]/[-Lcs] overrides
%    the default option

	cprintf(s.r,...
		'-mt',	1,...
		'-Lcs',	'=',...
		'-Lr',	s.rlab,...
		'-Lh',	{'TBL.1'});
%% --- printing arrays that contain a predefined table

% an input cell array, which is already in the format of a complete table,
%    can be displayed using the auto-formatting option
%-------------------------------------------------------------------------------
% -it		assumes the input is a complete table	T|F	[F]
%		consisting of
%		- col labels and separators
%		- row labels and separators

	s.tbl={
%	row	row
%	labels	separators
	'tbl',	'|',		'C1','C2','C3','C4','C5','C6'	% col labels
	'=',	'=',		'=','-','=','sep','=','-'	% col separators
	'R1',	'sep',		1,2,3,4,5,6
	'R2',	'#',		1,2,3,4,5,6
	};

% raw output of the cell

	disp(s.tbl);
	cprintf(s.tbl);
%%

% output as a table

% NOTE:	single char col separators are repeated to the width of each col
%	multi  char col separators are printed according to the [-la] option
%	the [-ic] CELL width option may be applied

	cprintf(s.tbl,...
		'-it',	1);
%%

% any row [-Lr]/[-Lrs] and/or column option [-Lc]/[-Lcs] overrides
%    the default option

	cprintf(s.tbl,...
		'-it',	1,...
		'-Lr',	{'a','b'},...
		'-Lcs',	'*');
%% --- leading and trailing text

% text can be added before or after an array or table
%-------------------------------------------------------------------------------
% -L		leading  text				character string
% -T		trailing text				character string

	cprintf(s.r,...
		'-L',	'row\t');
%%
	cprintf(s.tbl,...
		'-it',	1,...
		'-L',	'my rows\t',...
		'-T',	'\teor');
%% CELL conversion options

% NOTE: the [-mt] option is used throughout for ease of visibility
%% --- 2D character arrays and ND cell arrays of ND character strings

% by default, Mx1 or ND character strings are expanded to 2D arrays
%    before conversion
%-------------------------------------------------------------------------------
% -nex		do NOT expand				T|F	[F]
%		use [-O] option display style

% a simple array of character strings

	s.sc=[
		'this is a'
		'  test   '
	];

% white spaces are kept and the [-la] option has NO effect

	cprintf(s.sc,...
		'-mt',	1,...
		'-la',	1);				% left align
%%

% a cell array of 1xN character strings

	s.cs={
		'this is'	'yet another'
		'test'		'    example'		% (2,2) has SPACEs!
	};

	cprintf(s.cs,...
		'-mt',	1);
%%

% NOTE:	the [-la] option has an effect EXCEPT for embedded SPACEs

	cprintf(s.cs,...
		'-mt',	1,...
		'-la',	1);				% left align
%%

% a cell array of mixed MxN character strings

	s.ndcs={'foo','GOO'.';'HOO'.','ioo'};

	disp(s.ndcs);
%%

% NOTE:	placeholders for 1xN character arrays

	cprintf(s.ndcs,...
		'-mt',	1,...
		'-la',	1,...				% left align
		'-d',	'\t');				% delimiter: TAB
%%
	s.ndsc={
		s.sc,		[s.sc;s.sc;s.sc],	s.sc
		s.sc,		'FOO',			[s.sc;s.sc]
		'XXX'.',	s.sc			[s.sc;s.sc]
	};

	disp(s.ndsc);
%%
	cprintf(s.ndsc,...
		'-mt',	1,...
		'-la',	1,...				% left align
		'-d',	'   ');				% delimiter: 3 SPACEs
%%

% a ND cell array of MxN character strings

	s.ndnd(:,:,1)=s.ndsc;
	s.ndnd(:,:,2)=s.ndsc;
	
	disp(s.ndnd);
%%
	cprintf(s.ndnd,...
		'-mt',	1,...
		'-dt',	'|');
%%

% do not expand non-1xN character strings

	cprintf(s.ndnd,...
		'-mt',	1,...
		'-nex',	1);
%% --- full mixed real and/or complex numeric arrays

% formatting options for numeric data (see section -formatting-)
%-------------------------------------------------------------------------------
% -n		real numbers				FS	['%g']
% -cr		real      part of complex numbers	FS	['%g']
% -ci		imaginary part of complex numbers	FS	['%+gi]
% -cd		delimiter between real/imaginary part	FS	[' ']
% -cla		left align imaginary parts		T|F	[F]

	s.m=[
		1,		2,		3+3i
		nan-4.2i,	nan,		inf
		-inf,		inf+3i*inf,	inf+5.5i
	];

	disp(num2cell(s.m));
%%

% unlike DISP, CPRINTF prints real numbers
%    at the position of the real part of the complex number
%    for reasons of consistency
% column widths of real and imaginary parts are adjusted to the respective max
%    width of each component (see c:2)

	cprintf(s.m,...
		'-mt',	1);
%%

% left align imaginary parts of complex numbers [-cla]

	cprintf(s.m,...
		'-mt',	1,...
		'-cla',	1);
%%

% using various numeric options

% NOTE:	only real numbers use the [-n] option
%	the [-cd] separator is applied to ALL numbers

	cprintf(s.m,...
		'-mt',	1,...
		'-n',	@(x) sprintf('%g=%s',x,num2hex(single(x))),...
		'-cr',	'R %g',...
		'-ci',	'I %+gi',...
		'-cd',	' / ');
%% --- sparse mixed real and/or complex numeric arrays

% sparse arrays are displayed in accordance with ML's DISP functionality
% col 1		indices
% col 2		values
%
% formatting options for numeric data (see section -formatting-)
%-------------------------------------------------------------------------------
% -s		sparse indices must take TWO args	FS	[('%g %g)']

	s.s=sparse(4,4);

% output of empty sparse arrays show max indices with a value of zero

	cprintf(s.s,...
		'-mt',	1);
%%
	s.s(1,1)=nan;
	s.s(1,2)=inf;
	s.s(1,3)=-inf;
	s.s(1,4)=pi;
	s.s(2,1)=4+3i;
	s.s(2,2)=5+4i*inf;
	s.s(2,3)=6+5i*nan;

	disp(s.s);
%%
	cprintf(s.s,...
		'-mt',	1);
%%

% using various options

	cprintf(s.s,...
		'-mt',	1,...
		'-ic',	1,...				% use individual col width
		'-n',	'    value: %6.2f',...		% real numbers
		'-cr',	'    value: %6.2f',...		% real part of complex numbers
		'-s',	@(x,y) sprintf('R=%d C=%d',x,y));
%% --- logical arrays

% formatting options for logical data (see section -formatting-)
%-------------------------------------------------------------------------------
% -t		TRUE  values				FS	['T']
% -f		FALSE values				FS	['F']

	cprintf(real(s.m)>4,...
		'-mt',	1);
%%
	cprintf(imag(s.m)<1,...
		'-mt',	1,...
		'-f',	'NO',...
		'-t',	'yes');

%% --- NaN/Inf

% formatting options for NaN/Inf (see section -formatting-)
%-------------------------------------------------------------------------------
% -N:		NaN      values				FS	[ML def]
% -I:		-Inf/Inf values				FS	[ML def]
%
% NOTE:	[-N] and/or [-I] options override ALL other options for numeric
%	data types for reasons of consistency

	cprintf(s.m,...
		'-mt',	1,...
		'-d',	'\t');				% use TAB as CELL separator
%%

% using various options

	cprintf(s.m,...
		'-mt',	1,...
		'-d',	'\t',...			% use TAB as CELL separator
		'-n',	'%5.2f',...			% real numbers
		'-N',	'#fff#',...			% NaNs, eg, for a spreadsheet
		'-I',	@(x) sprintf('%dI',sign(x)),...	% Infs
		'-ci',	'%5.2fi');			% imaginary part of complex numbers
%% --- hexadecimal representation of numeric input

% conversion to singles/doubles IEEE hexadecimal strings
%-------------------------------------------------------------------------------
% -hs		all numeric CELLs are converted to
%		   single hex				T|F	[F]
% -hd		   double hex				T|F	[F]

	cprintf(s.m,...
		'-mt',	1,...
		'-hs',	1);
%% cell arrays

% ALL CELLs of a cell array are formatted according to the above formatting rules
% if a CELL contains complex data, ALL numbers will be displayed using
%    the complex spacing mode (see above)
% if a CELL contains sparse data, it is converted to full
%    prior to conversion
% a cell array of any (mixed) content can be displayed as a table

	s.ss=sparse(1,1,-pi);
	s.cc={
		'this',		'is',	false,	2-3i,	nan
		'a cell',	pi,	1000,	s.ss,	true
	};

	cprintf(s.cc,...
		'-mt',	1);
%% --- empty CELLs and other data types

% formatting options for empty CELLs and other data types
%    (see section -formatting-)
%-------------------------------------------------------------------------------
% -E		empty CELLs				FS	['E(class)']
%		EXCEPT empty character strings
% -O		other data types			FS	[built-in]
%		including ND arrays of any type
%		default output according to the input
%		logical:	L(dim:size:class)
%		numeric:	N(dim:size:class.g.s.c)
%		      g:	is global
%		      s:	is sparse
%		      c:	is complex
%		cell   :	C(dim:size:class)
%		struct :	S(dim:size:class)
%		func   :	F(@function)
%		other  :	O(dim:size:class)
% -Or		only show content within brackets, eg,	T|F	[F]
%		to save function handle definitions
%		in an ASCII file

% NOTE:	the [-O] option will be applied to ALL other data types
%	the engine structure returns a function handle, which
%	may be used to check conversion (see section -macros-)
%
%		p=cprintf;	% returns the engine structure
%		r=p.other(magic(3)+4i)
%		r = N(2:3x3:double.c)
%
%	users can write their own functions, which must return
%	a 1xN character string
%
%		cprintf(...,'-O',@myother)

% NOTE:	difference of default output of various empty/other data

	cprintf([],...
		'-Lr',	{'empty array             []'},...
		'-mt',	1);

	cprintf({[]},...
		'-Lr',	{'cell with empty array {[]}'},...
		'-mt',	1);

	cprintf({},...
		'-Lr',	{'empty cell              {}'},...
		'-mt',	1);

	cprintf({{}},...
		'-Lr',	{'cell with empty cell  {{}}'},...
		'-mt',	1);
%%
	s.cc(1,1)={@(x) x};
	s.cc(1,3)={[]};
	s.cc(2,3)={sparse(3,3,-pi)};
	s.cc(2,1)={repmat(struct('a','b'),[2,3,3])};
	s.cc(3,1)={cell([3,3,3])};
	s.cc(3,2)={s.cc};
	s.cc(3,3)={repmat(struct('a',{'a'}),[2,2,3])};

	disp(s.cc);
%%
	cprintf(s.cc,...
		'-mt',	1);
%%

% using [-E]/[-O] options

	s.O=@(x) sprintf('o:%d/%s',ndims(x),class(x));

	cprintf(s.cc,...
		'-mt',	1,...
		'-E',	'empty',...			% EMPYT CELLs
		'-O',	s.O);				% OTHER CELL content
%% --- surrounding text of CELLs

% CELLs may be displayed with an additional surrounding text
%-------------------------------------------------------------------------------
% -C		surrounding of EACH CELL		FS	['%s']
%		including labels and separators
% -Ct		surrounding of table body CELLs only	FS	['%s']

% NOTE:	FS must only contain %s format specs

% embedding all CELLs

	cprintf(s.m,...
		'-mt',	1,...
		'-C',	'<%s>');
%%

% only embedding table body CELLs

	cprintf(s.m,...
		'-mt',	1,...
		'-Ct',	'{%s\t}');
%% ND arrays

% ND arrays are displayed according to ML's display rules
%    all arrays are converted to their 2D equivalent
%

	s.nda=reshape(1:2*4*2*2*2,[2,4,2,2,2]);

	disp(s.nda);
%%
	cprintf(s.nda,...
		'-mt',	1);
%% ND cell arrays of mixed content

% ND cell arrays are converted according to their CELL content
%    the cell array is converted to its 2D equivalent

	s.cmc=reshape(num2cell(1:2*2*2*2),[2,2,2,2]);
	s.cmc(1,1,1,1)={magic(3)};
	s.cmc(2,2,1,1)={s.ndcs};
	s.cmc(1:2,1:2,2,1)=s.ndcs;
	s.cmc(1,1,1,2)={'foo'};
	s.cmc(1:2,2,1,2)={false,true};
	s.cmc(1,1:2,2,2)={nan,-inf};
	s.cmc(2,1:2,2,2)={-100-2i,100+2i};

	disp(s.cmc);
%%

% NOTE:	display of vertical character strings
%	with proper spacing of other content
%	all numbers are displayed in complex style

	cprintf(s.cmc,...
		'-mt',	1,...
		'-f',	'NO',...
		'-t',	'YES',...
		'-N',	'a nan',...
		'-Ct',	'<%s>',...
		'-dt',	' | ',...
		'-d',	'\t');
%% --- page indexing

% a page index may be displayed, which adheres to the common ML display rules
%-------------------------------------------------------------------------------
% -nd		show page indices			T|F	[F]

% NOTE:	if the [-Lr] option is used, only data rows
%	must be labelled
%	page index rows are labelled 'page' by default

	cprintf(s.nda,...
		'-mt',	1,...
		'-nd',	1);
%%
	cprintf(s.cmc,...
		'-mt',	1,...
		'-nd',	1);
%% ND structures

% a (ND) structure is displayed in the form
%
%	fieldname separator fieldcontent
%
% any content of a field is formatted according to the CELL formatting rules
% if the table option [-mt] is used, both fieldname(s) and fieldcontent have
%    their own column(s)
% the layout of ND structures follows the general ML rules for ND data types
%    and the [-nd] option may be used

% an ND structure with unique fieldcontent reflecting the indices

	s.siz=[2,2,3,2];
	s.fh=cprintf;
	s.ns=num2cell(s.fh.comb(s.siz));
	s.st=struct('a','','bb','','ccc','');
	s.st=repmat(s.st,s.siz);
for	ii=1:size(s.ns,1)
	s.txt=sprintf('%d.',s.ns{ii,:});
	s.txt(end)='';
	s.st(s.ns{ii,:}).a=sprintf('%s-A',s.txt);
	s.st(s.ns{ii,:}).bb=sprintf('%s-B',s.txt);
	s.st(s.ns{ii,:}).ccc='XYZ'.';			% a 3x1 character string
end
	clear ii;

	disp(s.st);
	disp('structure with indices (1,1,3,2)');
	disp(s.st(1,1,3,2));
%%

% NOTE:	display of Mx1 character string

	cprintf(s.st(1,1,3,2),...
		'-mt',	1,...
		'-Lc',	{'fieldname','value'},...
		'-dt',	' :  ');
%%

% NOTE:	page indices 1 and 2 display the number of structures and
%	NOT the number of fieldnames and columns

	cprintf(s.st,...
		'-mt',	1,...
		'-nd',	1,...				% show page indices
		'-nex',	1,...				% do not expand char strings
		'-dt',	'|');				% table cell separator
%%

% show content of a LINE graphics handle

	s.fh=figure('visible','off');
	cprintf(struct(handle(line)),...
		'-mt',	1,...
		'-dt',	' :  ');
	delete(s.fh);					% clean up
%% --- anomalous structures

% anomalous structures of the form

	s.san=struct('a',{},'b',1);

% are valid MATLAB constructs

	s.san
	clear ans;					% clean up

% but yield an error if used like
%{
	val=s.san.a
??? Too many outputs requested.  Most likely cause is missing [] around
left hand side that has a comma separated list expansion.
%}

% they are displayed with a fieldcontent of [?]

	cprintf(s.san,...
		'-mt',	1,...
		'-Lc',	{'fieldname','value'},...
		'-dt',	' :  ');
%% ND arrays and ND cell arrays of other objects

% other objects are displayed according to the [-O] options

if	ispc
	s.fh=figure('visible','off');
	s.com=actxcontrol('mscal.Calendar',[20 20 300 300],s.fh);
	s.fun=@(x) x+10;

	s.O=@(x) sprintf('%s(%s )',class(x),sprintf(' %d',size(x)));

	cprintf([s.com;s.com],...
		'-mt',	1,...
		'-O',	s.O);				% show class:size
%%
%
	s.com={s.com,s.fun;s.fun,s.com};

	cprintf(s.com,...				% using the built-in other engine
		'-mt',	1);
%%
%
	cprintf(s.com,...
		'-mt',	1,...
		'-O',	s.O);				% show class:size

	delete(s.fh);
end
%% output files

% results can be written or appended to ASCII files
% file templates can be used to insert results and replace token templates
%    with user defined values
%-------------------------------------------------------------------------------
% -fa		append to output file			aname	[]
% -fc		create output file			cname	[]
%		precedence: [-fa] > [-fc]
%
% -fi		template file				tname	[]
%		may contain markers as single line
%		entries where the result is
%		inserted
%		for each run, only the FIRST marker
%		is being used
%		if no marker is found and the
%		[-fa] option is used, the result is
%		appended to the content of the template
%		file before creating the output file
% -fm		marker for the [-fi] option		marker	[]
%		must be a 1xN character string
% -fr		cell of token/value pairs		{t,v...}
%		must be 1xN character strings
%		all tokens are replaced with their
%		values in a [-fi] template file
%		before insertion of the result

% NOTE:	output files contain the exact result including
%	TAB characters
%	markers and tokens are replaced using the
%	regular expression engine
%	therefor, if special characters are used,
%	they must be marked as literal by \char
%
% the accompanying file CPTMPL.TXT is used for demonstration
% a file CP_00X.00X is created during the demonstration and
% removed afterwards

	s.ftmpl='cptmpl.txt';
	s.fout='cp_00X.00X';

	type(s.ftmpl);
%% --- creating output files

% using a return argument prevents the result from being displayed after
%    runtime
% a print message is always printed to confirm the writing process

	s.p=cprintf(s.m,...
		'-mt',	1,...
		'-fc',	s.fout);
	type(s.fout);
%% --- appending to files

	s.p=cprintf(s.r,...
		'-mt',	1,...
		'-fa',	s.fout,...
		'-dt',	'|');
	type(s.fout);
	delete(s.fout);
%% --- using a template file

	type(s.ftmpl);
	s.p=cprintf(s.m,...
		'-mt',	1,...
		'-fc',	s.fout,...
		'-fi',	s.ftmpl,...
		'-fm',	'\$MARK\$');			% note liteal \$

% NOTE:	only marker #1 is replaced

	type(s.fout);
%%

% to replace other markers, the newly created FOUT must be used

	s.p=cprintf(s.r,...
		'-mt',	1,...
		'-dt',	'|',...
		'-fc',	s.fout,...
		'-fi',	s.fout,...
		'-fm',	'\$MARK\$');
	type(s.fout);
	delete(s.fout);
%% --- replacing tokens in template files

	s.tok={
		'FO',		'file'
		'<DATE>',	datestr(clock)
		'output \$',	'result #'
	};
	s.p=cprintf(s.r,...
		'-mt',	1,...
		'-dt',	'|',...
		'-fc',	s.fout,...
		'-fi',	s.ftmpl,...
		'-fr',	s.tok,...
		'-fm',	'\$MARK\$');
	s.p=cprintf(s.m,...
		'-mt',	1,...
		'-dt',	'|',...
		'-fa',	s.fout,...
		'-fi',	s.fout,...
		'-fr',	s.tok,...
		'-fm',	'\$MARK\$');

% NOTE:	no marker exists in S.FOUT during the second call
%	the result is appended       with the [-fa] option
%	whereas an error would occur with the [-fc] option

	type(s.fout);
	delete(s.fout);
%% option processing
%% --- option precedence

% options are resolved in this order of precedence from lowest to highest
%
%	1) getpref					low
%	2) option file
%	3) option structure
%	4) command line					high
%
% NOTE:	if multiple options are set, the LAST values is used
%% --- option preferences across MATLAB sessions

% user may define preferred options, which are kept across MATLAB sessions, by
%
%	setpref('cprintf','opt',{'OPT1',VAL1,...,'OPTn',VALn});

% setting [-mt] options by default

	s.opref=getpref('cprintf');			% save old state
	setpref('cprintf','opt',{'-mt',1});

	cprintf(s.r);
%%

% modification of preferred options
%-------------------------------------------------------------------------------
% -p		do NOT use preferred options		[F]

% NOTE:	[-p] requires NO value

	cprintf(s.r,...
		'-p',...				% do NOT use preferred options
		'-d',	' | ');
%%

% in addition, any runtime option will override a preferred option

	cprintf(s.r,...
		'-mt',	0,...				% override preferred [-mt]
		'-d',	' | ');

if	~isempty(s.opref)
	setpref('cprintf','opt',s.opref.opt);		% reset
else
	rmpref('cprintf');				% clean up
end
%% --- reading and writing option files

% options can be retrieved from any user created ASCII file
% there are two possible file formats:
% 1)		user created files, which contain a sequence of option/value pairs
%		as if typed on the command line
%		if it is an m-file, the extension is not necessary, otherwise, it
%		must be added to identify the file
% 2)		files generated by CPRINTF, which contain a sequence of commands
%		that create an option structure
%		these files are m-files with a short help header,
%		which may be used as stand-alones as well to
%		retrieve the option structure
%
% read and create option files
%-------------------------------------------------------------------------------
% -ini		read options from file			iname	[]
%		format: user created or CPRINTF
% -sav		save options to   file			sname	[]
%		which are created from a [-ini] file
%		if this option is not used, a temporary
%		ASCII file is created in the current
%		folder and discarded immediately after
%		loading
%
%		see also -macros- below for possibilities
%		to read/create option files without running
%		CPRINTF
%
% NOTE:	if the [-ini] option is used without
%	the [-sav] option, a temporary files is
%	created and removed right after options
%	are read if the [-ini] and [-sav] use
%	the same file name, a new copy of the
%	file is saved without warning
%	only command line options are saved
%%

% user created ASCII files
% - option/values must be written as if typed at the command line, eg,
%		'-mt',1,'-Lc',{1,2,'a'}
% - the parsing/decoding of these files starts at the first occurrence of
%   an option
% - option names must be convertible to valid structure fieldnames
% - values may contain any regular, valid ML expression and must be
%   separated by at least one <,> from the option
% - comments may be used according the ML's common rules
% - white spaces, <;>, and multiple <,>s are discarded
% - newlines are discarded

% CPINI.TXT is an exemplary [-ini] file showing some syntax possibilities,
% which typically are less disorganized and chaotic
% - parsing starts at '-begin'
% - note several ML type comments starting with a <%>
% - non-CPRINTF options may be used; they are created but discarded at runtime
% - any valid ML command in between option/value pairs will be created during
%   runtime and results may be displayed in the command window if it is not
%   suppressed by a <;>

	s.fini='cpini.txt';
	s.fout=sprintf('cpini_%s',num2hex(rand(1,2)));	% a unique file name
							% extension not required

	type(s.fini);
%%

% load CPINI.TXT during runtime and save it as a CPRINTF option file

% NOTE:	correct formatting of table elements according to options

	cprintf(bitand(s.r,1)~=0,...
		'-mt',	1,...
		'-ini',	s.fini,...
		'-sav',	s.fout);
%%

% output and standalone use of CPRINTF created m-file

% NOTE:	correct formatting of multiline commands
%	positioning of comments
%	the option structure contains all option/value pairs

	disp('***** FILE CONTENT');
	type(s.fout);
%%
	disp('***** HELP SECTION');
	help(s.fout);
%%
	disp('***** OPTION STRUCTURE');
	s.val=feval(s.fout);				% typically: s.val=foo;
	disp(s.val);
%% --- option structures

% a predefined option structure may be used to easily enter options
%-------------------------------------------------------------------------------
% -opt		an option structure			os	[]

% NOTE:	an option structure consists of
%
%		s.option = value
%
%	where option is a valid CPRINTF option
%	without the [-] and a valid value
%	unused options are discarded
%	CPRINTF created [-ini] files produce option structures

	s.sopt.a='a';					% currently not an option
	s.sopt.Lc={1,'XYZ',2};
	s.sopt.Lcs='=';
	s.sopt.Lr=s.rlab;
	
	cprintf(s.r,...
		'-mt',	1,...
		'-opt',	s.sopt);
%%

% using a CPRINTF created ini file, which returns an option structure

	cprintf(s.r>12,...
		'-mt',	1,...
		'-opt',	feval(s.fout));			% typically: '-opt','foo'

	delete([s.fout,'.m']);				% to clean up we need the extension
%% engine options
%% --- TAB size and TAB replacement in CELLs

% conversion of TAB [\t] characters to SPACES for format specs
%    OTHER than the CELL separator [-d]
%-------------------------------------------------------------------------------
% -tab		number of SPACES/TAB			n	[8]
% -ntab		do NOT replace TABs in CELLs		T|F	[F]

% compare with output above

	cprintf(s.m,...
		'-mt',	1,...
		'-Ct',	'{%s\t}',...
		'-tab',	2);
%%
%

% NOTE:	replacement of CELL TABS with TAB size SPACES in the output

	s.stab={sprintf('a\t\t\ttab'),sprintf('b\t\ttab'),'c d';100,200,300};
	cprintf(s.stab,...
		'-mt',	1);
%%
%

% NOTE:	with the [-ntab] option, table colums have an additional length
%	of one char <\t> per TAB and the TABs are kept in the output
%	this may be required in other applications

	s.res=cprintf(s.stab,...
		'-mt',	1,...
		'-ntab',1);

	s.rres=strrep(cellstr(s.res),sprintf('\t'),'.');
	disp(sprintf('\noriginal output\n'));
	disp(s.res);
	disp(sprintf('\noriginal output with later TAB replacements\n'));
	disp(char(s.rres));
%% --- runtime printout

% runtime processing and timing of conversion stages can be print out
%-------------------------------------------------------------------------------
% -db		show processing of conversions		T|F	[F]
% -q		quiet mode				T|F	[F]
%		except for [-db] output ALL
%		runtime processing for
%		[-fc]/[-fa]/[-ini] options
%		as well as the display of the result
%		is suppressed even if the command
%		does NOT end with a <;>

	cprintf(s.m,...
		'-mt',	1,...
		'-db',	1);

% explanation of runtime output
%-------------------------------------------------------------------------------
% - numeric col 1:	 # of cells remaining
%			-# of cells converted
%			if the array contains complex numbers, the -# will
%			   be larger than the #
% - numeric col 2:	time [sec] used to convert

% NOTE:	after formatting, 0 c should be left to convert
%% output
%% --- macros

% several macros (function handles) are available if the engine structure
%    is retrieved, which may be helpful to prepare and/or organize often
%    used options
% if no output argument is selected, macros do NOT return the default
%    ANS
%-------------------------------------------------------------------------------
% .getpar	to create an option structure from
%		command line parameters
% .setpar	to create a CPRINTF options only
%		option structure from a user defined
%		option structure
% .ini		to read/create ini files
% .write	to write an output file
% .other	to test the [-O] display engine, which
%		must produce an output for any data
%		type and array size
% .comb		to create a page index sequence

% macros can be retrieved by a simple call to CPRINTF without arguments,
%    returns a structure with the engine parameters; see also -arguments-

	s.mac=cprintf;
%%

% GETPAR	OS = P.getpar(OPT1,VAL1,...,OPTn,VALn)

% NOTE:	options currently not used by CPRINTF are discarded
%	values, which would not be active during runtime
%	are discarded as well, eg, [-nd]

	s.opt=s.mac.getpar('-foo','goo','-Lc',{'a',1,'b'},'-mt',0,'-nd',[]);
	disp(s.opt);
%%

% SETPAR	OS = P.setpar(OS)

% NOTE:	same as GETPAR
%	the input may be a CPRINTF created [-sav] file

%	s.opt=s.mac.setpar('CPRINTF_savfile');
	s.opt=s.mac.setpar(s.val);				% [-ini] created option structure
	disp(s.opt);
%%

% INI		OS = P.ini(INIFILE,[SAVFILE],['-q',1,...,OPTn,VALn]);

% NOTE:	see section -option files- for rules regarding the creation of output files
%	to add options, SAVFILE must be set
%	to create a temporary file, SAVFILE can be set to ''

	s.opt=s.mac.ini(s.fini);
	disp(s.opt);
%%

% WRITE		[ARG,...] = P.write(A,OPT1,VAL1,...,OPTn,VALn)

% NOTE:	mimicks the command
%
%		[arg,...] = cprintf(A,WOPT,WVAL,OPT1,VAL1,...,OPTn,VALn)
%
%	one of the WOPTs must be a [-fc]/[-fa] file creation command
%	the [-q] silent processing flag is turned on [T] by default
%	the [-mt] table option is turned off [F] by default and
%	must be set if this feature is required

	s.fwrite='CP_00X.00X';
	s.mac.write(s.r,'-fc',s.fwrite);
	type(s.fwrite);
	delete(s.fwrite);				% clean up
%%

% OTHER		OVAL = P.other(O1,...,On)

	s.oval={pi,magic(3)+1i,cell(2,2),sparse(2,2,3+1i),@(x) x+10};
	s.oval=s.mac.other(s.oval{:});
	disp(s.oval);
%%

% COMB		PIX = P.comb(ARRAYSIZE)

	s.pix=s.mac.comb(size(rand([2,2,1,2])));
	disp(s.pix);
%% --- arguments

% the engine structure is returned by a call to CPRINTF without input arguments
% this may be useful to access several CPRINTF functions
%
	s.p=cprintf;

%	s.p.other()	check other display output

% several outputs are available to import the result into other programs,
%    in particular, spreadsheets, which take cell input, eg,
%    XLSWRITE (see OUTPUT section)
% all output can be collected into a single structure
%-------------------------------------------------------------------------------
% -so		collect output				T|F	[F]

% NOTE:	the result is NOT displayed during runtime because of
%	1) output arguments and 2) the semicolon at the end

	[s.rs,s.rc,s.mc,s.par]=cprintf(s.m,...
		'-mt',	1);

% the result

	disp(s.rs);
%%

% the output prior to conversion

	disp(s.mc);
%%

% the output in a cell array with formatted CELLs

	disp(s.rc);
%%

% all output collected in a structure

	s.p=cprintf(s.m,...
		'-mt',	1,...
		'-so',	1);

% show engine parameters, input, and output

	disp(s.p);
%%

% same as RC above

	disp(s.p.cell);
%% --- input to xlswrite

% use RC as input to xlswrite
%{
	s.rc=strtrim(s.rc);	% if leading/trailing spaces need to be removed
	s.fnam='foo.xls';
	xlswrite(s.fnam,s.rc);
	winopen(s.fnam);
%}