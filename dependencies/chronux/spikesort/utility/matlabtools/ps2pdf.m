function ps2pdf(filename)
%PS2PDF            Converts a PostScript file into Adobe PDF format.
%   PS2PDF(FILENAME), where FILENAME is a PostScript file, will convert a 
%   .ps file to an Adobe .pdf file.  Although Matlab supports directly
%   printing to a PDF driver, the PS driver has added flexibility (e.g.,
%   it allows pages to be appended to an existing file).  PS2PDF eases
%   the process of ending up with a PDF document after using the PS
%   flexibility while creating a document.
%
%   PS2PDF deletes the PostScript file after the PDF file is created.
%
%   This function uses the GhostScript driver packaged with Matlab
%   for portability.
%
%   See also REPORT_ADDPAGE.

% Argument checking.
if (~exist(filename, 'file')),  error('File does not exist.');  end;
[p,file,ext] = fileparts(filename);   
filename = fullfile(p, file);

% Standard locations for Matlab's GhostScript directories
s = filesep;
gs_root = [matlabroot s 'sys' s 'ghostscript'];
gs_bin  = [gs_root s 'bin' s 'win32' s 'gs'];
gs_init = [gs_root s 'ps_files'];
gs_font = [gs_root s 'fonts'];

% Set up a system call to GhostScript
libraries  = [' -I' gs_init ';' gs_font ' '];
nointeract =  ' -dBATCH -dNOPAUSE ';
pdfoutput  = [' -sDEVICE=pdfwrite -sOutputFile=' filename '.pdf ' ];
[s,w] = system([gs_bin libraries nointeract pdfoutput   filename '.ps']);

% Delete the .ps file
delete([filename '.ps']);
