function varargout = save4spss(varargin)

% save4spss saves Matlab numerical matrices as ASCII data and 
% creates a corresponding SPSS-syntax-file for easy import in SPSS
%
% save4spss([DATA],[VARNAMES],[FILE]);
%          DATA: Data matrix (optional)
%          VARNAMES (optional): cell array of strings
%          FILE (optional): Path and file name of output file as a string 
%
% save4spss()     opens dialogue boxes to select a matrix,
%                 a cell array with variable names (optional),
%                 and specify a file name to save the matrix
%                 the SPSS syntax gets the same name
% save4spss(DATA) saves the matrix DATA as 'spss.dat' in the current directory
%                 and creates a SPSS syntax file spssdat.sps to import these data
%                 in SPSS
%                 variables (data columns) are named 'V1' to 'Vn'                  
% save4spss(DATA,VARNAMES) saves the matrix DATA as 'spss.dat' in the current directory
%                 and creates a SPSS syntax file spssdat.sps to import these data
%                 in SPSS
%                 variables (data columns) are named according to the elements of 
%                 cell array VARNAMES 
% save4spss(DATA,VARNAMES,FILE) saves the matrix DATA in the file given in FILE
%                 and creates a SPSS syntax file FILE.sps to import these data
%                 in SPSS
%                 variables (data columns) are named according to the elements of 
%                 cell array VARNAMES 
%
% after running save4spss, just open the generated SPSS syntax in SPSS and run it

% example:
% mystudy=[1 30 10000; 2 35 15000];
% varnames={'subject';'age';'income'}
% save4spss(mystudy,varnames,'mystudy.dat')
% creates an ASCII file mystudy.dat
% and a SPSS syntax mystudy.sps
% if you run this syntax in SPSS
% the data are imported, the variables named
% and the file is saved as mystudy.sav
%***********************************

%written by r.schleicher (at uni-koeln.de)

% get user input via dialogue boxes
if nargin==0 & exist('uigetVariable.m','file') & exist('uigetVariable.fig','file') 
   data=uigetVariable('isnumeric','Get data matrix to export');
   if data==-1
       return
   end
   varnames=uigetVariable('iscell','Get variable names');
   if ~iscell(varnames) | (iscell(varnames) & ~ischar(varnames{1}))
      varnames={};
    
   end   
   
   [filename,pathname]=uiputfile({'*.dat';'*.txt'},'Save ASCII-Output for SPSS Import');
   if filename==0 %'Cancel'
       return
   end
   filename=[pathname filename];
   syntaxfile=[filename(1:end-3) 'sps'];
end
%or get user input from command line:
%input data matrix
if (nargin>=1 & isnumeric(varargin{1}))
    data=varargin{1};
    filename='spss.dat';
    syntaxfile='spssdat.sps';
    varnames={};
end
%input cell array with variable names
if (nargin>1 & iscell(varargin{2}))
    varnames=varargin{2};
end

%alternative name for syntax file
if (nargin>2 & isstr(varargin{3}))
    filename=varargin{3};
    syntaxfile=[filename(1:end-3) 'sps'];
end

disp(['Output file: ' filename]);
disp(['SPSS Syntax file: ' syntaxfile]);

%***********************************
%write ASCII-Output
formatstring=[];
fid=fopen(filename,'w');
for i=1:size(data,2) %write varnames for each column
    if i<=size(varnames,1) & isstr(varnames{i}) %use input varnames
        fprintf(fid,'%s\t',varnames{i});
    else %write default names
        fprintf(fid,'%s\t',['V' num2str(i)]);
    end
    formatstring=[formatstring ' %f'];
end
fprintf(fid,'\r\n'); 
formatstring=[formatstring '\r\n'];

%write data
for i=1:size(data,1)
    fprintf(fid,formatstring,data(i,:));
end
fclose(fid);

%***********************************
% write SPSS-Syntax
% Remark: this function was written for the German version of SPSS
% if you have trouble with the decimal point/comma, play around with the
% import formats of SPSS and change the formatting commands that save4spss()
% writes in the syntax file.

fid=fopen(syntaxfile,'w');
fprintf(fid,'GET DATA  /TYPE = TXT\r\n');
fprintf(fid,' /FILE = ''%s''\r\n ' ,filename);
fprintf(fid,' /DELCASE = LINE\r\n');
fprintf(fid,' /DELIMITERS = "\\t "\r\n');
fprintf(fid,' /ARRANGEMENT = DELIMITED\r\n');
fprintf(fid,' /FIRSTCASE = 2\r\n');
fprintf(fid,' /IMPORTCASE = ALL\r\n');
fprintf(fid,' /VARIABLES =\r\n');
for i=1:size(data,2) %write varnames for each column
    if (i<=size(varnames,1)& isstr(varnames{i})) %use input varnames
        fprintf(fid,' %s COMMA12.4\r\n',varnames{i});
    else %write default names
        fprintf(fid,' %s COMMA12.4\r\n',['V' num2str(i)]);
    end
end

fprintf(fid,' .\r\n');
fprintf(fid,'CACHE.\r\n');
fprintf(fid,'EXECUTE.\r\n');

%these lines reformat COMMA-Variables to 'normal' numerical values

for i=1:size(data,2) %write varnames for each column
    if (i<=size(varnames,1)& isstr(varnames{i})) %use input varnames
        fprintf(fid,'FORMATS %s (F12.4) .\r\n',varnames{i});
    else %write default names
        fprintf(fid,'FORMATS %s (F12.4) .\r\n',['V' num2str(i)]);
    end
end
fprintf(fid,'EXECUTE.\r\n');
fprintf(fid,'SAVE OUTFILE=''%s'' \r\n',[filename(1:end-3) 'sav']);
fprintf(fid,'/COMPRESSED.');

fclose(fid);
