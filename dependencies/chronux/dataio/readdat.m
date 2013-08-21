function [data,sampling]=readdat
%
% function for reading *.mat, *.wav and *.au files
% *.mat files have to contain variables data and sampling in order for the
% returns to work properly
%
[fname,pname]=uigetfile('*.*');
strindex=findstr(fname,'.plx');
if ~isempty(strindex);
   error('This file format can be read with Plexon supplied routines in folder ReadingPLXandDDTfilesinMatlab');
end;
strindex=findstr(fname,'.ddt');
if ~isempty(strindex);
   error('This file format can be read with Plexon supplied routines in folder ReadingPLXandDDTfilesinMatlab');
end;
strindex=findstr(fname,'.nex');
if ~isempty(strindex);
   error('This file format can be read with Plexon supplied routines in folder HowToReadNexFilesInMatlab');
end;
strindex=findstr(fname,'.mat');
if ~isempty(strindex);
   eval(['load ' pname fname]);
end;
strindex=findstr(fname,'.wav');
if ~isempty(strindex);
   filename=[pname fname];
   [data,sampling]=wavread(filename);
end;
strindex=findstr(fname,'.au');
if ~isempty(strindex);
   filename=[pname fname];
   [data,sampling]=auread(filename);
end;
