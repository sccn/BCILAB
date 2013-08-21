% SB2_USEROPTIONS  User option specification for SPARSEBAYES
%
% OPTIONS = SB2_USEROPTIONS(parameter1, value1, parameter2, value2,...)
%
% OUTPUT ARGUMENTS:
% 
%	OPTIONS		An options structure to pass to SPARSEBAYES
% 
% INPUT ARGUMENTS:
% 
%	Optional number of parameter-value pairs to specify the following:
% 
%	ITERATIONS	Number of interations to run for.
% 
%	TIME		Time limit to run for, expressed as a space-separated 
%				string. e.g. '1.5 hours', '30 minutes', '1 second'.
% 
%	DIAGNOSTICLEVEL	Integer [0,4] or string to determine the verbosity of
%					diagnostic output.
%					0 or 'ZERO' or 'NONE'	No output
%					1 or 'LOW'				Low level of output
%					2 or 'MEDIUM'			etc...
%					3 or 'HIGH'
%					4 or 'ULTRA'
%
%	DIAGNOSTICFILE	Filename to write diagnostics to file instead of
%					the default stdout.
% 
%	MONITOR		Integer number: diagnostic information is output
%				every MONITOR iterations.
% 
%	FIXEDNOISE	True/false whether the Gaussian noise is to be fixed
%				(default: false.
%
%	FREEBASIS	Indices of basis vectors considered "free" and not 
%				constrained by the Bayesian prior (e.g. the "bias").
% 
%	CALLBACK	External function to call each iteration of the algorithm
%				(string). Intended to facilitate graphical demos etc.
% 
%	CALLBACKDATA	Arbitrary additional data to pass to the CALLBACK
%					function.
%
% EXAMPLE:
%
%	OPTIONS = SB2_UserOptions('diagnosticLevel','medium',...
%				  'monitor',25,...
%				  'diagnosticFile', 'logfile.txt');
%
% NOTES:
%
% Each option (field of OPTIONS) is given a default value in
% SB2_USEROPTIONS. Any supplied property-value pairs over-ride those
% defaults.
% 

%
% Copyright 2009, Vector Anomaly Ltd
%
% This file is part of the SPARSEBAYES library for Matlab (V2.0).
%
% SPARSEBAYES is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the Free
% Software Foundation; either version 2 of the License, or (at your option)
% any later version.
%
% SPARSEBAYES is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
% FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
% more details.
%
% You should have received a copy of the GNU General Public License along
% with SPARSEBAYES in the accompanying file "licence.txt"; if not, write to
% the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,
% MA 02110-1301 USA
%
% Contact the author: m a i l [at] m i k e t i p p i n g . c o m
%
function OPTIONS = SB2_UserOptions(varargin)

% Ensure arguments are supplied in pairs
% 
if rem(nargin,2)
  error('Arguments to SB2_UserOptions should be (property, value) pairs')
end
% Any options specified?
numSettings	= nargin/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set defaults
%
% Assume we will infer the noise in the Gaussian case
% 
OPTIONS.fixedNoise		= false;
%
% Option to allow subset of the basis (e.g. bias) to be unregularised
% 
OPTIONS.freeBasis		= [];
%
% Option to set max iterations to run for
%
OPTIONS.iterations		= 10000;
%
% Option to set max time to run for
%
OPTIONS.time			= 10000; % seconds
%
% Set options for monitoring and recording the algorithm's progress
% 
OPTIONS.monitor			= 0;
OPTIONS.diagnosticLevel	= 0;
OPTIONS.diagnosticFID	= 1; % stdout
OPTIONS.diagnosticFile_	= [];
%
% Option to call a function during each iteration (to create demos etc)
% 
OPTIONS.callback		= false;
OPTIONS.callbackFunc	= [];
OPTIONS.callbackData	= {};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parse string/variable pairs

for n=1:numSettings
  property_	= varargin{(n-1)*2+1};
  value		= varargin{(n-1)*2+2};
  switch upper(property_)
    %
   case 'FIXEDNOISE',
    OPTIONS.fixedNoise	= value;
    %
   case 'FREEBASIS',
    OPTIONS.freeBasis	= value;
    %
   case 'ITERATIONS',
    OPTIONS.iterations	= value;
	%
   case 'TIME',
    OPTIONS.time		= timeInSeconds(value);
    %
   case 'MONITOR',
    OPTIONS.monitor		= value;
    %
   case 'DIAGNOSTICLEVEL',
    OPTIONS.diagnosticLevel	= value;
    MAX_LEVEL	= 4;
    if ischar(value)
      switch upper(value)
       case {'ZERO','NONE'},
	OPTIONS.diagnosticLevel	= 0;
       case 'LOW',
	OPTIONS.diagnosticLevel	= 1;
       case 'MEDIUM',
	OPTIONS.diagnosticLevel	= 2;
       case 'HIGH',
	OPTIONS.diagnosticLevel	= 3;
       case 'ULTRA',
	OPTIONS.diagnosticLevel	= 4;
       otherwise,
	error('Unrecognised textual diagnostic level: ''%s''\n', value)
      end
    elseif isnumeric(value)
      if value>=0 & value<=MAX_LEVEL
	OPTIONS.diagnosticLevel		= value;  
      else
	error(['Supplied level should be integer in [0,%d], '...
	       'or one of ZERO/LOW/MEDIUM/HIGH/ULTRA'], MAX_LEVEL);
	
      end
    end
    %
   case 'DIAGNOSTICFILE',
    OPTIONS.diagnosticFile_	= value;  
    OPTIONS.diagnosticFID	= -1; % It will be opened later
    %
   case 'CALLBACK',
    OPTIONS.callback		= true;
    OPTIONS.callbackFunc	= value;
    if exist(OPTIONS.callbackFunc)~=2
      error('Callback function ''%s'' does not appear to exist\n', value)
    end
    %
   case 'CALLBACKDATA',
    OPTIONS.callbackData	= value;
   otherwise,
    error('Unrecognised user option: ''%s''\n', property_)
    %
  end
end

%%
%% Support function: parse time specification
%%
function s = timeInSeconds(value_)
%
[v_ r_]				= strtok(value_,' ');
v					= str2num(v_);
r_(isspace(r_))		= [];
switch upper(r_)
 case {'SECONDS', 'SECOND'}
  s	= v;
 case {'MINUTES', 'MINUTE'}
  s = 60*v;
 case {'HOURS', 'HOUR'}
  s = 3600*v;
 otherwise,
  error('Badly formed time string: ''%s''', value_)
end
% Returns time in seconds
% 