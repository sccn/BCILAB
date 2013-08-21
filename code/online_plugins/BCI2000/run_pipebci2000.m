function run_pipebci2000(varargin)
% Integrate BCILAB into the BCI2000 real-time system.
% run_pipebci2000(BCI2000Directory)
%
% This function makes BCILAB usable as a Signal Processing module within the BCI2000 real-time
% experimentation environment. BCI2000 makes available a host of EEG acquisition systems, as well as
% stimulus presentation modules, many of are not natively supported by BCILAB. You can obtain BCI2000
% by following the "Getting Started" guide at http://www.bci2000.org/wiki/index.php/Main_Page.
%
% While BCI2000 ships several signal processing modules of its own (written in C++), it comes with
% an interface to MATLAB, allowing to rapidly develop custom signal processing implementations. This
% interface (implemented by the MatlabSignalProcessing.exe module of BCI2000) by default looks for a
% bci_Construct.m file in the directory /path/to/bci2000/prog, and runs it (as well as several
% others, if present). The function run_pipebci2000 installs a custom bci_Construct.m in this
% directory which refers the MATLAB runtime to BCILAB. This action can be easily reverted again by
% removing the file manually. Note that this integration does not interfere with the recommended way
% of using BCI2000 with MATLAB code, which involves passing a separate startup directory as a
% command line argument to MatlabSignalProcessing.exe.
%
% After the integration has been set up, classifiers can be trained using the standard BCILAB
% procedures and saved to disk for use with BCI2000. To use these classifiers inside the BCI2000
% environment, start BCI2000 either via a batch script (a template / example script is in
% /path/to/bci2000/batch/), or via the BCI2000 launcher GUI (there selecting the
% MatlabSignalProcessing module to run through BCILAB).
%
% Parameters of the processing system can then be configured in the BCI2000 operator GUI (like all
% other BCI2000 modules) under BCILAB - in particular, the classifier to be used can be selected. By
% default, the most recently saved classifier in /path/to/bcilab/resources/models would be used.
%
% In:
%	BCI2000Directory: a file or struct that contains a predictive model as previously computed by bci_train
%                     (default: 'lastmodel')
%
% Examples:
%   After having executed this function (also found in the GUI under Online Analysis / Process Data within...),
%   follow the steps in tutorial_bci2000.m to run a sample classifier in BCI2000.
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2011-06-02

% declare the name of this component (shown in the menu)
declare_properties('name','BCI2000 (integrate only)');

% read options
arg_define(varargin, ...
    arg({'bci2000path','BCI2000Directory'}, '', [], 'BCI2000 install directory. This is the location where bci2000 was installed;\n after clicking OK, BCILAB will be usable from within BCI2000 whenever the MatlabSignalProcessing module is used.','type','char'));

if ~exist(bci2000path,'dir')
    error('The chosen directory does not exist.'); end

if ~exist([bci2000path filesep 'prog'],'dir')
    error('This directory does not appear to contain a BCI2000 installation,'); end

fname = [bci2000path filesep 'prog' filesep 'bci_Construct.m'];

if exist(fname,'file')
    % a startup.m file already exists in this location: check contents...
    try
        fid = fopen(fname,'r');
        contents = fread(fid);
        if ~isempty(strfind(char(contents'),'bcilab'))
            % apparently, this startup file relates to BCILAB; assume that BCILAB has already been integrated
            x = questdlg2('BCILAB has apparently already been integrated into this BCI2000 installation; rebuild the integration? Note that this function only allows add BCILAB support to BCI2000; the actual online run must be started from within BCI2000 itself.','Note','Yes','No','Cancel','No');
            if any(strcmpi(x,{'no','cancel'}))
                return;  end
        else
            % apparently, this is another startup file for some other toolbox...
            x = questdlg2('It appears that this BCI2000 installation has already another MATLAB-based processing system integrated into it; replace by BCILAB?','Warning','Yes','No','Cancel','No');
            if any(strcmpi(x,{'no','cancel'}))
                return;  end
        end
    catch
        warndlg2('This BCI2000 installation apparently contains a startup.m file from some MATLAB-based processing system, but the file cannot be read (possibly due to permissions).');
        try close(fid); catch,end
        return;
    end

    % rename it
    [p,n,x] = fileparts(fname);
    backups = dir([p filesep n '-backup*' x]);
    backupname = [p filesep n '-backup' num2str(length(backups)+1) x];
    try
        movefile(fname,backupname);
        disp(['Renamed the previous startup.m file to ' backupname]);
    catch
        errordlg2(['Failed to rename the existing file ' fname ' to ' backupname ', possibly due to permissions.']);
        return;
    end
    errordlg2('Successfully integrated BCILAB into BCI2000.','Note');
end

% now create a new startup.m
try
    global tracking; %#ok<TLEV>
    fid = fopen(fname,'w+');
    fprintf(fid,'function [parameters,states] = bci_Construct()\n');
    fprintf(fid,'global tracking;\n');
    fprintf(fid,'if isempty(tracking)\n');
    fprintf(fid,'    cd %s; bcilab %s; end\n',env_translatepath('bcilab:/'),tracking.configscript);
    fprintf(fid,'[parameters,states] = bci_ConstructReal();\n');
    fclose(fid);
catch
    try fclose(fid); catch,end
    errordlg2(['Could not create the file ' fname ', possibly due to permissions.']);
end
