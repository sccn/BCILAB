function env_compile_bcilab
% Compile a deployable BCILAB binary using the MATLAB compiler toolbox
%
% The resulting binary is found in bcilab/build/distrib/. This binary
% can be executed just like the startup function "bcilab".
%
% Notes:
%   The compilation is done by calling the deploytool function with a .prj file. This file has a few quirks:
%   * contains absolute paths (i.e. won't compile when placed in another directory)
%   * contains a 1:1 copy of the entire MATLAB path (i.e. won't compile when toolboxes move)
%   * files that only appear in eval() or feval() clauses are not found and will be omitted (unless explicitly noted somewhere)
%   * files in 'private' sub-directories that are only referenced in eval() or feval() must be added to the project to be included
%   * there is OS-specific stuff in the file
%
%   To bypass these quirks, BCILAB uses the following (currently somewhat arduous) process:
%   When a new release of BCILAB is made (new dependencies added, another MATLAB release), 
%   then for each operating system, we do the following (currently using MATLAB 2012a):
%   1. start the deploy tool as follows:
%      * cd /path/to/bcilab; bcilab
%      * cd .. (move up to the directory that contains bcilab.m)
%      * deploytool
%   2. create a new project named 'build.prj', and:
%      * add the bcilab.m file as the *main file*
%      * add the code directory as an additional file/directory
%      * add the file build/dependency_list as an additional file/directory
%      * add the path dependencies/eeglab-*/plugins/dipfit2.2/private as an additional file/directory
%      * add the path dependencies/eeglab-*/external/fieldtrip_partial/forward/private as an additional file/directory
%      * add the path dependencies/eeglab-*/external/fileio/private as an additional file/directory
%      * in Project/Settings, Toolboxes on Path, remove all toolboxes that are not needed (lengthen the compilation process unnecessarily)
%        The ones that are currently included are: Bioinformatics Toolbox, Curve Fitting Toolbox,
%        Global Optimization Toolbox, Image Procesing Toolbox, Instrument Control Toolbox,
%        Optimization Toolbox, Parallel Computing Toolbox, Signal Processing Toollbox, Statistics
%        Toolbox, Wavelet Toolbox
%      * quit the deploy tool
%   3. open the build.prj file in a text editor, and:
%      * substitute all occurrences of /your/path/to/bcilab/ or C:\your\path\to\bcilab\ (note the trailing slash) by bcilab:/ (always use a forward slash in the substitution)
%      * also substiute the one occurrence of the path without slash in the third line (after location=") by bcilab:/
%      * If you are on Windows and you would like to have console output (highly recommended), you also need to replace on the 
%        same line the string target="target.standalone.win" by target="target.standalone".
%      * in MATLAB, type lower(computer) to get the identifier for your architecture
%      * Save the file, and rename it to build/build-<arch>.reference (where <arch> is your architecture identifier, e.g. pcwin64)
%   4. Run env_compile_bcilab.
%   
%   When you try to compile the result on some machine, the toolbox will grab the architecture-specific version, substitute the actual 
%   toolbox installation directory, rename it to build.prj, and run deploytool on it.
%
% See also:
%   env_buildslave, deploytool
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2011-03-13


% these are our additional dependencies (just in case they are not caught by the compiler)...
special_dependencies = { ...
    ... % part of PropertyGrid
    'resolve_function', ...
    ... % filter functions
    'buttord','cheb1ord','cheb2ord','ellipord', ...
    ... % statistics functions 
    'betastat','binostat','chi2stat','copulastat','evstat','expstat','fstat', 'gamstat','geostat','gevstat','gpstat','hypgestat', ...
    'lognstat','mnstat','mvnstat','mvtstat','nbinstat','ncfstat','nctstat','ncx2stat','normstat','poissstat','raylstat','tstat', ...
    'unidstat','unifstat','wblstat','betapdf','binopdf','chi2pdf','copulapdf','evpdf','exppdf','fpdf','gampdf','geopdf','gevpdf', ...
    'gppdf','hypgepdf','lognpdf','mnpdf','mvnpdf','mvtpdf','nbinpdf','ncfpdf','nctpdf','ncx2pdf','normpdf','poisspdf','raylpdf', ...
    'tpdf','unidpdf','unifpdf','wblpdf',...
    ... % optimization functions 
    'optimset', 'fminsearch', 'fminbnd', 'fzero', 'lsqnonneg', 'optimget', ...
    ... % stuff that shows up inside eval() or feval() somewhere...
    'logpfun','ffun','cvx_s_banded','cvx_s_complex','cvx_s_diagonal','cvx_s_hankel','cvx_s_hermitian','cvx_s_lower_bidiagonal', ...
    'cvx_s_lower_hessenberg','cvx_s_lower_triangular','cvx_s_scaled_identity','cvx_s_skew_symmetric','cvx_s_sparse','cvx_s_symmetric', ...
    'cvx_s_symmetric_ut','cvx_s_toeplitz','cvx_s_tridiagonal','cvx_s_upper_bidiagonal','cvx_s_upper_hankel','cvx_s_upper_hessenberg','cvx_s_upper_triangular', ...
    'eeg_leadfield1','eeg_leadfield4','print','crossf','spectopo','timef','fmrib_fastr','fmrib_pas','plsTN','plsCGBT','plsCG','plsLBFGS','plsSB','plsBB', ...
    'potCat','potExpPow','potGauss','potLaplace','potLogistic','potSech2','potT','penAbs','penAbsSmooth','penLogSmooth','penNegLin','penNegQuad','penPow','penPowSmooth', ...
    'penQuad','penVB','penVBNorm','penVBnormsplit','penZero','phi','rigidbody','globalrescale','traditional','nonlin1','nonlin2','nonlin3','nonlin4','nonlin5', ...
    ... % vdpgm functions
    'mkopts_avdp','mkopts_bj','mkopts_bjrnd','mkopts_cdp','mkopts_csb','mkopts_vb','mkopts_vdp', ...
    ... % CVX keywords
    'dual','epigraph','expression','expressions','hypograph','ln','maximise','maximize','minimise','minimize','subject','variable','variables', ...
    ... % some mex identifiers
    'pnet','libsvmread','libsvmwrite','lbfgs_','find_existing_source_transpose','find_sources_complement_grid_fast_int_c','mexkernel','mexsinglekernel','bsb_close','bsb_open','bsb_read', ...
    'lsl_loadlib','lsl_inlet','lsl_outlet','lsl_streaminfo','lsl_xml_element', ...
    ... % more stuff
    'eeg_load_xdf','load_xdf',...
    ...% EEGLAB plugins 
    'eegplugin_eepimport', 'eegplugin_bdfimport', 'eegplugin_brainmovie', 'eegplugin_bva_io', 'eegplugin_ctfimport', 'eegplugin_dipfit', ...
    'eegplugin_erpssimport', 'eegplugin_fmrib', 'eegplugin_iirfilt', 'eegplugin_ascinstep', 'eegplugin_loreta', 'eegplugin_miclust', 'eegplugin_4dneuroimaging',...
    'eegplugin_xdfimport'};

try
    % create a file where we list all the dependencies that are pulled in by the loader
    % (so that they are are found by it when compiled)
    listname = env_translatepath('bcilab:/build/dependency_list.m');
    io_mkdirs(listname,{'+w','a'});
    disp(['Creating dependency list ' listname]);
    fnew = fopen(listname,'w+');
    fprintf(fnew,'exist;\n'); % this will make it exit immediately
    % go through all dependency paths and add contents of all env_exec files...
    for p = hlp_split(genpath(env_translatepath('dependencies:/')),pathsep)
        fname = [p{1} filesep 'env_exec.m'];
        if exist(fname,'file')
            try
                % env_exec file: read contents and append to fnew...
                fid = fopen(fname,'r');
                contents = char(fread(fid))'; %#ok<FREAD>
                try
                    fwrite(fnew,[contents ';']);
                catch
                    disp_once('Error writing to the dependency list.');
                end
                fclose(fid);
            catch
                if fid ~= -1
                    try fclose(fid); catch, end
                end
                disp(['Error reading from the dependency loader ' fname]);
            end
        end
    end
    % also recursively add all functions in the code directory...
    allfiles = cellfun(@dir,hlp_split(genpath(env_translatepath('bcilab:/code')),pathsep),'UniformOutput',false);
    allfiles = vertcat(allfiles{:});
    filenames = {allfiles(~cellfun('isempty',strfind({allfiles.name},'.m'))).name};
    identifiers = cellfun(@(s)s(1:end-2),filenames,'UniformOutput',false);
    fprintf(fnew,'%s;\n',identifiers{:});
    fprintf(fnew,'%s;\n',special_dependencies{:});
    fclose(fnew);
catch 
    if fnew ~= -1
        try fclose(fnew); catch, end
    end
    disp('Error during creation of the dependency list.');
end

% make writable for all
try
    fileattrib(listname,'+w','a');    
catch
    disp('Error setting permissions for the dependency list.');
end

% and add to the path, so that it is properly referenced
if ~exist(listname,'file')
    error(['MATLAB did not find the dependency list ' listname]); end

% copy the build/build-<arch>.reference file to build.prj and substitute 
% bcilab:/ by the actual path to the toolbox
try
    srcfile = fopen(env_translatepath(['bcilab:/build/build-' lower(computer) '.reference']),'r');
    dstfile = fopen(env_translatepath('bcilab:/build.prj'),'w+');
    content = fread(srcfile, Inf, 'uint8=>char')';
    content = strrep(content,'bcilab:/',env_translatepath('bcilab:/'));
    fwrite(dstfile,content);
    fclose(srcfile);
    fclose(dstfile);
catch e
    disp('Encountered a problem trying to copy the build/build-<arch>.reference file to the build.prj file.');
    disp(['Reason: ' e.message]);
end

% now run the deploytool...
disp('Compiling BCILAB...');
clear functions;
deploytool('-build',env_translatepath('bcilab:/build.prj'));
disp('Done.');

% also try to build the documentation (note: this runs in parallel to the deployment process)
disp('Generating documentation...');
try
    cd(env_translatepath('bcilab:/'));
    m2html('mfiles','code','html','build/docs','recursive','on');
catch e
    disp('Error while generating documentation.');
    env_handleerror(e);
end
