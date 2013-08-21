% runamica11() - Perform AMICA -- adaptive ICA using multiple models with shared components
%              Schedules a run on the cluster using num_procs (keyword, default=8) processors.
% Usage:
%         >> runamica11(file,outdir,chans,frames,'Key1',Value1,...);
%
% Inputs:
%
%   file        file name (including path) of file containing floating point data 
%   outdir      name of directory to write output (does not have to exist)
%   chans       number of channels in the data file
%   frames      number of frames in the data file
%
% Optional keywords:
%
%   indir               optional input directory from which to load init
%   num_models          number of models to learn, default = 1
%   num_mix_comps       number of mixture components in source model, def=3
%   share_comps         flag to share components when num_models > 1, def=0
%   comp_thresh         correlation threshold to share component, def=0.99
%   share_start         iteration to start looking for shared components, def=100
%   share_int           number of iterations between sharing checks, def=100
%   numprocs            number or processors (slots) to use, def=8
%   max_iter            maximum number of iterations to perform, def=2500
%   lrate               initial learning rate for natural gradient, def=0.1
%   lratefact           multiplicative factor by which to decrease lrate, def=0.5
%   minlrate            lrate after which to stop, def=1e-8
%   rholrate            initial lrate for shape parameters, def=0.05
%   rho0                initial shape parameter value, def=1.5
%   minrho              minimum shape parameter value, def=1.0
%   maxrho              maximum shape parameter value, def=2.0
%   rholratefact        multiplicative factor by which to dec rholrate, def=0.5
%   do_newton           flag for newton method, default = 1 (do newton)
%   newt_start          for newton method, iter at which to start newton, def=50
%   newtrate            for newton method, lrate for newton iterations, def=1.0
%   newt_ramp           for newton method, number of iter to ramp up to newtrate, def=10
%   writestep           iteration interval between output writes, def=10
%   write_nd            flag to write history of component update norms, def=1
%   write_llt           flag to write model log likelihoods of time points, def=1
%   do_reject           flag for doing rejection of time points, def=0
%   numrej              for rejection, number of rejections to perform, def=3
%   rejsig              for rejection, number of standard dev of likelihood
%                           below which to reject data
%   rejstart            for rejection, iteration at which to start reject, def=3
%   rejint              for rejection, iteration interval between reject, def=3
%   max_threads         maximum number of threads to use on a node, def=4
%   decwindow           moving average window to detect likelihood decrease, def=1
%   update_A            flag to update mixing matrices, def=1
%   update_c            flag to update model centers, def=1
%   update_gamma        flag to update model probabilities, def=1
%   update_alpha        flag to update source mixture proportions, def=1
%   update_mu           flag to update source mixture mixture locations, def=1
%   update_sbeta        flag to update source mixture scales, def=1
%   invsigmax           maximum value of inverse scale parameters, def=100.0
%   invsigmin           minimum value of inverse scale parameters, def=1e-8
%   do_rho              flag to update shape parameters, def=1
%   load_rej            flag to load LLt to get rejections from, def=0
%   load_param          flag to load parameters, def=0
%   pcakeep             for PCA reduction, number of components to keep, def=chans
%   doscaling           flag to rescale unmixing matrix rows to unit norm, def=1
%   scalestep           iteration interval at which to rescale unmixing rows, def=1
%   qsub                ['on'|'off'] use qsub to run in parallel (default) or run on 
%                       local machine
%   use_queue           name of the Grid Engine queue to use, def='all.q'
%
% Disabled:
%   kurt_start          for ext. infomax, iter to start kurtosis calc, def=3
%   num_kurt            for ext. infomax, number of kurtosis calc, def=5
%   kurt_int            for ext. infomax, iteration interval between calc, def=1
%   load_comp_list      flag to load component assignment list, def=0
%   do_mean             flag to remove mean from data, def=1
%   do_sphere           flag to sphere data before ica, def=1
%   doPCA               flag to to PCA dimensionalit reduction, def=0
%
% Outputs:
%   
%   To load output use the function loadmodout() after job ends:
%                         
%       mod = loadmodout11(outdir);
%
%   mod is a structure containing the output components and density models. mod.A(:,:,h) is the components for model h. 
%   mod.varord(:,h) is the index order of the components in variance order, mod.Lht is the likelihood of time
%   points for each model (if set), mod.LL is the history of the log likelihood over iterations, mod.c(:,h)
%   is the center for model h, mod.W(:,:,h) is the unmixing matrix for model h, and mod.S is the sphering matrix.
%                       
% See also: loadmodout11()
%
%

function modres = runamica11(file,outdir,chans,frames,varargin)

if nargin < 2
    if ~isnumeric(file)
        help runamica11;
        return;
    else
        chans = size(file,1);
        frames = size(file,2);
    end
else
    if ~isnumeric(chans)
        varargin = [{chans frames} varargin]; end
    chans = size(file,1);
    frames = size(file,2);
end

modres = [];

%%%%%%%%%%%%%%%%%%%%%%% Set up keyword default values %%%%%%%%%%%%%%%%%%%%%%%%%


basepath = [fileparts(mfilename('fullpath')) filesep]; % directory where this .m file is located
if exist([basepath 'amica11'],'file')
    AMDIR = basepath;
else
    AMDIR = '/data/common/amica/';
end
if exist([basepath 'mpich-local'],'dir')
    MPI_BIN = [basepath 'mpich-local/bin/'];
else
    MPI_BIN = '/home/jason/mpich2-1.3.2-install/bin/';
end
OMPI_BIN = '/opt/openmpi/bin/';

if ~system('qconf -sp mpich > /dev/null')
    use_pe = 'mpich';
    disp('Using mpich parallel environment.');
    RUN_LOCALLY = 0;
elseif ~system('qconf -sp orte > /dev/null')
    use_pe = 'orte';
    disp('Using orte parallel environment.');
    RUN_LOCALLY = 0;
elseif ~system('qconf -sp make > /dev/null')
    use_pe = 'make';
    disp('Using make parallel environment.');
    RUN_LOCALLY = 0;
else
    disp('No recognized parallel environment found. Run qconf -spl to get a list of available environments and use keyword use_pe.');
    disp('Running locally with maximum of 4 threads.');
    RUN_LOCALLY = 1;
    max_threads = 4;
end

qsub = 'on';
use_queue = 'all.q';
numprocs = 8;
fix_init = 0;
block_size = 128;
num_models = 1;
num_mix_comps = 3;
pdftype = 0;
max_iter = 2500;
share_comps = 0;
comp_thresh = 0.98;
share_start = 100;
share_int = 100;
lrate = 0.1;
minlrate = 1e-8;
lratefact = 0.5;
rholrate = 0.05;
rho0 = 1.5;
minrho = 1.0;
maxrho = 2.0;
rholratefact = 0.5;
kurt_start = 3;
num_kurt = 5;
kurt_int = 1;
do_newton = 1;
newt_start = 50;
newtrate = 1.0;
newt_ramp = 10;
do_reject = 0;
numrej = 3;
rejsig = 3.0;
rejstart = 2;
rejint = 3;
if RUN_LOCALLY == 1
    max_threads = 4;
else
    max_threads = 999;
end
writestep = 10;
write_nd = 1;
write_LLt = 1;
decwindow = 1;
maxdecs = 3;
update_A = 1;
update_c = 1;
update_gamma = 1;
update_alpha = 1;
update_mu = 1;
update_sbeta = 1;
invsigmax = 100.0;
invsigmin = 1e-8;
do_rho = 1;
load_comp_list = 0;
load_rej = 0;
load_W = 0;
load_c = 0;
load_gm = 0;
load_alpha = 0;
load_mu = 0;
load_sbeta = 0;
load_rho = 0;
load_param = 0;
do_mean = 1;
do_sphere = 1;
doPCA = 1;
pcakeep = chans;
pcadb = 30.0;
bytesize = 4;
doscaling = 1;
scalestep = 1;

%
%%%%%%%%%% Collect keywords and values from argument list %%%%%%%%%%%%%%%
%
       
   for i = 1:2:length(varargin) % for each Keyword
       Keyword = varargin{i};
       Value = varargin{i+1};

      if ~isstr(Keyword)
         fprintf('runamica(): keywords must be strings')
         return
      end
      Keyword = lower(Keyword); % convert upper or mixed case to lower

      if strcmp(Keyword,'num_models')
         if isstr(Value)
            fprintf('runamica(): num_models must be a positive integer');
            return
         else
            num_models = Value;
         end
      elseif strcmp(Keyword,'indir')
         if ~isstr(Value)
            fprintf('runamica(): indir must be a string');
            return
         elseif ~isempty(Value)
            indir = Value;
            load_param = 1;
            load_W = 1;
            load_c = 1;
            load_gm = 1;
            load_alpha = 1;
            load_mu =1;
            load_sbeta = 1;
            load_rho = 1;
            load_comp_list = 1;
         end
      elseif strcmp(Keyword,'num_mix_comps')
         if isstr(Value)
            fprintf('runamica(): num_mix_comps must be a positive integer');
            return
         else
            num_mix_comps = Value;
         end
      elseif strcmp(Keyword,'pdftype') 
         if isstr(Value)
            fprintf('runamica(): pdftype should be either 0, 1, 2, or 3');
            return
         else
            pdftype = Value;
         end
      elseif strcmp(Keyword,'outdir') 
         if ~isstr(Value)
            fprintf('runamica(): outputdir must be a string');
            return
         else 
            outdir = Value;
         end
      elseif strcmp(Keyword,'max_iter')
         if isstr(Value)
            fprintf('runamica(): max_iter must be a positive integer');
            return
         else
            max_iter = Value;
         end
      elseif strcmp(Keyword,'lrate')
         if isstr(Value)
            fprintf('runamica(): lrate must be a positive number');
            return
         else
            lrate = Value;
         end
      elseif strcmp(Keyword,'minlrate')
         if isstr(Value)
            fprintf('runamica(): minlrate should be lrate to stop at');
            return
         else
            minlrate = Value;
         end
      elseif strcmp(Keyword,'lratefact')
         if isstr(Value)
            fprintf('runamica(): lratefact must be a positive number');
            return
         else
            lratefact = Value;
         end
      elseif strcmp(Keyword,'rholrate')
         if isstr(Value)             
            fprintf('runamica(): rholrate must be a number');
            return
         else
            rholrate = Value;
         end
      elseif strcmp(Keyword,'rho0')
         if isstr(Value)
            fprintf('runamica(): rho0 must be a number');
            return
         else
            rho0 = Value;
         end
      elseif strcmp(Keyword,'minrho')
         if isstr(Value)
            fprintf('runamica(): minrho must be a number');
            return
         else
            minrho = Value;
         end
      elseif strcmp(Keyword,'maxrho')
         if isstr(Value)
            fprintf('runamica(): maxrho must be a number');
            return
         else
            maxrho = Value;
         end
      elseif strcmp(Keyword,'rholratefact')
         if isstr(Value)
            fprintf('runamica(): rholratefact must be a number');
         else 
            rholratefact = Value;
         end
      elseif strcmp(Keyword,'invsigmax')
         if isstr(Value)
            fprintf('runamica(): invsigmax must be a number');
         else 
            invsigmax = Value;
         end
      elseif strcmp(Keyword,'invsigmin')
         if isstr(Value)
            fprintf('runamica(): invsigmin must be a number');
         else 
            invsigmin = Value;
         end
      elseif strcmp(Keyword,'fix_init')
         if isstr(Value)
            fprintf('runamica(): fix_init should be 0 or 1');
            return
         else             
            fix_init = Value;
         end
      elseif strcmp(Keyword,'do_newton')
         if isstr(Value)
            fprintf('runamica(): do_newton should be 0 or 1');
            return
         else             
            do_newton = Value;
         end
      elseif strcmp(Keyword,'update_a')
         if isstr(Value)
            fprintf('runamica(): update_a should be 0 or 1');
            return
         else             
            update_A = Value;
         end
      elseif strcmp(Keyword,'update_c')
         if isstr(Value)
            fprintf('runamica(): update_c should be 0 or 1');
            return
         else             
            update_c = Value;
         end
      elseif strcmp(Keyword,'update_gamma')
         if isstr(Value)
            fprintf('runamica(): update_gamma should be 0 or 1');
            return
         else             
            update_gamma = Value;
         end
      elseif strcmp(Keyword,'update_alpha')
         if isstr(Value)
            fprintf('runamica(): update_alpha should be 0 or 1');
            return
         else             
            update_alpha = Value;
         end
      elseif strcmp(Keyword,'update_mu')
         if isstr(Value)
            fprintf('runamica(): update_mu should be 0 or 1');
            return
         else             
            update_mu = Value;
         end
      elseif strcmp(Keyword,'update_sbeta')
         if isstr(Value)
            fprintf('runamica(): update_sbeta should be 0 or 1');
            return
         else             
            update_sbeta = Value;
         end
      elseif strcmp(Keyword,'do_rho')
         if isstr(Value)
            fprintf('runamica(): do_rho should be 0 or 1');
            return
         else             
            do_rho = Value;
         end
      elseif strcmp(Keyword,'load_param')
         if isstr(Value)
            fprintf('runamica(): load_param should be 0 or 1');
            return
         else
            if Value == 1
                load_W = 1;
                load_c = 1;
                load_gm = 1;
                load_alpha = 1;
                load_mu =1;
                load_sbeta = 1;
                load_rho = 1;  
                load_comp_list = 1;
                load_param = 1;
            else
                load_W = 0;
                load_c = 0;
                load_gm = 0;
                load_alpha = 0;
                load_mu =0;
                load_sbeta = 0;
                load_rho = 0;  
                load_comp_list = 0;
                load_param = 0;
            end
         end         
      elseif strcmp(Keyword,'load_rej')
         if isstr(Value)
            fprintf('runamica(): load_rej should be 0 or 1');
            return
         else
            load_rej = Value;
         end         
      elseif strcmp(Keyword,'write_nd')
         if isstr(Value)
            fprintf('runamica(): write_nd should be 0 or 1');
            return
         else             
            write_nd = Value;
         end
      elseif strcmp(Keyword,'write_llt')
         if isstr(Value)
            fprintf('runamica(): write_llt should be 0 or 1');
            return
         else             
            write_LLt = Value;
         end
      elseif strcmp(Keyword,'newtrate')
         if isstr(Value)
            fprintf('runamica(): newtrate must be a number');
            return
         else
            newtrate = Value;
         end
      elseif strcmp(Keyword,'newt_start')
         if isstr(Value)
            fprintf('runamica(): newt_start must be a number');
            return
         else
            newt_start = Value;
         end
      elseif strcmp(Keyword,'newt_ramp')
         if isstr(Value)
            fprintf('runamica(): newt_ramp must be a number');
            return
         else
            newt_ramp = Value;
         end
      elseif strcmp(Keyword,'do_reject')
         if isstr(Value)
            fprintf('runamica(): do_reject must be 0 or 1');
            return
         else
            do_reject = Value;
         end
      elseif strcmp(Keyword,'numrej')
         if isstr(Value)
            fprintf('runamica(): numrej must be a number');
            return
         else
            numrej = Value;
         end
      elseif strcmp(Keyword,'qsub')
         if ~isstr(Value)
            fprintf('runamica(): qsub must be ''on'' or ''off''');
            return
         else
            qsub = Value;
         end
      elseif strcmp(Keyword,'rejstart')
         if isstr(Value)
            fprintf('runamica(): rejstart must be a number');
            return
         else
            rejstart = Value;
         end
      elseif strcmp(Keyword,'rejsig')
         if isstr(Value)
            fprintf('runamica(): rejsig must be a number');
            return
         else
            rejsig = Value;
         end
      elseif strcmp(Keyword,'rejint')
         if isstr(Value)
            fprintf('runamica(): rejint must be a number');
            return
         else
            rejint = Value;
         end
         
      elseif strcmp(Keyword,'use_queue')
         if ~isstr(Value)
            fprintf('runamica(): use_queue argument must be a queue name (string)');
            return
         else
            use_queue = Value;
         end
         
      elseif strcmp(Keyword,'use_pe')
         if ~isstr(Value)
            fprintf('runamica(): use_pe argument must be a parallel environment name (string)');
            return
         elseif ~isempty(Value)
            use_pe = Value;
         end
         
      elseif strcmp(Keyword,'share_comps')
         if isstr(Value)
            fprintf('runamica(): share_comps must be 0 or 1');
            return
         else
            share_comps = Value;
         end
      elseif strcmp(Keyword,'share_start')
         if isstr(Value)
            fprintf('runamica(): share_start must be a number');
            return
         else
            share_start = Value;
         end
      elseif strcmp(Keyword,'comp_thresh')
         if isstr(Value)
            fprintf('runamica(): comp_thresh must be a number');
            return
         else
            comp_thresh = Value;
         end
      elseif strcmp(Keyword,'share_int')
         if isstr(Value)
            fprintf('runamica(): share_int int must be a number');
            return
         else
            share_int = Value;
         end                  
      elseif strcmp(Keyword,'max_threads')
         if isstr(Value)
            fprintf('runamica(): max_threads must be a number');
            return
         else
             if RUN_LOCALLY == 0
                 max_threads = Value;
             elseif RUN_LOCALLY == 1
                 max_threads = min(max_threads,Value);
             end
         end
      elseif strcmp(Keyword,'writestep')
         if isstr(Value)
            fprintf('runamica(): writestep must be a number');
            return
         else
            writestep = Value;
         end         
      elseif strcmp(Keyword,'decwindow')
         if isstr(Value)
            fprintf('runamica(): decwindow must be a number');
            return
         else
            decwindow = Value;
         end         
      elseif strcmp(Keyword,'maxdecs')
         if isstr(Value)
            fprintf('runamica(): maxdecs must be a number');
            return
         else
            maxdecs = Value;
         end         
      elseif strcmp(Keyword,'numprocs')
         if isstr(Value)
            fprintf('runamica(): numprocs must be a number');
            return
         else
            numprocs = Value;
         end         
      elseif strcmp(Keyword,'pcakeep')
         if isstr(Value)
            fprintf('runamica(): pcakeep must be a number');
            return
         else
            pcakeep = Value;
         end         
      elseif strcmp(Keyword,'bytesize')
         if isstr(Value)
            fprintf('runamica(): bytesize must be a number');
            return
         else
            bytesize = Value;
         end         
      elseif strcmp(Keyword,'doscaling')
         if isstr(Value)
            fprintf('runamica(): doscaling must be a number');
            return
         else
            doscaling = Value;
         end         
      else
         fprintf(['runamica(): unknown flag: ' Keyword])
         return
      end
   end
%
%%%%%%%%%%%%%%%%%%%%% Check keyword values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%


%%%%%%%%%%%%%%%%%%%%% create the param file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if outdir(end) ~= filesep
    outdir(end+1) = filesep;
end
system(['mkdir ' outdir]);

if isnumeric(file)
    filename = ['tmpdata.fdt' ];
    data2write = file;
    file = fullfile(outdir,filename);
    chans = size(data2write,1);
    frames = size(data2write,2);
else 
    data2write = [];
end;             

% write data on disk
if ~isempty(data2write)
    disp(['Writing data file: ' fullfile(outdir,filename)]);
    floatwrite(data2write, fullfile(outdir, filename));
end;             


fid = fopen([outdir 'input.param'],'w');
if fid < 1
    errordlg('cannot create file in outdir','Bad Input','modal');
    return;
end    
fprintf(fid,'files %s\n',file);
fprintf(fid,'outdir %s\n',outdir);
fprintf(fid,'block_size %d\n',block_size);
fprintf(fid,'num_models %d\n',num_models);
fprintf(fid,'num_mix_comps %d\n',num_mix_comps);
fprintf(fid,'pdftype %d\n',pdftype);
fprintf(fid,'max_iter %d\n',max_iter);
fprintf(fid,'num_samples 1\n');
fprintf(fid,'data_dim %d\n',chans);
fprintf(fid,'field_dim %d\n',frames);
fprintf(fid,'field_blocksize 1\n');
fprintf(fid,'share_comps %d\n',share_comps);
fprintf(fid,'share_start %d\n',share_start);
fprintf(fid,'comp_thresh %f\n',comp_thresh);
fprintf(fid,'share_iter %d\n',share_int);
fprintf(fid,'lrate %f\n', lrate);
fprintf(fid,'minlrate %e\n', minlrate);
fprintf(fid,'lratefact %f\n', lratefact);
fprintf(fid,'rholrate %f\n', rholrate);
fprintf(fid,'rho0 %f\n', rho0);
fprintf(fid,'minrho %f\n', minrho);
fprintf(fid,'maxrho %f\n', maxrho);
fprintf(fid,'rholratefact %f\n',rholratefact);
fprintf(fid,'kurt_start %d\n',kurt_start);
fprintf(fid,'num_kurt %d\n',num_kurt);
fprintf(fid,'kurt_int %d\n',kurt_int);
fprintf(fid,'do_newton %d\n',do_newton);
fprintf(fid,'newt_start %d\n',newt_start);
fprintf(fid,'newt_ramp %d\n',newt_ramp);
fprintf(fid,'newtrate %f\n', newtrate);
fprintf(fid,'do_reject %d\n',do_reject);
fprintf(fid,'numrej %d\n',numrej);
fprintf(fid,'rejsig %f\n',rejsig);
fprintf(fid,'rejstart %d\n',rejstart);
fprintf(fid,'rejint %d\n',rejint);
fprintf(fid,'max_threads %d\n',max_threads);
fprintf(fid,'writestep %d\n',writestep);
fprintf(fid,'write_nd %d\n',write_nd);
fprintf(fid,'write_LLt %d\n',write_LLt);
fprintf(fid,'decwindow %d\n',decwindow);
fprintf(fid,'max_decs %d\n',maxdecs);
fprintf(fid,'fix_init %d\n',fix_init);
fprintf(fid,'update_A %d\n',update_A);
fprintf(fid,'update_c %d\n',update_c);
fprintf(fid,'update_gm %d\n',update_gamma);
fprintf(fid,'update_alpha %d\n',update_alpha);
fprintf(fid,'update_mu %d\n',update_mu);
fprintf(fid,'update_beta %d\n',update_sbeta);
fprintf(fid,'invsigmax %f\n',invsigmax);
fprintf(fid,'invsigmin %f\n',invsigmin);
fprintf(fid,'do_rho %d\n',do_rho);
fprintf(fid,'load_rej %d\n',load_rej);
fprintf(fid,'load_W %d\n',load_W);
fprintf(fid,'load_c %d\n',load_c);
fprintf(fid,'load_gm %d\n',load_gm);
fprintf(fid,'load_alpha %d\n',load_alpha);
fprintf(fid,'load_mu %d\n',load_mu);
fprintf(fid,'load_beta %d\n',load_sbeta);
fprintf(fid,'load_rho %d\n',load_rho);
fprintf(fid,'load_comp_list %d\n',load_comp_list);
if load_param == 1
    fprintf(fid,'indir %s\n',indir);
end
fprintf(fid,'do_mean %d\n',do_mean);
fprintf(fid,'do_sphere %d\n',do_sphere);
fprintf(fid,'doPCA %d\n',doPCA);
fprintf(fid,'pcakeep %d\n',pcakeep);
fprintf(fid,'pcadb %f\n',pcadb);
fprintf(fid,'byte_size %d\n',bytesize);
fprintf(fid,'doscaling %d\n',doscaling);
fprintf(fid,'scalestep %d\n',scalestep);
fclose(fid);

if RUN_LOCALLY == 1 || strcmp(qsub,'off')
    system([AMDIR 'amica11 ' outdir 'input.param']);
    modres = loadmodout11(outdir);
else    
    % create the qsub file
    fid = fopen([outdir 'qsub.sh'],'w');
    
    fprintf(fid,'#!/bin/bash\n');
    fprintf(fid,'#$ -cwd\n');
    if use_queue
        fprintf(fid,['#$ -q ' use_queue '\n']);
    end
    fprintf(fid,'#$ -j y\n');
    fprintf(fid,'#$ -S /bin/bash\n');
    
    fprintf(fid,['#$ -pe ' use_pe ' %d\n'],numprocs);
    
    if strcmp(use_pe,'mpich') || strcmp(use_pe,'mpich2')
        fprintf(fid,...
            [MPI_BIN 'mpirun -np $NSLOTS -machinefile $TMP/machines ' AMDIR 'amica11 ' ...
            outdir 'input.param']);
    else
        fprintf(fid,...
            [OMPI_BIN 'mpirun ' AMDIR 'amica11_ompi ' outdir 'input.param']);
    end
    fclose(fid);
    
    system(['ssh `qconf -ss` -n qsub ' outdir 'qsub.sh > ' outdir 'lastqsubid']);
    
    % get the qsub job id
    fid = fopen([outdir 'lastqsubid'],'r');
    str = fscanf(fid,'%s');
    fclose(fid);
    qsubid = sscanf(str(8:end),'%d');
    
    disp('Please disregard any Warnings about X11 forwarding.');disp(sprintf(''));
    disp(['qsub id = ' int2str(qsubid)]);
    
    
    disp('You can watch the last 25 lines of the output with:');
    disp(['   $ watch tail -n 25 qsub.sh.o' int2str(qsubid)]);
    disp(['or in Matlab:']);
    disp(['   >> system(''watch tail -n 25 qsub.sh.o' int2str(qsubid) ''')']);
    disp(['You can load the output as it''s running by typing:']);
    disp(['   >> mod = loadmodout11(''' outdir ''')']);
end
