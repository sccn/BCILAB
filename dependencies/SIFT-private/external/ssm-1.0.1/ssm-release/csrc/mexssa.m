function mexssa(varargin)

%%%% LAPACK on Unix %%%%
if isunix
    options = [varargin {'-Ddgemm#dgemm_' ...
                         '-Ddgemv#dgemv_' ...
                         '-Ddgetrf#dgetrf_' ...
                         '-Ddgetri#dgetri_' ...
                         '-Ddsytrf#dsytrf_' ...
                         '-Ddsytri#dsytri_' ...
                         '-Ddpotrf#dpotrf_' ...
                         '-Ddpotri#dpotri_' ...
                         '-Ddsyev#dsyev_' ...
                         '-Ddgeev#dgeev_' ...
                         '-Ddgesv#dgesv_' ...
                         '-Ddgels#dgels_'}];
else options = varargin; end

%%%% Compile core source files %%%%
mex(options{:}, '-c', 'ssa_matlab.c', 'ssa_kalman.c', 'ssa_smooth.c', 'ssa_misc.c', 'mt19937ar.c');

%%%% Output directory %%%%
options     = [options {'-outdir' '../@ssmodel/private'}];

%%%% Object files and libraries %%%%
if ispc
    if verLessThan('matlab', '7.5')
        objfiles = {'ssa_matlab.obj' 'ssa_kalman.obj' 'ssa_smooth.obj' 'ssa_misc.obj' 'mt19937ar.obj' 'libmwlapack.lib'};
    else
        objfiles = {'ssa_matlab.obj' 'ssa_kalman.obj' 'ssa_smooth.obj' 'ssa_misc.obj' 'mt19937ar.obj' 'libmwlapack.lib' 'libmwblas.lib'};
    end
else objfiles = {'ssa_matlab.o' 'ssa_kalman.o' 'ssa_smooth.o' 'ssa_misc.o' 'mt19937ar.o'};
end

%%%% Build MATLAB mex source files %%%%
mex(options{:}, 'kalman_int_c.c', objfiles{:});
mex(options{:}, 'loglik_int_c.c', objfiles{:});
mex(options{:}, 'loglik_grad_int_c.c', objfiles{:});
mex(options{:}, 'fastsmo_int_c.c', objfiles{:});
mex(options{:}, 'statesmo_int_c.c', objfiles{:});
mex(options{:}, 'disturbsmo_int_c.c', objfiles{:});
mex(options{:}, 'simsmo_int_c.c', objfiles{:});
mex(options{:}, 'sample_int_c.c', objfiles{:});
mex(options{:}, 'signal_int_c.c', objfiles{:});
mex(options{:}, 'weights_int_c.c', objfiles{:});
% gauss_int_c.c

%%%% Cleanup %%%%
if ispc, delete('*.obj'); else delete('*.o'); end
