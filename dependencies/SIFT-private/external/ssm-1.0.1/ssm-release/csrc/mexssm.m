function mexssm(varargin)
mex(varargin{:}, '-outdir', '../@ssmat', 'getmat_c.c');
