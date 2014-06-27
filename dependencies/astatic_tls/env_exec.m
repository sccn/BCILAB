eig(randn(32)*randn(32)); % calling this causes MATLAB to load the BLAS library (without it BLAS may get loaded too late and would fail on Linux)
