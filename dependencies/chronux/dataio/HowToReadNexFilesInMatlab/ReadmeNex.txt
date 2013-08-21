This package allows you to read data from .nex files.

The following functions are provided:

nex_info(filename) - displays information about a nex file
nex_ts(filename, varname) - gets the timestamps for the specified variable
nex_int(filename, varname) - gets the intervals for the specified interval variable
nex_cont(filename, varname) - gets the A/D data for the specified continuous variable
nex_wf(filename, varname) - gets the waveform values for the specified continuous variable
nex_marker(filename, varname) - gets marker data for the specified marker variable

See comments in the corresponding *.m files for more info.
The file test_nex.m contains a sample script that calls all the functions
against the provided Nex file text.nex.
