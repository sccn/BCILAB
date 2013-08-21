This package allows you to read data from .plx and .ddt files.

The following functions are provided:

plx_info(filename, all) - displays information about a plx file

plx_ts(filename, channel, unit) - gets the timestamps for the specified
	DSP channel and unit

plx_event_ts(filename, channel) - gets the timestamps for the specified
	external channel

plx_wave(filename, channel, unit) - gets the waveforms for the specified
	DSP channel and unit

plx_ad(filename, channel, unit) - gets the A/D data for the specified
	A/D channel


ddt(filename) - reads A/D data from a .ddt file

