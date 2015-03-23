This is a MATLAB importer for .xdf files as created by the LabRecorder. To use it, simply copy it into a folder that is in your MATLAB path (or add this folder to the path). Then, in MATLAB type the line "doc load_xdf" to see the documentation.

Changelog:
1.14: Can handle broken files (fast-forwards to next boundary chunk), also eeg_load_xdf now handles event insertion 
      into segmented (interrupted) files properly.
1.13: Intelligently handles streams with mis-specified sampling rate.
