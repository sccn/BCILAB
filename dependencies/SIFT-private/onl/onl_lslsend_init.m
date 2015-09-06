function outlet = onl_lslsend_init(streamname,samplingrate,numchans,datatype)
if nargin<2
    samplingrate = 0;
end
if nargin<3
    numchans = 1;
end
if nargin<4
    datatype = 'cf_string';  % can also be cft_float32, etc...
end
    
% instantiate the library
lib = lsl_loadlib();
info = lsl_streaminfo(lib,streamname,'MATLAB-serialized',numchans,samplingrate,datatype,'whatevergd456546');
outlet = lsl_outlet(info);
