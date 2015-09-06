function inlet = onl_lslreceive_init(streamname,timeout,maxbuffered,chunksize)

if nargin<2
    timeout = Inf; end
if nargin<3
    maxbuffered = []; end
if nargin<4
    chunksize = []; end
    
lib = lsl_loadlib();
result = {};

tmr = tic;
while isempty(result) && toc(tmr) < timeout
    result = lsl_resolve_byprop(lib,'name',streamname); end
if isempty(result)
    inlet = [];
else
    inlet = lsl_inlet(result{1},maxbuffered,chunksize); 
end
