function ndiscarded = onl_lsl_bufmemcheck(~, ~, stream_inlet,mem_limit)

mibAvail = hlp_getAvailableMemory('MiB');
ndiscarded = 0;

while mibAvail < mem_limit
    % discard samples
    stream_inlet.pull_sample(-1,true);
    ndiscarded = ndiscarded + 1;
    mibAvail = hlp_getAvailableMemory('MiB');
end

