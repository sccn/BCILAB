function onl_sendLSL(data,outlet)
% send a data structure over an LSL outlet

str = hlp_serialize(data);
outlet.push_sample({str});
