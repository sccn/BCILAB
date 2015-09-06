function onl_lslsend(data,outlet)
% send a data structure over an LSL outlet
% tic;
tmp= hlp_serialize(data);
% toc
% tic;
outlet.push_sample({tmp});
% toc

