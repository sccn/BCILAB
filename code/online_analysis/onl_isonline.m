function tf = onl_isonline
% Test whether the calling code is running as part of the online processing (using onl_predict)
tf = hlp_resolve('is_online',false);