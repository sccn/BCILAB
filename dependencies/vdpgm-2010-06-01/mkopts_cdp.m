function opts = mkopts_cdp(K);

opts.algorithm = 'cdp';
opts.use_kd_tree = 0;
opts.initial_K = K;
opts.do_greedy = 0;
opts.do_split = 0;
opts.do_merge = 0;
opts.do_sort = 0;

% Local Variables: ***
% mode: matlab ***
% End: ***
