function opts = mkopts_csb(T);

opts.algorithm = 'csb';
opts.use_kd_tree = 0;
opts.initial_K = T;
opts.do_greedy = 0;
opts.do_split = 0;
opts.do_merge = 0;
opts.do_sort = 0;

% Local Variables: ***
% mode: matlab ***
% End: ***
