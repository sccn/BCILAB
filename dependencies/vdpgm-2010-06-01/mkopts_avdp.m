function opts = mkopts_avdp();

opts.algorithm = 'vdp';
opts.use_kd_tree = 1;
opts.initial_K = 1;
opts.do_greedy = 1;
opts.do_split = 0;
opts.do_merge = 0;
opts.do_sort = 1;
opts.initial_depth = 4;
opts.max_target_ratio = 0.1;
opts.recursive_expanding_depth = 2;
opts.recursive_expanding_threshold = 0.1;
opts.recursive_expanding_frequency = 3;

% Local Variables: ***
% mode: matlab ***
% End: ***
