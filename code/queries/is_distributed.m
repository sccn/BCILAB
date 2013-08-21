function x = is_distributed(x)
x = is_1d_distributed(x) || is_Nd_distributed(x) || is_discrete(x) || is_custom_distributed(x);
