function [f c] = setlinear(f, alpha)

%@SSFUNC/SETLINEAR Update the linear approximation.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

n   = size(alpha, 2);
c   = cell(1, length(f.f));
for i = 1 : length(f.f)
    horzmask    = f.horzmask(:, i);
    vertmask    = f.vertmask(:, i);
    dsubvec     = zeros(nnz(horzmask)*nnz(vertmask), n);
    c{i}        = zeros(nnz(vertmask), n);
    for t = 1 : n
        dot             = f.df{i}(alpha(horzmask, t), t);
        c{i}(:, t)      = f.f{i}(alpha(horzmask, t), t) - dot*alpha(horzmask, t);
        dsubvec(:, t)   = dot(:);
    end
    fmmask                      = false(size(f.ssmat));
    fmmask(vertmask, horzmask)  = true;
    f.ssmat                     = setdvec(f.ssmat, dsubvec, fmmask);
end
