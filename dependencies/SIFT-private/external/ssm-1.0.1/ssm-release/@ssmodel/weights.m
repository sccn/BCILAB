function [wt_fil wt_smo] = weights(model, n, t0, com, varargin)

%@SSMODEL/WEIGHTS Compute weight functions for Kalman state filter and smoother.
%   [wtfil wtsmo] = WEIGHTS(model, n, t0[, com])
%   WEIGHTS(..., optname1, optvalue1, optname2, optvalue2, ...)
%       n is the data sequence length. (total number of time points)
%       t0 is the target time point.
%       Set com to true to calculate weights for components, else weights for
%           elements of the state vector is calculated.
%       wtfil is the filter weights. (m * n*p)
%       wtsmo is the smoother weights. (m * n*p)
%       (these weights operate on stacked vector of data y, (y(:)))
%       if com == false then alpha(:, t) == wt*y(:).
%       if com == true then
%           if p > 1 then ycom(:, t, i) == wt(:,:, i)*y(:) for each component i.
%           if p == 1 then ycom(:, t) == wt*y(:).
%           where ycom is the output of signal(alpha, model).

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 4, com = false; end

opt     = setopt([], varargin{:});

Zdyn    = ~issta(model.Z);
if opt.usec
    %% Weights %%
    if nargout > 1
        [wt_fil wt_smo] = weights_int_c(n, t0, getmat_c(model.H), ~issta(model.H), getmat_c(model.Z), Zdyn, ...
            getmat_c(model.T), ~issta(model.T), getmat_c(model.R), ~issta(model.R), getmat_c(model.Q), ~issta(model.Q), ...
            getmat_c(model.P1), com, model.mcom, opt.tol, opt.tol, opt.inv, true);
    else
        wt_fil          = weights_int_c(n, t0, getmat_c(model.H), ~issta(model.H), getmat_c(model.Z), Zdyn, ...
            getmat_c(model.T), ~issta(model.T), getmat_c(model.R), ~issta(model.R), getmat_c(model.Q), ~issta(model.Q), ...
            getmat_c(model.P1), com, model.mcom, opt.tol, opt.tol, opt.inv, true);
    end
else
    Zmat    = getmat(model.Z);

    %% Weights %%
    if nargout > 1, [omega omegaalpha] = weights_int(2, n, t0, ~issta(model.H), Zdyn, ~issta(model.T), ~issta(model.R), ~issta(model.Q), getmat(model.H), Zmat, getmat(model.T), getmat(model.R), getmat(model.Q), model.P1.mat, opt.tol);
    else omega = weights_int(1, n, t0, ~issta(model.H), Zdyn, ~issta(model.T), ~issta(model.R), ~issta(model.Q), getmat(model.H), Zmat, getmat(model.T), getmat(model.R), getmat(model.Q), model.P1.mat, opt.tol);
    end

    %% Data postprocessing %%
    m       = size(model.P1.mat, 1);
    wt_a    = [omega{:} zeros(m, model.p*(n-t0+1))];
    if nargout > 1, wt_alpha = [omegaalpha{:}]; end

    if com
        if Zdyn, if opt.usec, Z = Zmat(:,:, t0); else Z = Zmat{t0}; end, else Z = Zmat; end
        ncom        = length(model.mcom) - 1;
        if nargout > 1
            wt_fil  = zeros(model.p, n*model.p, ncom);
            wt_smo  = zeros(model.p, n*model.p, ncom);
            for i = 1 : ncom
                mask            = model.mcom(i)+1:model.mcom(i+1);
                wt_fil(:,:, i)  = Z(:, mask)*wt_a(mask, :);
                wt_smo(:,:, i)  = Z(:, mask)*wt_alpha(mask, :);
            end
        else
            wt_fil   = zeros(model.p, n*model.p, ncom);
            for i = 1 : ncom
                mask            = model.mcom(i)+1:model.mcom(i+1);
                wt_fil(:,:, i)  = Z(:, mask)*wt_a(mask, :);
            end
        end
        if model.p == 1
            wt_fil  = permute(wt_fil, [3 2 1]);
            if nargout > 1, wt_smo = permute(wt_smo, [3 2 1]); end
        end
    else
        wt_fil   = wt_a;
        if nargout > 1, wt_smo = wt_alpha; end
    end
end

