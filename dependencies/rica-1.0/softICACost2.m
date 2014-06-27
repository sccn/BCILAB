function [cost,grad] = softICACost2(theta, x, params)
% unpack weight matrix
W = reshape(theta, params.numFeatures, params.n);

% project weights to norm ball (prevents degenerate bases)
Wold = W;
W = l2rowscaled(W, 1);

% Forward Prop
h = W*x;
r = W'*h;

% Sparsity Cost
K = sqrt(params.epsilon + h.^2);
sparsity_cost = params.lambda * sum(sum(K));
K = 1./K;

% Coherence cost
cdiff = W*W' - eye(length(W));
coherence_cost = 0.5 * sum(sum(cdiff.^2));
W3grad = 2 * cdiff * W;

% Reconstruction Loss and Back Prop
diff = (r - x);
reconstruction_cost = 0.5 * sum(sum(diff.^2));
outderv = diff;

% compute the cost comprised of: 1) sparsity and 2) reconstruction
cost = sparsity_cost + reconstruction_cost + coherence_cost;

% Backprop Output Layer
W2grad = outderv * h';

% Baclprop Hidden Layer
outderv = W * outderv;
outderv = outderv + params.lambda * (h .* K);

W1grad = outderv * x';
Wgrad = W1grad + W2grad' + W3grad;

% unproject gradient for minFunc
grad = l2rowscaledg(Wold, W, Wgrad, 1);
grad = grad(:);