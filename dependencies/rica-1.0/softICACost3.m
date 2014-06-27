function [cost,grad] = softICACost3(theta, x, params)

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

% Reconstruction Loss and Back Prop
diff = (r - x);
reconstruction_cost = params.theta * 0.5 * sum(sum(diff.^2));
outderv = diff;

% Location constraint
location_cost = 0;
A = W';
location_grad = zeros(size(A));
for c=1:length(params.subspaces)
    subterm = params.subspaces{c} * A(:,c);
    location_cost = location_cost + A(:,c)' * subterm;
    location_grad(:,c) = 2 * subterm;
end

% compute the cost comprised of: 1) sparsity, 2) reconstruction, 3) location
cost = sparsity_cost + reconstruction_cost + location_cost*size(x,2);
fprintf('%.2f\t%.2f\t%.2f\n',sparsity_cost,reconstruction_cost,location_cost*size(x,2));

% Backprop Output Layer
W2grad = outderv * h';

% Baclprop Hidden Layer
outderv = W * outderv;
outderv = outderv + params.lambda * (h .* K);

W1grad = outderv * x';
Wgrad = W1grad + W2grad' + location_grad';

% unproject gradient for minFunc
grad = l2rowscaledg(Wold, W, Wgrad, 1);
grad = grad(:);
