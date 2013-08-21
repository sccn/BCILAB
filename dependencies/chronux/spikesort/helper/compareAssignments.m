function [score, scoreIndep, p] = compareAssignments(assigns1, assigns2, showTables)

% COMPAREASSIGNMENTS   Computes a quality measure of the similarity between assignments.
%
% [score, scoreIfIndependent] = compareAssignments(assignments1, assignments2)
%     The inputs, 'assignments1' and 'assignments2', must be two column vectors of
%     the same length where each row contains an integer category label for the
%     corresponding sample.  The integer labels used in the assignment vectors need
%     have no intrinsic meaning (in particular, e.g., category 1 in 'assignments1'
%     has no relationship to category 1 in 'assignments2').
%
%     The first output, 'score', is a scalar between 0 and 1 that measures the
%     similarity between the two classifications.  A 'score' of 1 implies perfect
%     correspondance, ignoring actual labels.  For example, if all samples in
%     'assignments1' are labelled by 1 and relabelled as 2 in 'assignments2', the
%     'score' would be 1.  Deviations from this correspondance are penalized in a
%     a fashion that recognizes category splitting/merging and penalizes these less
%     than completely random redistribution.
%
%     The algorithm is motivated by a Chi^2 two-way classification; however, here we
%     return a similarity score rather than simply testing the hypothesis that the
%     classifications are independent.  The expected score if the classifications
%     were independent is returned as the second output, 'scoreIfIndependent', with
%     the standard Chi^2 two-way p-value returned as an optional third output (this
%     requires the statistics toolbox).  This p-value represents the probability that
%     the two assignments were independent.
%
%     Conceptually (though not computationally), the algorithm considers all N*(N-1)
%     pairs of data samples and counts pairs that cosegregate, where a pair of samples
%     is defined as cosegregating if they either share the same category in both
%     assignments or if they do not share category in either assignment.  For example,
%     consider the following assignments:
%            sample #          assignments1         assignments2
%               1                   1                     2
%               2                   1                     2
%               3                   2                     3
%               4                   1                     3
%     The pairs (1,2) and (1,3) cosegregate while the pair (1,4) does not (since they
%     share a label in 'assignments1' but not in 'assignments2').  'score' is the fraction
%     of pairs that cosegregate between the two assignments.
%
%     (An optional third boolean input argument 'showTables' (default 0) produces a graphical
%     output with the contingency table, conditional probabilities and marginals for the
%     assignments.  The 'score' described above is calculated efficiently using these matrices).

if ((size(assigns1, 2) > 1) || (size(assigns2, 2) > 1) || (size(assigns1,1) ~= size(assigns2, 1)))
    error('Error in assignment vectors.  The first two inputs must be column vectors of equal length.');
end

if ((nargin < 3) || (showTables == 0))    % if we're not doing graphics, this is more memory efficient.
    assigns1 = sortassignments(assigns1);
    assigns2 = sortassignments(assigns2);
    showTables = 0;
end

s = warning('MATLAB:divideByZero', 'off');

numSamples = size(assigns1, 1);
numCategories1 = length(unique(assigns1));
numCategories2 = length(unique(assigns2));

%  Construct classification table and marginals
joint = full(sparse(assigns1, assigns2, 1, max(assigns1), max(assigns2))) ./ numSamples;
marginal1 = sum(joint, 2);
marginal2 = sum(joint, 1);

% This somewhat cryptic expression computes the score described above.  i'll comment it
% later to explain.
score = (2 * joint(:)' * joint(:)) - sum(sum(joint' * joint)) - sum(sum(joint * joint'));
score = 1 + (numSamples / (numSamples - 1)) * score;

% Now get the score expected if the classifications were independent; we do this by
% reconstructing a joint under the assumption of independent classifications (i.e.,
% p(x,y) = p(x)p(y)) and then using the same mystery expression to find the score.
jointIndep = (marginal1 * marginal2);
scoreIndep = (2 * jointIndep(:)' * jointIndep(:)) ...
             - sum(sum(jointIndep' * jointIndep)) - sum(sum(jointIndep * jointIndep'));
scoreIndep = 1 + (numSamples / (numSamples-1)) * scoreIndep;

% if a p-value was requested, compute Chi^2
if (nargout > 2)
    X2 = numSamples .* (((joint - jointIndep).^2)./jointIndep);  % chi^2
    X2(isnan(X2)) = 0;  % (clean up divide by zeros)
    X2 = sum(X2(:));
    df = (numCategories1 - 1) * (numCategories2 - 1);  % degrees of freedom
    p = 1 - chi2cdf(X2,df);
end

% Optional graphical output
if (showTables)
    % construct conditional tables
    oneGivenTwo = joint ./ repmat(marginal2, [size(joint,1), 1]);
    oneGivenTwo(find(isnan(oneGivenTwo))) = 0;  % (deal with divide by zeros)
    twoGivenOne = joint ./ repmat(marginal1, [1, size(joint,2)]);
    twoGivenOne(find(isnan(twoGivenOne))) = 0; % (deal with divide by zeros)

    figure;
    subplot(2,2,1);  imagesc(joint);
    title('Two-Way Classification Table'); ylabel('Assignments 1'); xlabel('Assignments 2');
    subplot(2,2,2);  imagesc(oneGivenTwo);
    title('Assignments 1  given  Assignments 2'); ylabel('Assignments 1'); xlabel('Assignments 2');
    subplot(2,2,3);  imagesc(twoGivenOne);
    title('Assignments 2  given  Assignments 1'); ylabel('Assignments 1'); xlabel('Assignments 2');
    subplot(4,2,6);  bar(marginal1); axis tight;
    title('Assignments 1 Marginal');
    subplot(4,2,8);  bar(marginal2); axis tight;
    title('Assignments 2 Marginal');
    pixval on;
end

warning(s);