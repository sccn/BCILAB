function [p q P Q] = armadegree(y, varargin)

%ARMADEGREE Select the best fit ARMA or SARMA model to data.
%   [p q] = ARMADEGREE(y[, mean, mr])
%   [p q P Q] = ARMADEGREE(y, s[, mean, mr, ms])

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 2 || islogical(varargin{1})
    % Determine ARMA degree
    if nargin < 3, mr = 3; else mr = varargin{2}; end
    if nargin < 2, mean = false; else mean = varargin{1}; end

    BICr    = zeros(mr^2, 3);
    i       = 1;
    for p = 0 : mr
        for q = 0 : mr
            arma                = ssmodel('arma', p, q, mean);
            [arma logL output]  = estimate(y, arma, 0.1);
            BICr(i, :)          = [output.BIC p q];
            i                   = i + 1;
        end
    end
    [m i]   = min(BICr(:, 1));
    p       = BICr(i, 2);
    q       = BICr(i, 3);
    BICr    = sortrows(BICr, 1);
    fprintf(1, 'BIC\t\tp\tq\n');
    fprintf(1, '%6.4f\t%d\t%d\n', BICr');
    fprintf(1, '\n');
else
    % Determine SARMA degree
    if nargin < 5, ms = 1; else ms = varargin{4}; end
    if nargin < 4, mr = 3; else mr = varargin{3}; end
    if nargin < 3, mean = false; else mean = varargin{2}; end
    s       = varargin{1};

    BICs    = zeros(ms^2, 3);
    BICr    = zeros(mr^2, 3);

    p       = 3;
    q       = 0;

    i       = 1;
    for P = 0 : ms
        for Q = 0 : ms
            sarma               = ssmodel('sarima', p, 0, q, P, 0, Q, s, mean);
            [sarma logL output] = estimate(y, sarma, 0.1);
            BICs(i, :)          = [output.BIC P Q];
            i                   = i + 1;
        end
    end
    [m i]   = min(BICs(:, 1));
    P       = BICs(i, 2);
    Q       = BICs(i, 3);

    i       = 1;
    for p = 0 : mr
        for q = 0 : mr
            sarma               = ssmodel('sarima', p, 0, q, P, 0, Q, s, mean);
            [sarma logL output] = estimate(y, sarma, 0.1);
            BICr(i, :)          = [output.BIC p q];
            i                   = i + 1;
        end
    end
    [m i]   = min(BICr(:, 1));
    p       = BICr(i, 2);
    q       = BICr(i, 3);
    BICr    = sortrows(BICr, 1);
    fprintf(1, 'Regular part:\n');
    fprintf(1, 'BIC\t\tp\tq\n');
    fprintf(1, '%6.4f\t%d\t%d\n', BICr');
    fprintf(1, '\n');

    i       = 1;
    for P = 0 : ms
        for Q = 0 : ms
            sarma               = ssmodel('sarima', p, 0, q, P, 0, Q, s, mean);
            [sarma logL output] = estimate(y, sarma, 0.1);
            BICs(i, 1)          = output.BIC;
            i                   = i + 1;
        end
    end
    BICs    = sortrows(BICs, 1);
    P       = BICs(1, 2);
    Q       = BICs(1, 3);
    fprintf(1, 'Seasonal part:\n');
    fprintf(1, 'BIC\t\tP\tQ\n');
    fprintf(1, '%6.4f\t%d\t%d\n', BICs');
end

