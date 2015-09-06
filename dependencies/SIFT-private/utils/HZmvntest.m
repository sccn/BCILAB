function [HZmvntest] = HZmvntest(X,c,alpha)
%HZMVNTEST. Henze-Zirkler's Multivariate Normality Test.
% Henze and Zirkler (1990) introduce a multivariate version of the univariate
% There are many tests for assessing the multivariate normality in the 
% statistical literature (Mecklin and Mundfrom, 2003). Unfortunately, there
% is no known uniformly most powerful test and it is recommended to perform
% several test to assess it. It has been found that the Henze and Zirkler
% test have a good overall power against alternatives to normality.
%
% The Henze-Zirkler test is based on a nonnegative functional distance that
% measures the distance between two distribution functions. This distance
% measure is,
%                           _ 
%                          /              
%          D_b(P,Q)K  =   /  |P(T) - Q(t)|^2 f_b(t)dt
%                       _/   
%                    R^d
%
% where P(t) [the characteristic function of the multivariate normality
% distribution] and Q(t) [the empirical characteristic function] are the
% Fourier transformations of P and Q, and f_b is the weight or kernel
% function. Its normal density is,
%
%                                        -|t|^2
%          f_b(t) = (2*pi*b^2)^-p/2*exp(--------), t E R^p
%                                         2*b^2
%
% where p = number of variables and |t| = (t't)^0.5.
%
% The smoothing parameter b depends on n (sample size) as,
%
%                       1     2p + 1
%            b_p(n) = ------ (-------)^1/(p+4)*n^1/(p+4) 
%                     sqrt(2)    4
%
% The Henze-Zirkler statistic is approximately distributed as a lognormal.
% The lognormal distribution is used to compute the null hypothesis 
% probability.
%                      _n _  _n _ 
%                  1   \     \            b^2
%        W_n,b = ----- /_ _  /_ _ exp(- ------- Djk) - 2(1+b^2)^-p/2*
%                 n^2  j=1   k=1           2
%
%                 _n _
%              1  \             b^2
%             --- /_ _ exp(- --------- Dj) + (1+2b^2)^-p/2
%              n  j=1         2(1+b^2)
%
%        HZ = n*(4  1{S is singular} + W_n,b 1{S is nonsingular})
%
% where  W_n,b = weighted L^2-distance
%        Djk = (X_j - X_k)*inv(S)*(X_j - X_k)'
%        Dj = (X - MX)*inv(S)*(X - MX)'
%        X = data matrix
%        MX = sample mean vector
%        S = covariance matrix (normalized by n)
%        HZ = Henze-Zirkler statistic test (it is known in literature as
%             T_n,b)
%        1{.} = stands for the indicator function
%
% According to Henze-Wagner (1997), this test has the desirable properties
% of,
%    --affine invariance
%    --consistency against each fixed nonnormal alternative distribution
%    --asymptotic power against contiguous alternatives of order n^-1/2
%    --feasibility for any dimension and any sample size
%
% If the data is multivariate normality, the test statistic HZ is 
% approximately lognormally distributed. It proceeds to calculate the mean,
% variance and smoothness parameter. Then, mean and variance are 
% lognormalized and the z P-value is estimated.
%
% Also, for all the interested people we provide the q_b,p(alpha) = 
% HZ(1-alpha)-quantile of its lognormal distribution (critical value) 
% having expectation E(T_b(p)) and variance var(T_b(p)).
%
% Syntax: function HZmvntest(X,alpha) 
%      
% Inputs:
%      X - data matrix (size of matrix must be n-by-p; data=rows,
%          indepent variable=columns) 
%      c - covariance normalized by n (=1, default)) or n-1 (~=1)
%  alpha - significance level (default = 0.05)
%
% Output:
%        - Henze-Zirkler's Multivariate Normality Test
%
% Example: From the Table 11.5 (Iris data) of Johnson and Wichern (1992,
%          p. 562), we take onlt the Isis setosa data to test if it has a 
%          multivariate normality distribution using the Doornik-Hansen
%          omnibus test. Data has 50 observations on 4 independent 
%          variables (Var1 (x1) = sepal length; Var2 (x2) = sepal width; 
%          Var3 (x3) = petal length; Var4 (x4) = petal width. 
%
%                      ------------------------------
%                        x1      x2      x3      x4       
%                      ------------------------------
%                       5.1     3.5     1.4     0.2     
%                       4.9     3.0     1.4     0.2     
%                       4.7     3.2     1.3     0.2     
%                       4.6     3.1     1.5     0.2     
%                       5.0     3.6     1.4     0.2     
%                       5.4     3.9     1.7     0.4     
%                        .       .       .       .      
%                        .       .       .       .      
%                        .       .       .       .      
%                       5.1     3.8     1.6     0.2     
%                       4.6     3.2     1.4     0.2     
%                       5.3     3.7     1.5     0.2     
%                       5.0     3.3     1.4     0.2     
%                      ------------------------------
%
% Total data matrix must be:
% You can get the X-matrix by calling to iris data file provided in
% the zip as 
%          load path-drive:irisetdata 
% 
% Calling on Matlab the function: 
%          HZmvntest(X)
%
% Answer is:
%
% Henze-Zirkler's Multivariate Normality Test
% -------------------------------------------------------------------
% Number of variables: 4
% Sample size: 50
% -------------------------------------------------------------------
% Henze-Zirkler lognormal mean: -0.2794083
% Henze-Zirkler lognormal variance: 0.1379069
% Henze-Zirkler statistic: 0.9488453
% P-value associated to the Henze-Zirkler statistic: 0.0499536
% With a given significance = 0.050
% -------------------------------------------------------------------
% Data analyzed do not have a normal distribution.
% -------------------------------------------------------------------
%
% Created by A. Trujillo-Ortiz, R. Hernandez-Walls, K. Barba-Rojo, 
%             and L. Cupul-Magana
%             Facultad de Ciencias Marinas
%             Universidad Autonoma de Baja California
%             Apdo. Postal 453
%             Ensenada, Baja California
%             Mexico.
%             atrujo@uabc.mx
%
% Copyright. December 5, 2007.
%
% To cite this file, this would be an appropriate format:
% Trujillo-Ortiz, A., R. Hernandez-Walls, K. Barba-Rojo and L. Cupul-Magana.
%   (2007). HZmvntest:Henze-Zirkler's Multivariate Normality Test. A MATLAB
%   file. [WWW document]. URL http://www.mathworks.com/matlabcentral/
%   fileexchange/loadFile.do?objectId=17931
%
% References:
% Henze, N. and Zirkler, B. (1990), A Class of Invariant Consistent Tests
%      for Multivariate Normality. Commun. Statist.-Theor. Meth., 19(10):
%      3595–3618.
% Henze, N. and Wagner, Th. (1997), A New Approach to the BHEP tests for 
%      multivariate normality. Journal of Multivariate Analysis, 62:1-23.
% Johnson, R. A. and Wichern, D. W. (1992), Applied Multivariate Statistical
%      Analysis. 3rd. ed. New-Jersey:Prentice Hall.
% Mecklin, C. J. and Mundfrom, D. J. (2003), On Using Asymptotic Critical
%      Values in Testing for Multivariate Normality. InterStat: Statistics
%      on Internet. http://interstat.statjournals.net/YEAR/2003/abstracts/
%      0301001.php
%

if (nargin < 3),
    alpha = 0.05; %default
end

if (alpha <= 0 || alpha >= 1),
    fprintf(['Warning: Significance level error; must be 0 <'...
        ' alpha < 1 \n']);
    return;
end

if (nargin < 2),
    c = 1;
end

if (nargin < 1),
    error('Requires at least one input argument.');
    return;
end

if c == 1  %covariance matrix normalizes by (n) [=default]
    S = cov(X,1);
else   %covariance matrix normalizes by (n-1)
    S = cov(X);
end

[n,p] = size(X);

difT = (X - repmat(mean(X),n,1)); 

Dj = diag(difT*inv(S)*difT');  %squared-Mahalanobis' distances

Y = X*inv(S)*X';

%Djk=[];
%for j = 1:n,
%    for k = 1:n,
%        d = [Y(j,j)- 2*Y(j,k)+Y(k,k)];
%        Djk = [Djk;d];
%    end
%end

%Djk = reshape(Djk,n,n);

Djk = - 2*Y' + diag(Y')*ones(1,n) + ones(n,1)*diag(Y')'; %this procedure was
                       %taken into account in order to '..avoiding loops and
                       %it (the file) runs much faster' we thank to Johan
                       %(J.D.) for it valuabe comment (15/12/08) to improve 
                       %it.

b = 1/(sqrt(2))*((2*p + 1)/4)^(1/(p + 4))*(n^(1/(p + 4))); %smoothing
                                                                 %parameter
if (rank(S) == p),    
    HZ = n * (1/(n^2) * sum(sum(exp( - (b^2)/2 * Djk))) - 2 *...
        ((1 + (b^2))^( - p/2)) * (1/n) * (sum(exp( - ((b^2)/(2 *...
        (1 + (b^2)))) * Dj))) + ((1 + (2 * (b^2)))^( - p/2)));
else
    HZ = n*4;
end

wb = (1 + b^2)*(1 + 3*b^2); % = 1 + 4*b^2 + 3*b^4

a = 1 + 2*b^2;

mu = 1 - a^(- p/2)*(1 + p*b^2/a + (p*(p + 2)*(b^4))/(2*a^2)); %HZ mean

si2 = 2*(1 + 4*b^2)^(- p/2) + 2*a^( - p)*(1 + (2*p*b^4)/a^2 + (3*p*...
    (p + 2)*b^8)/(4*a^4)) - 4*wb^( - p/2)*(1 + (3*p*b^4)/(2*wb) + (p*...
    (p + 2)*b^8)/(2*wb^2)); %HZ variance

pmu = log(sqrt(mu^4/(si2 + mu^2))); %lognormal HZ mean
psi = sqrt(log((si2 + mu^2)/mu^2)); %lognormal HZ variance

P = 1 - logncdf(HZ,pmu,psi);%P-value associated to the HZ statistic

%(1-alpha)-quantile of the lognormal distribution
%q = mu*(1+si2/mu^2)^(-1/2)*exp(norminv(1-alpha)*sqrt(log(1+si2/mu^2))); 

disp(' ')
disp('Henze-Zirkler''s Multivariate Normality Test')
disp('-------------------------------------------------------------------')
fprintf('Number of variables: %i\n', p);
fprintf('Sample size: %i\n', n);
disp('-------------------------------------------------------------------')
fprintf('Henze-Zirkler lognormal mean: %3.7f\n', pmu);
fprintf('Henze-Zirkler lognormal variance: %3.7f\n', psi);
fprintf('Henze-Zirkler statistic: %3.7f\n', HZ);
fprintf('P-value associated to the Henze-Zirkler statistic: %3.7f\n', P);
fprintf('With a given significance = %3.3f\n', alpha);
disp('-------------------------------------------------------------------')
if P >= alpha;
    disp('Data analyzed have a normal distribution.');
else
    disp('Data analyzed do not have a normal distribution.');
end
disp('-------------------------------------------------------------------')

return,