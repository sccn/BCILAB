% GMMBayes Toolbox - Gaussian mixture model learning and classification
% Version 1.0
%
% Class learning interface
%   gmmb_create          - Construct new classifier with GMM PDFs
%                             Calls gmmb_em, gmmb_fj or gmmb_gem.
%
% Classification interface
%   gmmb_pdf             - multiclass and -variate (complex) GMM PDF
%   gmmb_weightprior     - multiply PDF values with constant priors
%   gmmb_normalize       - scale all row sums to one
%                             (useful in Bayesian classifier)
%   gmmb_decide          - make the classification decisions
%   gmmb_fracthresh      - threshold PDF values according to given distribution
%                             density quantile, i.e., find outliers
%
% Legacy Classification interface
%   gmmb_classify        - Run the classifier with bayesS struct.
%                             (Legacy interface)
%
% Estimation functions
%   gmmb_em              - The basic EM algorithm.
%   gmmb_fj              - The Figueiredo-Jain algorithm.
%   gmmb_gem             - Wrapper for the Vlassis code, Greedy EM algorithm.
%
% EM initializer functions
%   gmmb_em_init_cmeans1 - c-means clustering, uniform weight and covariance
%   gmmb_em_init_cmeans2 - c-means clustering, cluster weight and covariance
%   gmmb_em_init_fcm1    - Fuzzy c-means clustering, the Toolbox v0.1 method
%
% Helper functions
%   gmmb_cmvnpdf         - (Real and complex range) multivariate Gaussian pdf
%   gmmb_covfixer        - Force matrix to a valid covariance matrix.
%   gmmb_cmeans          - Simple c-means clustering
%   gmmb_mkcplx          - generate complex data with Gaussian distribution
%   gmmb_hist            - Create histS structure from user supplied data.
%   gmmb_generatehist    - Create histS structure from generated data,
%                             based on bayesS struct.
%   gmmb_lhood2frac      - map PDF values to density quantiles
%   gmmb_frac2lhood      - map density quantiles to PDF threshold values
%
% General functions
%   gmmb_version         - Return version string.
%   getargs              - Parse variable argument list into a struct.
%   warning_wrap         - Wrapper to allow Matlab R13 style warning calls
%                             in Matlab R12
%
% Vlassis code
%   gmmbvl_demo1d
%   gmmbvl_demo2d
%   gmmbvl_ellipse
%   gmmbvl_em
%   gmmbvl_em_gauss
%   gmmbvl_em_init_km
%   gmmbvl_em_step
%   gmmbvl_em_step_partial
%   gmmbvl_kmeans
%   gmmbvl_mixgen
%   gmmbvl_plot2
%   gmmbvl_rand_split
%   gmmbvl_sqdist
%   
% Demos
%   gmmb_demo01          - A simple demo presenting FJ-learning and
%                             Bayesian classification
%
% References:
%   [1] Duda, R.O., Hart, P.E, Stork, D.G, Pattern Classification,
%    2nd ed., John Wiley & Sons, Inc., 2001.
%   [2] Bilmes, J.A., A Gentle Tutorial of the EM Algorithm and its
%    Application to Parameter Estimation for Gaussian Mixture and Hidden
%    Markov Models
%    International Computer Science Institute, 1998
%   [3] Figueiredo, M.A.T., Jain, A.K., Unsupervised Learning on
%    Finite Mixture Models, IEEE transactions of pattern analysis and
%    machine intelligence, vol.24, no3, March 2002
%   [4] Vlassis, N., Likas, A., A Greedy EM Algorithm for Gaussian Mixture
%    Learning, Neural Processing Letters 15, Kluwer Academic Publishers, 2002.
%    http://carol.wins.uva.nl/~vlassis/research/learning/index_en.html
%   [5] Paalanen, P., Kamarainen, J.-K., Ilonen, J., Kälviäinen, H.,
%    Feature Representation and Discrimination Based on Gaussian Mixture Model
%    Probability Densities - Practices and Algorithms, Research Report 95,
%    Lappeenranta University of Technology, Department of Information
%    Technology, 2005.
%
% Authors:
%    Joni Kamarainen <Joni.Kamarainen@lut.fi>
%    Pekka Paalanen <pekka.paalanen@lut.fi>
%
% Acknowledgements:
%    All the gmmbvl_* functions and files are not written by the authors of
%    the Toolbox.
%    They are written by Dr. Nikos Vlassis and Sjaak Verbeek
%    unless otherwise noted in each file and are intellectual property of
%    their writers. The code is included "as is" [4], except the file and
%    function names are changed.
%
%    The Figueiredo-Jain algorithm implementation is partly based
%    on the code published in http://www.lx.it.pt/~mtf/
%
% Copyright:
%
%   The GMMBayes Toolbox is Copyright (C) 2003, 2004
%   by Joni Kamarainen and Pekka Paalanen except
%   the gmmbvl_*.m files which are copyrighted by their
%   respective authors.
%
%   The software package is free software; you can redistribute it
%   and/or modify it under terms of GNU General Public License as
%   published by the Free Software Foundation; either version 2 of
%   the license, or any later version. For more details see licenses
%   at http://www.gnu.org
%
%   The software package is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%   See the GNU General Public License for more details.
%
%   As stated in the GNU General Public License it is not possible to
%   include this software library or even parts of it in a proprietary
%   program without written permission from the owners of the copyright.
%   If you wish to obtain such permission, you can reach us by mail:
%
%      Department of Information Processing
%      Lappeenranta University of Technology
%      P.O. Box 20 FIN-53851 Lappeenranta
%      FINLAND
%
%  and by e-mail:
%      
%      joni.kamarainen@lut.fi
%      pekka.paalanen@lut.fi
%
%  Please, if you find any bugs contact authors.
%
%  Project home page: http://www.it.lut.fi/project/gmmbayes/
%
%   $Name:  $ $Revision: 1.4 $  $Date: 2005/04/14 10:33:34 $
%
