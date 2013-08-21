************************************************************************
          RecPF: Reconstruction from partial Fourier data
************************************************************************

Copyright (C) 2008 Junfeng Yang, Yin Zhang and Wotao Yin

-- please note: this version is from February 2009
   There is a more recent version from May 2009, available at
   http://www.caam.rice.edu/~optimization/L1/RecPF/

1). Get Started
===================

       Run "demo_RecPF" to see how to use "RecPF" to reconstruct MR images from
    partial Fourier data.


2). Introduction
====================
   
       RecPF stands for (image/signal) Reconstruction from Partial Fourier data, 
    which may be contaminated with Gaussian noise, i.e.

                       f_p = F_p*ubar + omega,                              (1)

    where ubar is an original MR image, F_p is a partial Fourier operator, omega 
    is random noise and f_p is the (noisy) partial Fourier coefficients.

    RecPF reconstructs ubar by solving the following TVL1-L2 model:

              min_u aTV*TV(u) + aL1*||Psi'u||_1 + .5*||F_p*u - f_p||^2,     (2)

    where aTV, aL1 > 0 are regularization parameters and Psi is a sparsifying 
    basis. In "RecPF.m", one can set either aTV = 0 or aL1 = 0, but not both.


3). Details about "RecPF.m" 
=============================

    RecPF is used as:

           [U,Out] = RecPF(m,n,aTV,aL1,picks,f_p,PsiT,Psi,opts,varargin)

    Inputs:
         m, n      --- image size
         aTV, aL1  --- parameters in model (2);
         picks     --- a vector with its components the indices where Fourier 
                       coefficients are taken;
         f_p       --- partial Fourier data, a vector;
         PsiT      --- sparsifying basis (Haar wavelet as default)
         Psi       --- inverse Discrete wavelet transform
         opts      --- contains parameters for algorithm
                     * opts.mit_inn: maxium inner iteration number for each beta {default 30}
                     * opts.mit_out: maxium outer iteration number {default 10} 
                     * opts.tol_inn: inner iteration error tolerance {1.e-3}
                       (this option is effective only when "opts.stc = 2", see below)
                     * opts.tol_rel_inn: tolerance of relative change in an inner iteration {5.e-2}
                     * opts.tol_rel_out: tolerance of relative change in an outer iteration {1.e-2}
                       (the above two options are effective only when "opts.stc = 1", see below)
                     * opts.beta0:   initial penalty parameter {2^5}
                     * opts.beta_max: final penalty parameter  {2^15}
                     * opts.beta_rate: increase rate of beta
                     * opts.idisp: 0 or nonzero, display inner iteration info or not
                     * opts.recordf: 0 or 1, keep (or not) history of function values, TV and fidelity;
                     * opts.U0: starting poing {default "a least squares solution"} 
                     * opts.wbarflag: 1 (on) or 0 (off), controls wait bar.
                     * opts.stc: 1 (stopping criterion based on Relative Change) 
                              or 2 (stopping criterion based on optimality conditions) {default 1}
                     * opts.TVtype: 1 or 2, corresponding to anisotropic/isotropic TV in (2).
        varargin{1} --- a m x n matrix contains local weights of TV. Specifically, the local weigthed 
                      TV is discretized as:
                         TV(u) = sum_i w_i ||D_i u||. 
                      If all w_i == 1, then it is just the isotropic discritization of normal TV. 
                      We require all w_i > 0.
        varargin{2} --- true image (used to compute relative errors when it is present)
    Outputs:
            U   --- reconsctructed image
            Out --- a structrue contains
                    * Out.iter: total iteration number
                    * Out.Inner: inner iteration numbers for each beta
                    * Out.ftrue: function values at each iteration 
                    * Out.Fvalvscpu: function values decrease with respect to CPU time 
                    * Out.TVhist: total variation changes with iteration 
                    * Out.FIDhist: fidelity history
                          {the above four subfields are present only when opts.recordf = 1}
                    * Out.ITvsERR: relative error vs. iteration number 
                    * Out.CPUvsERR: CPU time vs. relative error 
                     {the above two subfields are present only when the original image is present}

    M-files in "utilities":

         dctPhi        --- DCT as sparsifying basis;
         identityPhi   --- Identity operator;
         MRImask       --- used to generate mask where Fourier coefficients are taken;
         snr           --- to compute Signal-to-Noise ratio;
         Wavedb1Phi    --- Haar wavelet transform;
         funcval       --- to compute function values of objective in (2).

4). IMPORTANT NOTES

    a). About penalty parameter "beta". A suitable final value of beta, i.e. 
    beta_max (a subfield of opts), depends on the range of intensity values 
    of the original image ubar. In the demo, the original image is scaled 
    into [0,1]. In this circumstance, we used default value opts.beta_max 
    = 2^15. If the original image has intensity values between 0 and 255, 
    a suitable value of beta_max would be much smaller. If, in this case, 
    one does not change the beta sequence accordingly, RecPF could be slow.

    b). The residue based on optimality conditions is NOT scalling invariant. 
    As a result, if the final value of beta is inappropriate, RecPF could be 
    slow if users happen to stop RecPF by optimality, i.e. "opts.stc = 2". 

    We suggest that users scale the intensity values of an input image to 
    [0,1] and set beta_max = 2^15 (as given in the default setting of RecPF)
    before calling RecPF.

5). Reference
====================

        For algorithmic details, such as continution on penalty parameters and 
     optimality conditions, see references:

        J. Yang, Y. Zhang and W. Yin, "A fast TVL1-L2 minimization algorithm 
    for signal reconstruction from partial Fourier data", Tech. Report 08-27, 
    CAAM, Rice University.

6). Contact Information
=======================

    RecPF is available at: http://www.caam.rice.edu/~optimization/L1/RecPF/

    Comments or suggestions? Please feel free to e-mail:

    Junfeng Yang,  Dept. Math., Nanjing Univ.,  <jfyang2992@yahoo.com.cn>


7).  Copyright Notice
=======================

       RecPF is free software; you can redistribute it and/or modify it under 
    the terms of the GNU General Public License as published by the Free 
    Software Foundation; either version 3 of the License, or (at your option) 
    any later version.

    This program is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details at
    <http://www.gnu.org/licenses/>. 
