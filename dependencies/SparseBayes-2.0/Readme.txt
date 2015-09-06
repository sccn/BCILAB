SPARSEBAYES Matlab Toolbox

Version 2.00

-- INTRODUCTION --

"SparseBayes" is a package of Matlab functions designed to implement
an efficient learning algorithm for "Sparse Bayesian" models.

The "Version 2" package is an expanded implementation of the algorithm
detailed in:
    
    Tipping, M. E. and A. C. Faul (2003). "Fast marginal likelihood
    maximisation for sparse Bayesian models." In C. M. Bishop and
    B. J. Frey (Eds.), Proceedings of the Ninth International Workshop
    on Artificial Intelligence and Statistics, Key West, FL, Jan 3-6.

This paper, the accompanying code, and further information regarding
Sparse Bayesian models and the related "relevance vector machine" may
be obtained from:

http://www.relevancevector.com


-- CODE STATUS --

March 2009: The code is currently at Version 2.0, and has seen limited
testing under Matlab Version 7.4.


-- GETTING STARTED --

The SparseBayes distribution comes with a basic user manual: see
SB2_Manual.pdf.

In summary, there are a number of files supplied in the SB2
distribution, but the user should only need to call a subset of them.

To set up options and parameters:

   SB2_UserOptions.m
   SB2_ParameterSettings.m

The main algorithm

    SparseBayes.m

Type "help SparseBayes" etc. at the Matlab prompt for further details
on how to use these functions. 

There is a also a simple demonstration program, showing how
SparseBayes may be used, in:

    SparseBayesDemo.m
	

-- LICENCE --

Note that the "SparseBayes" software is supplied subject to version 2
of the "GNU General Public License" (detailed in the file "licence.txt").

SPARSEBAYES is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation; either version 2 of the License, or (at your
option) any later version.

SPARSEBAYES is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with SPARSEBAYES in the accompanying file "licence.txt"; if not,
write to the Free Software Foundation, Inc., 51 Franklin St, Fifth
Floor, Boston, MA 02110-1301 USA


-- ACKNOWLEDGEMENTS --

The author would like to thank Mark Hatton, Anita Faul, Ian Nabney,
Arnulf Graf and Gavin Cawley for their assistance in producing this
code.


--

Mike Tipping
www.relevancevector.com
m a i l [at] m i k e t i p p i n g . c o m

This README file (Readme.txt) was created on 13-Mar-2009 at  8:58 AM.
