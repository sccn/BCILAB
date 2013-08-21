This is a reference implementation of the artifact subspace reconstruction (ASR) method for artifact removal.
This code is published under the UC Copyright Notice (reproduced below).

Using the core implementation directly
======================================

The reference implementation may require the Signal Processing toolbox or alternatively a pre-computed 
IIR filter kernel if spectrally weighted statistics should be used (this is an optional feature; the filter 
kernel can be set to A=1 and B=1 to disable spectral weighting).

In all cases the signal that is passed to the method should be zero-mean, i.e., first high-pass filtered if necessary. 
In the following code samples the data is randomly generated with zero mean.

The general calling convention of the method is:

1)  Calibrate the parameters using the asr_calibrate function and some reasonably clean 
    calibration data, e.g., resting-state EEG, as in the following code snippet. 
    The recommended data length is ca. 1 minute or longer, the absolute minimum would be ca. 15 seconds. 
    There are optional parameters as documented in the function, in particular a tunable threshold parameter 
    governing the aggressiveness of the cleaning (although the defaults are sensible for testing purposes).

    calibdata = randn(20,10000);          % simulating 20-channel, 100-second random data at 100 Hz
    state = asr_calibrate(calibdata,100)  % calibrate the parameters of the method (sampling rate 100Hz)


2a) Apply the processing function to data that shall be cleaned, 
    either in one large block (entire recording) as follows...

    rawdata = randn(20,1000000);                  % simulating random data to process
    cleandata = asr_process(rawdata,100,state);   % apply the processing to the data (sampling rate 100Hz); see documentation of the function for optional parameters.

    Also note that for offline processing the method may be calibrated on the same data that should also be 
    cleaned, as long as the fraction of contaminated data is lower than a theoretical upper bound of 
    50% (a conservative empirical estimate would be approx. 30%). For extremely contaminated data one may 
    extract one or more clean segments from the data and calibrate on those.


2b) ... or alternatively apply the processing function to data in an online/incremental (chunk-by-chunk) 
    fashion, as follows:

    while 1
        newchunk = randn(20,50);                               % here using 0.5-second chunks; the chunk length is near-arbitrary and may vary, but transitions will be sharper for very short chunks.
        [cleanchunk,state] = asr_process(newchunk,100,state);  % apply the processing to the data and update the filter state (as in MATLAB's filter())
    end


License Terms
=============

UC Copyright Notice
This software is Copyright © 2013 The Regents of the University of California. All Rights Reserved.

Permission to copy, modify, and distribute this software and its documentation for educational,
research and non-profit purposes, without fee, and without a written agreement is hereby granted,
provided that the above copyright notice, this paragraph and the following three paragraphs appear
in all copies.

Permission to make commercial use of this software may be obtained by contacting:
Technology Transfer Office
9500 Gilman Drive, Mail Code 0910
University of California
La Jolla, CA 92093-0910
(858) 534-5815
invent@ucsd.edu 

This software program and documentation are copyrighted by The Regents of the University of
California. The software program and documentation are supplied "as is", without any accompanying
services from The Regents. The Regents does not warrant that the operation of the program will be
uninterrupted or error-free. The end-user understands that the program was developed for research
purposes and is advised not to rely exclusively on the program for any reason.

IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF
THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF THE UNIVERSITY OF CALIFORNIA HAS BEEN ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE. THE UNIVERSITY OF CALIFORNIA SPECIFICALLY DISCLAIMS ANY WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE. THE SOFTWARE PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF
CALIFORNIA HAS NO OBLIGATIONS TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
MODIFICATIONS.
