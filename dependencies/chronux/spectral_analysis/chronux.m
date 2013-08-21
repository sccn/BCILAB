function chronux
% This library performs time-frequency analysis (mostly using the
% multi-taper spectral estimation method) of univariate and multivariate
% data, both for continuous processes such as LFP/EEG and for point
% processes such as spike times. Point process can either be stored as
% times or as a binned process of counts. The routines in this library
% are named differently for the three cases. For calculations
% that can be performed for each of the three data types, we use suffixes
% c, pb, or pt to refer to continuous, point process binned counts, or
% point process times. For example, the spectrum calculation is performed
% mtspectrumc for continuous processes, mtspectrumpb for a binned point
% process, and mtspectrumpt for a point process consisting of times. There
% are also routines for calculating hybrid quantities involving one continuous
% and one point process. These are suffixed in a similar manner. For
% example, coherencycpb calculates the coherency between a binned point process
% and a continuous process. 
% 
% Certain variables are used repeatedly in this library.
%
% DATA
% data in most cases can be univariate or multivariate, and either point process, 
% or continuous.
%
%      Continuous data: Continuous data is assumed to be a matrix with 
%                       dimensions samples x channels/trials.
%
%      Point Process: A single time series of spike times can be in the form of 
%                     a column vector.
%                     Multichannel/trial spike time data is not amenable to this 
%                     storage format, since there are generally different 
%                     number of spikes in each channel/trial. Instead, 
%                     multichannel/trial spike data is stored in a structure 
%                     array. A structure is a matlab data object with various 
%                     fields. These fields contain the elements
%                       e.g. The command data=struct('times',[]); creates an empty 
%                            structure with field 'times'. Similarly, the command
%                            data=struct('times',[1 2 3]); creates the structure with
%                            the field 'times' containing integers 1, 2, and 3. 
%        
%                     We can also have a structure array (or an array of structures)
%                     defined for example, by
%                     data(1)=struct('times',rand(1,100)); and
%                     data(2)=struct('times',rand(1,200));
%                     This is a 2 dimensional structure array where the
%                     first field is a 100 dimensional random vector, and
%                     the second field is a 200 dimensional random vector.
%                     This format allows storage of multichannel point
%                     process times in a single variable data.
%                     
%                     The above holds for point processes stored as times.
%                     If instead, the point processes are binned, then one
%                     can use a matrix to represent them 
%                     
%
%      Summary: data - array of continuous data with dimensions time x channels
%                      structural array of spike times with dimensions
%                               equal to the number of channels
%                      1d array of spike times as a column vector
%                      array of binned spike counts with dimensions time x channels
%
% PARAMETERS:
% These are various parameters used in the spectral calculations. Since
% these parameters are used by most routines in Chronux, they are stored in
% a single structure params. The fields of params are
%
% tapers : precalculated tapers from dpss or in the one of the following
%          forms: 
%          (1) A numeric vector [TW K] where TW is the
%              time-bandwidth product and K is the number of
%              tapers to be used (less than or equal to
%              2TW-1). 
%          (2) A numeric vector [W T p] where W is the
%              bandwidth, T is the duration of the data and p 
%              is an integer such that 2TW-p tapers are used. In
%              this form there is no default i.e. to specify
%              the bandwidth, you have to specify T and p as
%              well. Note that the units of W and T have to be
%              consistent: if W is in Hz, T must be in seconds
%              and vice versa. Note that these units must also
%              be consistent with the units of params.Fs: W can
%              be in Hz if and only if params.Fs is in Hz.
%              The default is to use form 1 with TW=3 and K=5
%
%
% pad:   (padding factor for the FFT) - optional (can take values -1,0,1,2...). 
%         -1 corresponds to no padding, 0 corresponds to padding
%         to the next highest power of 2 etc.
%			  e.g. For N = 500, if PAD = -1, we do not pad; if PAD = 0, we pad the FFT
%			       to 512 points, if pad=1, we pad to 1024 points etc.
%			       Defaults to 0.
%
% Fs:sampling frequency.optional (default 1)
%
%
% fpass: frequencies in an fft calculation can range from 0 to Fs/2 where
%        Fs is the sampling frequency. Sometimes it may be useful to
%        compute fourier transforms (and resulting quantities like the
%        spectrum over a smaller range of frequencies). This is specified
%        by fpass, which can be in the form [fmin fmax] where fmin >=0 and
%        fmax<=Fs/2. optional (default [0 Fs/2])
%
% err=[errtype p] calculates theoretical error bars (confidence levels)
%                 when errtype=1 and jackknife error bars when errchk=2. In each case, the
%                 error is calculated at a p value specified by p. -
%                 optional (default 0)
%
% trialave: trialave controls whether or not to average over channels/trials for
%           multichannel/trial analyses. trialave=0 (default) implies no trial
%           averaging, trialave=1 implies that the quantity of interest is averaged
%           over channels/trials. optional (default 0)
% 
% Other parameters are discussed in individual routines as and when they
% are used.
