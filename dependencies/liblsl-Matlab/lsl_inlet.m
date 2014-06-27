classdef lsl_inlet < handle
	% A stream inlet.
	% Inlets are used to receive streaming data (and meta-data) from the lab network.
    %
    % Note:
    %   Most functions of the inlet can throw an error of type 'lsl:lost_error' if the stream
    %   read by this inlet has been lost irretrievably. Those functions that take a timeout can in
    %   addition throw an error of type 'lsl:timeout_error'.

    properties (Hidden)
        LibHandle = 0;      % this is a handle to the liblsl library (to call its functions)
        InletHandle = 0;    % this is a handle to an lsl_inlet object within the library.        
        ChannelCount = 0;   % copy of the inlet's channel count
        IsString = 0;       % whether this is a string-formatted inlet
        
        CorrectionBuffer = []; % buffer of time-correction values
        TimestampBuffer = [];  % buffer of time-stamps associated with the correction values
        RegressionCoeff = [];       % cached coefficients for linear regression
    end
    
    methods
        
        function self = lsl_inlet(info, maxbuffered, chunksize, recover)
            % Inlet = lsl_inlet(Streaminfo, MaxBuffered, ChunkSize, Recover)
            % Construct a new stream inlet from a resolved stream info.
            %
            % In:
            %   Streaminfo : A resolved stream info object (as coming from one of the lsl_resolve_* functions).
            %
            %   MaxBuffered : Optionally the maximum amount of data to buffer (in seconds if there is a nominal
            %				  sampling rate, otherwise x100 in samples). Recording applications want to use a fairly
            %				  large buffer size here, while real-time applications would only buffer as much as
            %				  they need to perform their next calculation. (default: 360)
            %   
            %   ChunkSize : Optionally the maximum size, in samples, at which chunks are transmitted
            %				(0 corresponds to the chunk sizes used by the sender).
            %				Recording applications can use a generous size here (leaving it to the network how
            %				to pack things), while real-time applications may want a finer (perhaps 1-sample) granularity.
            %               (default 0)
            %
            %   Recover : Try to silently recover lost streams that are recoverable (=those that that have a source_id set).
            %			  In all other cases (recover is false or the stream is not recoverable) functions may throw a
            %			  lost_error if the stream's source is lost (e.g., due to an app or computer crash).
            %             (default: 1)
            
            if ~exist('maxbuffered','var') || isempty(maxbuffered) maxbuffered = 360; end %#ok<*SEPEX>
            if ~exist('chunksize','var') || isempty(chunksize) chunksize = 0; end
            if ~exist('recover','var') || isempty(recover) recover = 1; end            
            self.LibHandle = info.LibHandle;
            self.ChannelCount = info.channel_count();
            self.IsString = strcmp(info.channel_format(),'cf_string');
            self.InletHandle = lsl_create_inlet(info.LibHandle,info.InfoHandle,maxbuffered,chunksize,recover);
        end
        
        
        function delete(self)
            % Destroy the inlet when it is deleted.
            % The inlet will automatically disconnect if destroyed.
            
            lsl_destroy_inlet(self.LibHandle,self.InletHandle);
        end

        
        function result = info(self, timeout)
            % Retrieve the complete information of the given stream, including the extended description.
            % Fullinfo = info(Timeout)
            %
            % Can be invoked at any point of the stream's lifetime.
            %
            % In:
            %   Timeout : Timeout of the operation, in seconds (default: 60).
            %
            % Out:
            %   Fullinfo : the full information, including description field, of the stream.
            
            if ~exist('timeout','var') || isempty(timeout) timeout = 60; end            
            result = lsl_streaminfo(self.LibHandle,lsl_get_fullinfo(self.LibHandle, self.InletHandle, timeout));
        end
		
        
        function open_stream(self, timeout)
            % Subscribe to the data stream from this moment onwards.
            % open_stream(Timeout)
            %
            % All samples pushed in at the other end from this moment onwards will be queued and 
		    % eventually be delivered in response to pull_sample() or pull_chunk() calls. 
		    % Pulling a sample without some preceding open_stream is permitted (the stream will then 
            % be opened implicitly).
            %
            % In:
            %   Timeout : Timeout of the operation (default: 60).
            %
            % Out:
            %   Fullinfo : the full information, including description field, of the stream.
            
            if ~exist('timeout','var') || isempty(timeout) timeout = 60; end            
            lsl_open_stream(self.LibHandle, self.InletHandle, timeout);
        end
        
        
        function close_stream(self)
            % Drop the current data stream.
            % close_stream()
            %
            % All samples that are still buffered or in flight will be dropped and the source will halt its buffering of data for this inlet.
            % If an application stops being interested in data from a source (temporarily or not) but keeps the outlet alive, it should
            % call close_streams() to not pressure the source outlet to buffer unnecessarily large amounts of data.
            
            lsl_close_stream(self.LibHandle, self.InletHandle);
        end
       
        
        function result = time_correction(self,timeout,postproc,winlen)
            % Retrieve an estimated time correction offset for the given stream.
            %
            % The first call to this function takes several miliseconds until a reliable first estimate is obtained.
            % Subsequent calls are very fast (and rely on periodic background updates).
            % The raw estimates should have a standard deviation under 1 ms (empirically within +/-0.5 ms).
            %
            % In:
            %   Timeout : Timeout to acquire the first time-correction estimate (default: 60).
            %
            %   PostProcessing : type of post-processing to apply; can be one of the following:
            %                    * 'raw' : return the raw time-correction estimates; these tend to 
            %                              have an average jitter of up to +/- 0.5ms (default)
            %                    * 'median' : return the median of measurements in a given time
            %                                 window (most robust; note that this becomes inaccurate
            %                                 for time windows > 1 minute if there is drift)
            %                    * 'linear' : return a linear fit in the given time window (fast and 
            %                                 accurate but sensitive to rare outliers, e.g., under
            %                                 network load spikes)
            %                    * 'trimmed' : a trimmed linear estimator; the two observations with
            %                                  the largest positive (resp. negative) deviation will 
            %                                  be removed
            %                    * 'robust' : return a robust linear fit in the given time window 
            %                                 (slowest but best tradeoff between accuracy and
            %                                 robustness)
            %                    * 'zero' : return zero (i.e., dismiss time-correction)
            %                    note: if you switch the post-processing you need to wait for up to
            %                          5 seconds until the estimator has updated
            %
            %   PostProcessingWindow : time window for post-processing, in seconds (default: 60)
            %
            % Out:
            %   Offset : The time correction estimate. If the first estimate cannot be obtained
            %            within the alloted time, the result is NaN.
            
            if ~exist('timeout','var') || isempty(timeout) timeout = 60; end
            if ~exist('postproc','var') || isempty(postproc) postproc = 'raw'; end
            if ~exist('winlen','var') || isempty(winlen) winlen = 60; end
            result = lsl_time_correction(self.LibHandle, self.InletHandle,timeout);            
            % do post-processing if requested (and don't attempt to post-process NaN values)
            if ~isnan(result) && ~strcmp(postproc,'raw')
                t0 = lsl_local_clock(self.LibHandle);
                % insert new measurement into buffer (but skip identical measurements to save
                % space/time since the time-correction values in LSL will only refresh every few
                % seconds)
                if isempty(self.CorrectionBuffer) || result ~= self.CorrectionBuffer(end)
                    self.CorrectionBuffer(end+1) = result;
                    self.TimestampBuffer(end+1) = t0;
                    % for each insert, check if we can remove any now-outdated values
                    remove_mask = self.TimestampBuffer < t0-winlen;
                    if any(remove_mask)
                        self.CorrectionBuffer(remove_mask) = [];
                        self.TimestampBuffer(remove_mask) = [];
                    end
                    % update regression model
                    switch postproc
                        case 'linear'
                            if length(self.CorrectionBuffer) > 1
                                % perform linear regression                        
                                self.RegressionCoeff = self.CorrectionBuffer / [ones(1,length(self.TimestampBuffer)); self.TimestampBuffer];
                            else
                                self.RegressionCoeff = [self.CorrectionBuffer(end) 0];
                            end
                        case 'trimmed'
                            if length(self.CorrectionBuffer) > 3
                                % perform trimmed linear regression (remove max/min value before linear fitting)
                                [dummy,maxidx] = max(self.CorrectionBuffer); %#ok<ASGLU>
                                [dummy,minidx] = min(self.CorrectionBuffer); %#ok<ASGLU>
                                indices = 1:length(self.CorrectionBuffer);
                                indices(indices==maxidx | indices==minidx) = [];
                                self.RegressionCoeff = self.CorrectionBuffer(indices) / [ones(1,length(indices)); self.TimestampBuffer(indices)];
                            else
                                self.RegressionCoeff = [self.CorrectionBuffer(end) 0];
                            end                            
                        case {'robust','robust_linear'}
                            if length(self.CorrectionBuffer) > 2
                                winsor_threshold = 0.0002;    % typical standard deviation of 0.2 ms
                                % perform linear fitting under Laplacian measurement noise (the below
                                % calculation tolerates 10-20% corruption by major outliers)
                                self.RegressionCoeff = robust_fit([ones(length(self.TimestampBuffer),1) self.TimestampBuffer']/winsor_threshold, self.CorrectionBuffer'/winsor_threshold,0.001,100);
                            else
                                self.RegressionCoeff = [self.CorrectionBuffer(end) 0];
                            end
                    end
                end
                % predict the time-correction
                switch postproc
                    case 'zero'
                        result = 0;
                    case 'median'
                        result = median(self.CorrectionBuffer);
                    case {'linear','trimmed','robust','robust_linear'}
                        result = self.RegressionCoeff(1) + self.RegressionCoeff(2) * t0;
                    otherwise
                        error('The given post-processing option is not valid: %s',postproc);
                end
            end
        end
        
        function [data,timestamp] = pull_sample(self,timeout,binary_blobs)
            % Pull a sample from the inlet and read it into an array of values.
            % [SampleData,Timestamp] = pull_sample(Timeout)
            %
            % In:
            %   Timeout : The timeout for this operation, if any. (default: 60)
            %
            % Out:
            %   SampleData : The sample's contents. This is either a numeric vector (type: double) if
            %                the stream holds numeric contents, or a cell array of strings if the stream
            %                is string-formatted.
            %
            %                Note: this is empty if the operation times out.
            %
            %   Timeout : Timeout for the operation. (default: 60)
            %
            %   BinaryBlobs : If true, strings will be received as uint8
            %                 arrays that may contain the \0 character.
            %                 Otherwise strings will be received as char,
            %                 which cannot contain \0's. (default: false)
            %
            % Notes:
            %   It is not a good idea to set the timeout to a very large value since MATLAB would 
            %   otherwise become unresponsive - better call pull_sample in a loop until it succeeds.
            
            if ~exist('timeout','var') || isempty(timeout) timeout = 60; end
            if self.IsString
                if ~exist('binary_blobs','var') || isempty(binary_blobs) binary_blobs = false; end
                if binary_blobs
                    [data,timestamp] = lsl_pull_sample_buf(self.LibHandle,self.InletHandle,self.ChannelCount,timeout);
                else
                    [data,timestamp] = lsl_pull_sample_str(self.LibHandle,self.InletHandle,self.ChannelCount,timeout);
                end
            else
                [data,timestamp] = lsl_pull_sample_d(self.LibHandle,self.InletHandle,self.ChannelCount,timeout);
            end
        end
        
        
        function [chunk,timestamps] = pull_chunk(self)
            % Pull a chunk of numeric samples and their timestamps from the inlet.
            % [ChunkData,Timestamps] = pull_chunk()
            %
            % This function obtains a chunk of data from the inlet; the chunk contains all samples
            % that have become available since the last chunk/sample was requested. Note that the 
            % result may be empty. For each returned sample there is also a timestamp being
            % returned.
            %
            % Out:
            %   ChunkData : The chunk contents; this is a MxN matrix with one column per returned
            %               sample (and as many rows as the stream has channels).
            %
            %   Timestamps : A vector of timestamps for the returned samples.
            
            [chunk,timestamps] = lsl_pull_chunk_d(self.LibHandle,self.InletHandle,self.ChannelCount);
        end
        
           
        
        % --- internal functions ---
 
        function h = get_libhandle(self)
            % get the library handle (e.g., to query the clock)
            h = self.LibHandle;
        end
        
        function x = robust_fit(A,y,rho,iters)
            % Perform a robust linear regression using the Huber loss function.
            % x = robust_fit(A,y,rho,iters)
            %
            % Input:
            %   A : design matrix
            %   y : target variable
            %   rho : augmented Lagrangian variable (default: 1)
            %   iters : number of iterations to perform (default: 1000)
            %
            % Output:
            %   x : solution for x
            %
            % Notes:
            %   solves the following problem via ADMM for x:
            %     minimize 1/2*sum(huber(A*x - y))
            %
            % Based on the ADMM Matlab codes also found at:
            %   http://www.stanford.edu/~boyd/papers/distr_opt_stat_learning_admm.html
            %
            %                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
            %                                2013-03-04
            
            if ~exist('rho','var')
                rho = 1; end
            if ~exist('iters','var')
                iters = 1000; end
            Aty = A'*y;
            L = sparse(chol(A'*A,'lower')); U = L';
            z = zeros(size(y)); u = z;
            for k = 1:iters
                x = U \ (L \ (Aty + A'*(z - u)));
                d = A*x - y + u;
                z = rho/(1+rho)*d + 1/(1+rho)*max(0,(1-(1+1/rho)./abs(d))).*d;
                u = d - z;
            end
        end
    end
       
end
