function out_signal = bci_Process(in_signal)
% Handle a BCI2000 Process event: map a raw input chunk to a raw output chunk.
% Out-Signal = bci_Process(In-Signal)
%
% In:
%   In-Signal  : incoming block of raw data [Channels x Samples]
%
% Out:
%   Out-Signal : block of BCI output data (i.e. predictions)
%
% See also:
%   www.bci2000.org/wiki/index.php/Programming_Reference:GenericFilter_Class
%
%
%                               Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                               2010-09-08

global meta_data in_samples out_samples next_outputs outrate maxcpu outform;


% This implementation makes heavy use of sample replication, as we assume that the BCI signal will be sampled at a significantly lower rate (e.g. 25Hz or even 0.2Hz in some cases)
% than the underling EEG signal (e.g. 1000Hz).

% our index into the signal block
blockidx = 1;

% while we are not done with the block...
while 1
    t0 = tic;

    % the samples that we need to append to the stream
    to_append = {};
    
    % while we have samples to read (from the input) and repeat (in the output) until the next output sample should be computed...
    while round(in_samples*(outrate/meta_data.srate)) <= out_samples
        % write (likely repeat) the next output into the out_signal
        out_signal(:,blockidx) = next_outputs{1};
        % get the next sample & mark for feeding to the BCI
        to_append{end+1} = in_signal(:,blockidx);
        in_samples = in_samples+1;
        blockidx = blockidx+1;
        % we are done with this block
        if blockidx > size(in_signal,2)
            % append and quit
            onl_append('bci2000_stream',double([to_append{:}]));
            return;
        end
    end

    % append the block that we accumulated to the stream...
	onl_append('bci2000_stream',double([to_append{:}]));
 
    % we are not yet done with the block, but need to compute an output sample at this point
    t1 = tic;

    % drop the next_output that we have just used up (over the last couple of samples)
    next_outputs = next_outputs(2:end);
    if ~isempty(next_outputs)
        % if there are some further queued next outputs (e.g., from the logic that can decide whether we are late and need to extrapolate (instead of compute) some samples),
        % then pick the first one of these, instead of making an actual prediction
    else
        prediction = onl_predict('bci2000_model',outform);
        % this prediction will be our next output (for the next couple of samples in the underlying signal's rate)
        next_outputs = {prediction};
    end
    % count up the number of computed (or estimated) samples
    out_samples = out_samples+1;
    
    % make sure that we are not exceeding our CPU budget
    computation_time = toc(t1);
    if maxcpu < 1
        pause(computation_time * (1 - maxcpu) / maxcpu); end
    
    % check how much processing time is left for this output sample
    time_spent = toc(t0);
    time_left = 1/outrate - time_spent;
    if time_left < 0
        % we exceeded our time budget: replicate samples some time ahead
        for k=1:ceil(-time_left*outrate)
            next_outputs{end+1} = next_outputs{end}; end
    end
end
