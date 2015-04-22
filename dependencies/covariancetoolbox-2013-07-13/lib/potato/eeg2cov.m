function COV = eeg2cov(Signal,window,overlap,method_cov)

    if nargin < 4
        method_cov = 'scm';
    end
    
    [Ne Nt] = size(Signal);
    
    Step = window-overlap;
    index = 1:Step:Nt-window-1;
    I = length(index);
    
    COV = zeros(Ne,Ne,I);

    for i=1:I
        COV(:,:,i) = covariances(Signal(:,index(i):index(i)+window),method_cov);
    end