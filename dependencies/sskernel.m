function [y,t,optw,W,C,confb95,yb] = sskernel(x,tin,W)
% [y,t,optw,W,C,confb95,yb] = sskernel(x,t,W)
%
% Function `sskernel' returns an optimized kernel density estimate 
% using a Gauss kernel function.
%
% Examples:
% >> x = 0.5-0.5*log(rand(1,1e3)); t = linspace(0,3,1000);
% >> [y,t,optw] = sskernel(x,t);
% This example produces a vector of kernel density estimates, y, at points
% specified in a vector t, using an optimized bandwidth, optw (a standard 
% deviation of a normal density function).
% 
% >> sskernel(x);
% By calling the function without output arguments, the estimated density 
% is displayed along with 95% bootstrap confidence intervals.
%
% Input arguments:
% x:    Sample data vector. 
% tin (optinal):
%       Points at which estimation are computed. Please use fine resolution
%       to obtain a correct optimal bandwidth.
% W (optinal): 
%       A vector of kernel bandwidths. 
%       If W is provided, the optimal bandwidth is selected from the 
%       elements of W.
%       * Do not search bandwidths smaller than a sampling resolution of data.
%       If W is not provided, the program searches the optimal bandwidth
%       using a golden section search method. 
%
% Output arguments:
% y:    Estimated density
% t:    Points at which estimation was computed.
%       The same as tin if tin is provided. 
%       (If the sampling resolution of tin is smaller than the sampling 
%       resolution of the data, x, the estimation was done at smaller
%       number of points than t. The results, t and y, are obtained by 
%       interpolating the low resolution sampling points.)
% optw: Optimal kernel bandwidth.
% W:    Kernel bandwidths examined. 
% C:    Cost functions of W.
% conf95:
%       Bootstrap confidence intervals.
% yb:   Booststrap samples.
%
% 
% Usage:
% >> [y,t,optw] = sskernel(x);
% When t is not given in the input arguments, i.e., the output argument t 
% is generated automatically.
%
% >> W = linspace(0.01,1,20);
% >> [y,t,optw] = sskernel(x,t,W);
% The optimal bandwidth is selected from the elements of W.
%
% >> [y,t,optw] = sskernel(x,t,0.1);
% If the density estimate with a given bandwidth, simply put a scalar value
% as W. The computation is faster than the built-in function, ksdensity.
%
% >> [y,t,optw,confb95,yb] = sskernel(x);
% This additionally computes 95% bootstrap confidence intervals, confb95.
% The bootstrap samples are provided as yb.
% 
%
% Optimization principle:
% The optimal bandwidth is obtained as a minimizer of the formula, 
% sum_{i,j} \int k(x - x_i) k(x - x_j) dx  -  2 sum_{i~=j} k(x_i - x_j), 
% where k(x) is the kernel function, according to
%
% Hideaki Shimazaki and Shigeru Shinomoto
% Kernel Bandwidth Optimization in Spike Rate Estimation 
% Journal of Computational Neuroscience 2010
% http://dx.doi.org/10.1007/s10827-009-0180-4
%
% The above optimization is based on a principle of minimizing 
% expected L2 loss function between the kernel estimate and an unknown 
% underlying density function. An assumption is merely that samples 
% are drawn from the density independently each other. 
%
% For more information, please visit 
% http://2000.jukuin.keio.ac.jp/shimazaki/res/kernel.html
%
% See also SSVKERNEL, SSHIST
% 
% Bug fix
% 131004 fixed a problem for large values
%
% Hideaki Shimazaki 
% http://2000.jukuin.keio.ac.jp/shimazaki

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters Settings
x = reshape(x,1,numel(x));

if nargin == 1
    T = max(x) - min(x);
    [mbuf,nbuf,dt_samp] = find( sort(diff(sort(x))),1,'first');
    tin = linspace(min(x),max(x), min(ceil(T/dt_samp),1e3));
    t = tin;
    x_ab = x( logical((x >= min(tin)) .*(x <= max(tin))) ) ;
else
    T = max(tin) - min(tin);    
    x_ab = x( logical((x >= min(tin)) .*(x <= max(tin))) ) ;
    [mbuf,nbuf,dt_samp] = find( sort(diff(sort(x_ab))),1,'first');

    if dt_samp > min(diff(tin))
        t = linspace(min(tin),max(tin), min(ceil(T/dt_samp),1e3));
    else
        t = tin;
    end
end

dt = min(diff(t));

% Create a finest histogram
y_hist = histc(x_ab,t-dt/2);
L = length(y_hist);
N = sum(y_hist);
y_hist = y_hist/N/dt;   %density

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute a Cost Function

if nargin >= 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global search
C = zeros(1,length(W));
C_min = Inf;

for k = 1: length(W)
	w = W(k);     
    [C(k) yh] = CostFunction(y_hist,N,w,dt);
    
    if C(k) < C_min
        C_min = C(k);
        optw = w;
        y = yh;
    end
end

else
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Golden section search on a log-exp scale
% Initialize
Wmin = 2*dt; Wmax = 1*(max(x) - min(x));

tol = 10^-5; 
phi = (sqrt(5) + 1)/2;        %golden ratio
%logexp = @(x) log(1+exp(x));
%ilogexp = @(x) log(exp(x)-1);

%a = Wmin; b = Wmax;
a = ilogexp(Wmin); b = ilogexp(Wmax);

c1 = (phi-1)*a + (2-phi)*b;
c2 = (2-phi)*a + (phi-1)*b;

f1 = CostFunction(y_hist,N,logexp(c1),dt);
f2 = CostFunction(y_hist,N,logexp(c2),dt);

k = 1;
while abs(b-a) > tol*(abs(c1)+abs(c2)) && k <= 20
	if (f1 < f2)    
        b = c2;
        c2 = c1;

        c1 = (phi - 1)*a + (2 - phi)*b;
        
        f2 = f1;
        [f1 yh1] = CostFunction(y_hist,N,logexp(c1),dt);
        
        W(k) = logexp(c1);
        C(k) = f1;
        optw = logexp(c1);
        y = yh1./sum(yh1.*dt);  %make the final output a density
    else
        a = c1;
        c1 = c2;
        
        c2 = (2 - phi)*a + (phi - 1)*b;
        
        f1 = f2;
        [f2 yh2] = CostFunction(y_hist,N,logexp(c2),dt);
        
        W(k) = logexp(c2);
        C(k) = f2;
        optw = logexp(c2);
        y = yh2./sum(yh2.*dt);
    end
    
    k = k + 1;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bootstrap Confidence Intervals
if nargout == 0 || nargout >= 6
    nbs = 1e3;        %number of bootstrap samples
    yb = zeros(nbs,length(tin));

    for i = 1: nbs,
        %y_histb = poissrnd(y_hist*dt*N)/dt/N;
    
        idx = ceil(rand(1,N)*N);
        xb = x_ab(idx);
        y_histb = histc(xb,t-dt/2)/dt/N;
    
        yb_buf = fftkernel(y_histb,optw/dt);
        yb_buf = yb_buf / sum(yb_buf*dt);
        
        yb(i,:) = interp1(t,yb_buf,tin);
    end

    ybsort = sort(yb);
    y95b = ybsort(floor(0.05*nbs),:);
    y95u = ybsort(floor(0.95*nbs),:);
    
    confb95 = [y95b; y95u];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Return results
y = interp1(t,y,tin);
t = tin;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display results
if nargout == 0
            hold on;

            line([t; t],[y95b; y95u]...
                    ,'Color',[7 7 7]/8,'LineWidth',1 );
            plot(t,y95b,'Color',[7 7 7]/9,'LineWidth',1);
            plot(t,y95u,'Color',[7 7 7]/9,'LineWidth',1);

            plot(t,y,'Color',[0.9 0.2 0.2],'LineWidth',2);

            grid on;
            ylabel('density');
            set(gca,'TickDir','out');  
else
    if nargin >= 4
        if strcmp(option,'Visible')
            hold on; 
            
            if nargout >= 6
                line([t; t],[y95b; y95u]...
                    ,'Color',[7 7 7]/8,'LineWidth',1 );
                plot(t,y95b,'Color',[7 7 7]/9,'LineWidth',1);
                plot(t,y95u,'Color',[7 7 7]/9,'LineWidth',1);
            end

            plot(t,y,'Color',[0.9 0.2 0.2],'LineWidth',1);

            grid on;
            ylabel('density');
            set(gca,'TickDir','out'); 
            
        end
    end
end


function [C yh] = CostFunction(y_hist,N,w,dt)
yh = fftkernel(y_hist,w/dt);  %density

%formula for density
C = sum(yh.^2)*dt - 2* sum(yh.*y_hist)*dt...
        + 2*1/sqrt(2*pi)/w/N; 
C = C * N* N;

%formula for rate
%C = dt*sum( yh.^2 - 2*yh.*y_hist + 2/sqrt(2*pi)/w*y_hist );

    
function y = fftkernel(x,w)
% y = fftkernel(x,w)
%
% Function `fftkernel' applies the Gauss kernel smoother to 
% an input signal using FFT algorithm.
%
% Input argument
% x:    Sample signal vector. 
% w: 	Kernel bandwidth (the standard deviation) in unit of 
%       the sampling resolution of x. 
%
% Output argument
% y: 	Smoothed signal.
%
% MAY 5/23, 2012 Author Hideaki Shimazaki
% RIKEN Brain Science Insitute
% http://2000.jukuin.keio.ac.jp/shimazaki

L = length(x);
Lmax = max(1:L+3*w);
n = 2^(nextpow2(Lmax));

X = fft(x,n);

f = (0:n-1)/n;
f = [-f(1:n/2+1) f(n/2:-1:2)];

K = exp(-0.5*(w*2*pi*f).^2);

y = ifft(X.*K,n);

y = y(1:L);


function y = logexp(x) 
if x<1e2 
    y = log(1+exp(x));
else
    y = x;
end

function y = ilogexp(x)
%ilogexp = @(x) log(exp(x)-1);
if x<1e2
    y = log(exp(x)-1);
else
    y = x;
end
    


