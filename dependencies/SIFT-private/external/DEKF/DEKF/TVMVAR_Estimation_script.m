clear
clc
close all

%% Simulated time-varying MVAR model 
% Reference of the simulated model:
% [1] M. Winterhalder, B. Schelter, W. Hesse et al., “Comparison of linear
% signal processing techniques to infer directed interactions in multivariate 
% neural systems,” Signal Processing, vol. 85, no. 11, pp. 2137-2160, 2005.
% ---------> Example 4 (pp. 13-14)
% 
% Written by: Amir Omidvarnia

L = 5000;                                     % Number of time points
CH = 3;                                       % Number of channels
y = zeros(CH,L);                              % Simulated data: output of the time-varying MVAR model
p = 2;                                        % Model order

%% Define time-varying MVAR parameters
bb = sinc(linspace(pi/2+pi/4,5*pi,L));
b = (.8*(bb-min(bb))/(max(bb)-min(bb)))-.2;   % Parameter 'b': time-varying influence of channel 2 on channel 1
c = zeros(1,L);                               % Parameter 'c': time-varying influence of channel 3 on channel 1

%% Simulate the model
for n = p+1 : L
    if(n<=L/2)
        c(n) = (n/(L/2));
    else
        c(n) = (L-n)/(L/2);
    end
       
    y(1,n) = 0.59*y(1,n-1) - 0.2*y(1,n-2) + b(n)*y(2,n-1) + c(n)*y(3,n-1) + randn;
    y(2,n) = 1.58*y(2,n-1) - 0.96*y(2,n-2) + randn;
    y(3,n) = 0.6*y(3,n-1)  - 0.91*y(3,n-2) + randn;    
end

%% Time-varying MVAR parameter estimation using Dual Extended Kalman Filter (DEKF)
tic
[A Pe Pa] = DEKF(y,p,0.02,true);                                % Estimated time-varying parameters, A = [A1 A2 ... Ar]
t=toc/L

%% Plot parameters
figure, hold on, grid on
plot(squeeze(A(1,2,:)),'k','linewidth',2), plot(b,'r','linewidth',2)
plot(squeeze(A(1,3,:)),'b','linewidth',2), plot(c,'g','linewidth',2)

legend('Estimated b','b','Estimated c','c')
xlabel('Time (sample)','fontsize',14,'fontweight','bold')
ylabel('MVAR parameters','fontsize',14,'fontweight','bold')
set(gca,'fontsize',14,'fontweight','bold')

