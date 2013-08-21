% Local Regression and Likelihood, Figure 9.1.
%
% Hardle's Motorcycle accelaration dataset. Just a scatterplot!
%
% Author: Catherine Loader

load mcyc;
figure('Name','fig9_1: motorcycle scatterplot');
plot(time,accel,'.');
xlabel('Time');
ylabel('Acceleration');
