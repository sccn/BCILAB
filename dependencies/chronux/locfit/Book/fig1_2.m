% Local Regression and Likelihood, Figure 1.2.
%
% Spencer's mortality data, 15- and 21-point rules.
% Doesn't use locfit!

load spencer;

figure('Name','fig1_2a: Spencer mortality data');
plot(age,mortality,'o');
title('15-point rule');
xlabel('Age (Years)');
ylabel('Mortality Rate');
hold on;
plot(age,spence15(mortality));
hold off;
figure('Name','fig1_2b: Spencer mortality data');
plot(age,mortality,'o');
xlabel('Age (Years)');
ylabel('Mortality Rate');
title('21-point rule');
hold on;
plot(age,spence21(mortality));
hold off;
