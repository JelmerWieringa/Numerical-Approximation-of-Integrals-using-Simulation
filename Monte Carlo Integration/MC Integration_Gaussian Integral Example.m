clear

%% Monte Carlo Integration Example: Gaussian Integral
rng('default')   % Control random number generator

tic
n = 1000;   % Sample size
% Generate n realizations from N(0,1)
StandardNormal_Realizations = random('Normal',0,1,[n 1]);   

Exponential_Quadratic = exp(-0.5*StandardNormal_Realizations.^2);
Mean_Value = mean(Exponential_Quadratic);   % Expecation
MC_Integral_Estimate = sqrt(2*pi)*Mean_Value;   % Integral estimate
toc 

MSE = 1/(n*(n-1))*sum((sqrt(2*pi).*Exponential_Quadratic ...
    - sqrt(pi)).^2);


%% Convergence plots
% Plotting the convergence of the MC estimate:
Estimate_Convergence = cumsum(sqrt(2*pi).*Exponential_Quadratic)./(1:n)';
plot(1:n, Estimate_Convergence,'LineWidth',2,'Color','blue')
hold on;
yline(sqrt(pi),'--g','LineWidth',2)
title('Convergence of Estimate','FontSize',34,'Interpreter','latex')
xlabel('n','FontSize',21,'Interpreter','latex')
ylabel('Estimate','FontSize',21,'Interpreter','latex')
legend({'Monte Carlo Estimate', '$$\sqrt{\pi}$$'}, ...
    'FontSize',21,'Location','northeast','Interpreter','latex');
hold off;

% Plotting the error convergence of the estimate:
Sample_Variance = cumsum((sqrt(2*pi).*Exponential_Quadratic ...
    - sqrt(pi)).^2)./(0:n-1)';
MSE_Convergence = Sample_Variance./(1:n)';
plot(1:n,MSE_Convergence,'LineWidth',2,'Color','blue')
hold on;
Orange = [240 100 10]/256;
plot(1:n, 1./sqrt(1:n),'LineWidth',2,'Color',Orange)
title('MSE of Monte Carlo Estimate','FontSize',34,'Interpreter','latex')
xlabel('n','FontSize',21,'Interpreter','latex')
ylabel('MSE','FontSize',21,'Interpreter','latex')
legend({'MSE','$$n^{-1/2}$$'},'FontSize',21,'Location','northeast', ...
    'Interpreter','latex');
hold off;