clear

%% Monte Carlo Integration Example: Identity Function
rng('default')   % Control random number generator

tic
n = 1000;   % Sample size
% Generate n realizations from Unif(0,1)
StandardUniform_Realization = random('Uniform',0,1,[n 1]);
Mean_Value = mean(StandardUniform_Realization);   % Empirical average
MC_Integral_Estimate = Mean_Value;   % Integral estimate
toc

MSE = 1/(n*(n-1))*sum((StandardUniform_Realization - 0.5).^2);


%% Creating a plot
plot(StandardUniform_Realization,'ro','MarkerFaceColor','red')
% 'ro' specificies red circles without connecting lines
hold on;
% Selecting points below the (0,0)-(1,1) diagonal
Below_Diagonal = StandardUniform_Realization < (1:n)'/n; 
% Colour points below diagonal as green filled circles ('go')
plot(find(Below_Diagonal),StandardUniform_Realization(Below_Diagonal), ...
    'go','MarkerFaceColor','green')   
Grey = [0.5 0.5 0.5];  % [0.5 0.5 0.5] yields a grey color
plot([0 n], [0 1],'Color',Grey,'LineWidth',2)
title('Monte Carlo Integration $\phi(x)=x$', ...
    'FontSize',27,'Interpreter','latex');
xlabel('n','FontSize',21,'Interpreter','latex');
ylabel('Unif$(0,1)$ Realization','FontSize',21,'Interpreter','latex');
hold off;


%% Convergence plots
Estimate_Convergence = cumsum(StandardUniform_Realization)./(1:n)';
plot(1:n, Estimate_Convergence)

Sample_Variance = cumsum((StandardUniform_Realization - 0.5).^2)./(0:n-1)';
MSE_Convergence = Sample_Variance./(1:n)';

plot(1:n,log(MSE_Convergence))
hold on;
plot(1:n,log(1./sqrt(1:n)))
hold off;