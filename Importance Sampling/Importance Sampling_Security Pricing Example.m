clear

%% Importance Sampling Example: Security Pricing
%%% Plotting the dividend function
x_Values = -5:0.1:25;
d = sin( (2*pi)/21.*x_Values ) .* cos( (2*pi)/exp(1).*x_Values );

plot(x_Values,d,'LineWidth',2,'Color','magenta')
hold on;
title('$d(x) = \sin(\frac{2\pi}{21}x) \cdot \cos(\frac{2\pi}{e}x)$', ...
    'FontSize',34,'Interpreter','latex')
xlabel('x','FontSize',21,'Interpreter','latex')
ylim([-2 2])
ylabel('d','FontSize',21,'Interpreter','latex')
hold off;

%%% Plotting the trajectory of the dividend random variable
rng('default')   % Control random number generator
% Generate 101 realizations from Gamma(1.2,0.5)
X = random('Gamma',1.2,0.5,[1 101]);
% Realizations of dividend random variable
d_X = abs(sin( (2*pi)/21.*X ) .* cos( (2*pi)/exp(1).*X ));

plot(-5:0.1:5,d_X,'LineWidth',2,'Color','magenta')
hold on;
title(['$d(X) = |\sin(\frac{2\pi}{21}X) \cdot \cos(\frac{2\pi}{e}X)|,' ...
    '$ $X \sim$ Gamma$(1.2,0.5)$'],'FontSize',34,'Interpreter','latex')
xlabel('Time Period (Days)','FontSize',21,'Interpreter','latex')
ylabel('Dividend (in currency units)','FontSize',21,'Interpreter','latex')
hold off;

%%% Computing the Weighted IS Estimate
rng('default')   % Control random number generator
tic
% Sequence of sample size from 1000 to 100,000 by steps of 1000
n_Values = 1000:1000:100000;   
% Next, create empty vectors to store values:
WS_Exponential = zeros(size(n_Values));
SampleVariance_Exponential = zeros(size(n_Values));

for i = 1:length(n_Values)
    N = n_Values(i);

    %%% WS: Exponential 
    % Generate n realizations from Expo(1)
    Exponential_N_Realization = random('Exponential',1,[N 1]);
    % Caculate g(X_i)
    g_N_Exponential = (1/1.05)*abs( ...
        sin( (2*pi)/21.*Exponential_N_Realization ) .* ...
        cos( (2*pi)/exp(1).*Exponential_N_Realization ));
    % Calculate weights
    w_N_Exponential = 2*Exponential_N_Realization.^(0.2);
    % Calculate normalization constants
    NormalizationConstants_N_Exponential = w_N_Exponential/sum( ...
        w_N_Exponential);
    % Calculate the weighted importance sampling estimator
    WS_Exponential(i) = sum(NormalizationConstants_N_Exponential.* ...
        g_N_Exponential); 
    % Sample Variance
    SampleVariance_Exponential(i) = (N/(N-1))*sum( ...
        (NormalizationConstants_N_Exponential.^2) ...
        .*(g_N_Exponential-WS_Exponential(i)).^2 );
end
toc

%%% Convergence plot
plot(n_Values,WS_Exponential,'LineWidth',2,'Color','blue')
hold on;
title('Convergence of $\hat{\theta}_{IS,w}$', ...
    'FontSize',34,'Interpreter','latex')
xlabel('Sample Size $n$','FontSize',21,'Interpreter','latex')
legend('$\hat{\theta}_{IS,w} \sim$ Expo$(\frac{1}{2})$', ...
    'Location','northeast','FontSize',30,'Interpreter','latex');
hold off;

%%% Sample Variance plot
% Apply a logarithmic transformation to improve comparison
plot(n_Values,log(sqrt(SampleVariance_Exponential)), ...
    'LineWidth',2,'Color','blue')
hold on;
Orange = [240 100 10]/256;
plot(n_Values, log(1./sqrt(n_Values)),'LineWidth',2,'Color',Orange)
plot(n_Values, log(1./n_Values),'LineWidth',2,'Color','green')
title('Standard Deviation of $\hat{\theta}_{IS,w}$', ...
    'FontSize',34,'Interpreter','latex')
xlabel('Sample Size $n$','FontSize',21,'Interpreter','latex')
ylabel('$\log(s^2_{IS,w})$','FontSize',21,'Interpreter','latex')
legend({'SD $\hat{\theta}_{IS,w}$', ...
    '$\mathcal{O}(n^{-1/2})$','$\mathcal{O}(n^{-1})$'}, ...
    'Location','northeast','FontSize',30,'Interpreter','latex');
hold off;