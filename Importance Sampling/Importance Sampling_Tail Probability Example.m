clear

%% Importance Sampling Example: Tail Probability
rng('default')   % Control random number generator

%% True Integral Value
mu = 0; sigma = 1;
% We take the following value as the true value for the integral 
True_Integral = 1 - cdf('Normal',pi,mu,sigma);  % Tail probability

%% Plotting the region of interest
x_Values = [-4:0.1:3.1, pi, 3.2:0.1:4];  % Sequence from -5 to 5 by 0.1
% Theoretical Standard Normal PDF
PD_Normal = makedist('Normal','mu',0,'sigma',1);
PDF_Normal = pdf(PD_Normal,x_Values);

plot(x_Values,PDF_Normal,'LineWidth',2,'Color','black')
hold on;
% Colouring the region under the curve after pi
AfterPi = x_Values >= pi;
patch([x_Values(AfterPi) flip(x_Values(AfterPi))], ...
    [zeros(size(PDF_Normal(AfterPi))) flip(PDF_Normal(AfterPi))], ...
    'magenta','EdgeColor','none')
xline(pi,'LineWidth',2,'Color','cyan')
text(3.2,0.15,'\pi','FontSize',27)
title('Region $Z \geq \pi$','fontsize',34,'Interpreter','latex')
xlabel('z','fontsize',21,'Interpreter','latex')
ylabel('PDF','fontsize',21,'Interpreter','latex')
legend({'$\mathcal{N}(0,1)$'},'Interpreter','latex','Location', ...
    'northeast','FontSize', 30);
hold off;

%% Zoom in on the coloured region in the tail
plot(x_Values,PDF_Normal,'LineWidth',2,'Color','black')
hold on;
% Coloring the region under the curve after pi 
AfterPi = x_Values >= pi;
patch([x_Values(AfterPi) flip(x_Values(AfterPi))], ...
    [zeros(size(PDF_Normal(AfterPi))) flip(PDF_Normal(AfterPi))], ...
    'magenta','EdgeColor','none')
xline(pi,'LineWidth',2,'Color','cyan')
xlim([3 4])
text(3.15,0.004,'\pi','FontSize',27)
title('Zoomed-In Region $Z \geq \pi$', ...
    'fontsize',34,'fontweight','bold','Interpreter','latex')
xlabel('z','fontsize',21,'Interpreter','latex')
ylabel('PDF','fontsize',21,'Interpreter','latex')
legend({'$\mathcal{N}(0,1)$'},'Interpreter','latex','Location', ...
    'northeast','fontsize', 30);
hold off;



%% Monte Carlo Integration 
tic
n = 10000;   % Sample size
% Generate n realizations from N(0,1)
StandardNormal_Realization = random('Normal',mu,sigma,[n 1]); 
% MC integral estimate: mean of indicator functions
MC_Integral_Estimate = mean(StandardNormal_Realization>=pi);
toc

% Variance of MC estimator (= MSE, by unbiasedness)
MSE_MCIntegration = 1/(n*(n-1))*sum( ...
    (StandardNormal_Realization>=pi - MC_Integral_Estimate).^2);


%% Importance Sampling: Exponential Proposal PDF
tic
% Generate n realizations from Expo(1)
Exponential_Realization = random('Exponential',1,[n 1]);
% Indicator function: Z >= pi
Select_Exponential = Exponential_Realization( ...
    Exponential_Realization>=pi,1);
% IS integral estimate: mean of exponential quadratic
Transform_Exponential = exp(-0.5*Select_Exponential.^2)./exp( ...
    -Select_Exponential);
IS_Exponential = (1/sqrt(2*pi))*(sum(Transform_Exponential)/n); 
toc

% Variance of IS estimator (= MSE, by unbiasedness)
MSE_Exponential = 1/(n-1)*sum( ...
    ((1/sqrt(2*pi))*Transform_Exponential - IS_Exponential).^2 );


%% Importance Sampling: Truncated Exponential Proposal PDF
tic
% Generate n realizations from Unif(0,1)
U = random('Uniform',0,1,[n 1]); 
% Transform to CDF of truncated exponential
TruncatedExponential_Realization = pi - log(1-U);
% IS integral estimate: mean of exponential quadratic
Transform_Truncated = exp(-0.5*TruncatedExponential_Realization.^2 ...
    )./exp(-TruncatedExponential_Realization + pi);
% Mean function can be used because we do not have any deleted zeros
IS_TruncatedExponential = (1/sqrt(2*pi))*mean(Transform_Truncated);  
toc

% Variance of MC estimator (= MSE, by unbiasedness)
MSE_TruncatedExponential = 1/(n-1)*sum( ...
    ((1/sqrt(2*pi))*Transform_Truncated - IS_TruncatedExponential).^2 );


%% Importance Sampling: Gamma Proposal PDF
tic
% Generate n realizations from Gamma(1,2)
Gamma_Realization = random('Gamma',1,2,[n 1]);
% Indicator function: Z >= pi
Select_Gamma = Gamma_Realization(Gamma_Realization>=pi,1);
% IS integral estimate: mean of exponential quadratic
Transform_Gamma = exp(-0.5*Select_Gamma.^2)./exp(-0.5*Select_Gamma);
IS_Gamma = (sqrt(2)/sqrt(pi))*(sum(Transform_Gamma)/n);
toc

% Variance of MC estimator (= MSE, by unbiasedness)
MSE_Gamma = 1/(n-1)*sum( ...
    ((sqrt(2)/sqrt(pi))*Transform_Gamma - IS_Gamma).^2 );


%% Importance Sampling: Weibull Proposal PDF
tic
alpha = 2; beta = sqrt(2);
% Generate n realizations from Weibull(2,sqrt(2))
Weibull_Realization = random('Weibull',alpha,beta,[n 1]);
% Indicator function: Z >= pi
Select_Weibull = Weibull_Realization(Weibull_Realization>=pi,1);
% IS integral estimate: mean of 1/Z
Transform_Weibull = 1./Select_Weibull;
IS_Weibull = (1/sqrt(2*pi))*(sum(Transform_Weibull)/n);
toc

% Variance of MC estimator (= MSE, by unbiasedness)
MSE_Weibull = 1/(n-1)*sum( ...
    ((1/sqrt(2*pi))*Transform_Weibull - IS_Weibull).^2 );


%% Tails Comparison
x_Values = -5:0.1:10;  % Sequence from -5 to 10 by steps of 0.1

% Theoretical Standard Normal PDF
PD_Normal = makedist('Normal','mu',0,'sigma',1);
PDF_Normal = pdf(PD_Normal,x_Values);
% Theoretical Standard Exponential PDF
PD_Exponential = makedist('Exponential',1);
PDF_Exponential = pdf(PD_Exponential,x_Values);
% Theoretical Standard Truncated Exponential PDF
PD_Exponential_Truncated = truncate(PD_Exponential,pi,inf);
PDF_Exponential_Truncated = pdf(PD_Exponential_Truncated,x_Values);
% Theoretical Standard Weibull PDF
PD_Gamma = makedist('Gamma',1,2);
PDF_Gamma= pdf(PD_Gamma,x_Values);
% Theoretical Standard Weibull PDF
PD_Weibull = makedist('Weibull',2,sqrt(2));
PDF_Weibull= pdf(PD_Weibull,x_Values);

plot(x_Values,PDF_Normal,'LineWidth',2,'Color','black')
hold on;
plot(x_Values,PDF_Exponential,'LineWidth',2,'Color','blue')
plot(x_Values,PDF_Exponential_Truncated,'LineWidth',2,'Color','cyan')
plot(x_Values,PDF_Gamma,'LineWidth',2,'Color','magenta')
plot(x_Values,PDF_Weibull,'LineWidth',2,'Color','red')

title('Comparison of PDFs','FontSize',27,'Interpreter','latex')
xlabel('z','FontSize',21,'Interpreter','latex')
ylabel('PDF','FontSize',21,'Interpreter','latex')
legend({'$\mathcal{N}(0,1)$','Expo$(1)$','$\mathcal{T}$Expo$(1,\pi)$', ...
    'Gamma$(1,2)$','Weibull$(2,\sqrt{2})$'}, ...
    'Location','northeast','FontSize',30,'Interpreter','latex');
hold off;


%% Likelihood Ratio Comparison
plot(x_Values,PDF_Normal./PDF_Exponential,'LineWidth',2,'Color','blue')
hold on;
plot(x_Values,PDF_Normal./PDF_Exponential_Truncated, ...
    'LineWidth',2,'Color','cyan')
plot(x_Values,PDF_Normal./PDF_Gamma,'LineWidth',2,'Color','magenta')
plot(x_Values,PDF_Normal./PDF_Weibull,'LineWidth',2,'Color','red')

title('Likelihood Ratios','FontSize',27,'Interpreter','latex')
xlim([0 5])
xlabel('z','FontSize',21,'Interpreter','latex')
ylabel('PDF','FontSize',21,'Interpreter','latex')
legend({'Expo$(1)$','$\mathcal{T}$Expo$(1,\pi)$', ...
    'Gamma$(1,2)$','Weibull$(2,\sqrt{2})$'}, ...
    'Location','northeast','FontSize',30,'Interpreter','latex');
hold off;








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot over sample size
rng('default')   % Control random number generator

tic
% Sequence of sample size from 1000 to 100,000 by steps of 1000
n_Values = 1000:1000:100000;   
% Next, create empty vectors to store values:
MC_Integral_Estimates = zeros(size(n_Values));
IS_Exponential = zeros(size(n_Values));
IS_TruncatedExponential = zeros(size(n_Values));
IS_Gamma = zeros(size(n_Values));
IS_Weibull = zeros(size(n_Values));

MSE_MCIntegration = zeros(size(n_Values));
MSE_Exponential = zeros(size(n_Values));
MSE_TruncatedExponential = zeros(size(n_Values));
MSE_Gamma = zeros(size(n_Values));
MSE_Weibull = zeros(size(n_Values));

for i = 1:length(n_Values)
    N = n_Values(i);

    %%% MC Integration
    % Generate n realizations from N(0,1)
    StandardNormal_N_Realization = random('Normal',0,1,[N 1]);
    % MC integral estimate: mean of indicator functions
    MC_Integral_Estimates(i) = mean(StandardNormal_N_Realization>=pi);
    % Variance (=MSE) of estimator
    MSE_MCIntegration(i) = 1/(N-1)*sum( ...
        (StandardNormal_N_Realization>=pi - ...
        MC_Integral_Estimates(i)).^2 );

    %%% IS: Exponential 
    % Generate n realizations from Expo(1)
    Exponential_N_Realization = random('Exponential',1,[N 1]);
    % Indicator function: Z >= pi
    Select_N_Exponential = Exponential_N_Realization( ...
        Exponential_N_Realization>=pi,1);
    % IS integral estimate: mean of exponential quadratic
    Transform_N_Exponential = exp(-0.5*Select_N_Exponential.^2 + ...
        Select_N_Exponential);
    IS_Exponential(i) = (1/sqrt(2*pi))*(sum(Transform_N_Exponential)/N); 
    % Variance (=MSE) of estimator
    MSE_Exponential(i) = 1/(N-1)*sum( ...
        ((1/sqrt(2*pi))*Transform_N_Exponential - ...
        IS_Exponential(i)).^2 );

    %%% IS: Truncated Exponential
    % Generate n realizations from Unif(0,1)
    U = random('Uniform',0,1,[N 1]); 
    % Transform to CDF of truncated exponential
    TruncatedExponential_N_Realization = pi - log(1-U);
    % IS integral estimate: mean of exponential quadratic
    Transform_N_Truncated = exp(-0.5* ...
        TruncatedExponential_N_Realization.^2 ...
        + TruncatedExponential_N_Realization - pi);
    % Mean function can be used because we do not have any deleted zeros
    IS_TruncatedExponential(i) = (1/sqrt(2*pi))*mean( ...
        Transform_N_Truncated);  
    % Variance (=MSE) of estimator
    MSE_TruncatedExponential(i) = 1/(N-1)*sum( ((1/sqrt(2*pi))* ...
        Transform_N_Truncated - IS_TruncatedExponential(i)).^2 );
    
    %%% IS: Gamma
    % Generate n realizations from Gamma(1,2)
    Gamma_N_Realization = random('Gamma',1,2,[N 1]);
    % Indicator function: Z >= pi
    Select_N_Gamma = Gamma_N_Realization(Gamma_N_Realization>=pi,1);
    % IS integral estimate: mean of exponential quadratic
    Transform_N_Gamma = exp(-0.5*Select_N_Gamma.^2 + 0.5*Select_N_Gamma);
    IS_Gamma(i) = (sqrt(2)/sqrt(pi))*(sum(Transform_N_Gamma)/N);
    % Variance (=MSE) of estimator
    MSE_Gamma(i) = 1/(N-1)*sum( ...
        ((sqrt(2)/sqrt(pi))*Transform_N_Gamma - IS_Gamma(i)).^2 );

    %%% IS: Weibull
    alpha = 2; beta = sqrt(2);
    % Generate n realizations from Weibull(2,sqrt(2))
    Weibull_N_Realization = random('Weibull',alpha,beta,[N 1]);
    % Indicator function: Z >= pi
    Select_N_Weibull = Weibull_N_Realization( ...
        Weibull_N_Realization>=pi,1);
    % IS integral estimate: mean of 1/Z
    Transform_N_Weibull = 1./Select_N_Weibull;
    IS_Weibull(i) = (1/sqrt(2*pi))*(sum(Transform_N_Weibull)/N);
    % Variance (=MSE) of estimator
    MSE_Weibull(i) = 1/(N-1)*sum( ...
        ((1/sqrt(2*pi))*Transform_N_Weibull - IS_Weibull(i)).^2 );
end
toc


%%% Plot the estimates
plot(n_Values,MC_Integral_Estimates,'LineWidth',2,'Color','yellow')
hold on;
plot(n_Values,IS_Exponential,'LineWidth',2,'Color','cyan')
plot(n_Values,IS_TruncatedExponential,'LineWidth',2,'Color','magenta')
plot(n_Values,IS_Gamma,'LineWidth',2,'Color','green')
plot(n_Values,IS_Weibull,'LineWidth',2,'Color','blue')

title('Convergence of Estimators','FontSize',27,'Interpreter','latex')
yline(True_Integral,'LineWidth',2);
xlabel('Sample Size $n$','FontSize',21,'Interpreter','latex')
ylabel('Probability','FontSize',21,'Interpreter','latex')
legend({'MC Integration','IS $\sim$ Expo$(1)$', ...
    'IS $\sim$ $\mathcal{T}$Expo$(1,\pi)$','IS $\sim$ Gamma$(1,2)$', ...
    'IS $\sim$ Weibull$(2,\sqrt{2})$'}, ...
    'Location','northeast','FontSize',30,'Interpreter','latex');
hold off;



%%% Plot the estimates: Without Weibull
plot(n_Values,MC_Integral_Estimates,'LineWidth',2,'Color','yellow')
hold on;
plot(n_Values,IS_Exponential,'LineWidth',2,'Color','cyan')
plot(n_Values,IS_TruncatedExponential,'LineWidth',2,'Color','magenta')
plot(n_Values,IS_Gamma,'LineWidth',2,'Color','green')

title('Convergence of Estimators $\backslash$Weibull', ...
    'FontSize',27,'Interpreter','latex')
yline(True_Integral,'LineWidth',2);
xlabel('Sample Size $n$','FontSize',21,'Interpreter','latex')
ylabel('Probability','FontSize',21,'Interpreter','latex')
legend({'MC Integration','IS $\sim$ Expo$(1)$', ...
    'IS $\sim$ $\mathcal{T}$Expo$(1,\pi)$','IS $\sim$ Gamma$(1,2)$'}, ...
    'Location','northeast','FontSize',30,'Interpreter','latex');
hold off;


%%% MSE plot
% Apply a logarithmic transformation to improve comparison
plot(n_Values,log(MSE_MCIntegration),'LineWidth',2,'Color','yellow')
hold on;
plot(n_Values,log(MSE_Exponential),'LineWidth',2,'Color','cyan')
plot(n_Values,log(MSE_TruncatedExponential), ...
    'LineWidth',2,'Color','magenta')
plot(n_Values,log(MSE_Gamma),'LineWidth',2,'Color','green')
plot(n_Values,log(MSE_Weibull),'LineWidth',2,'Color','blue')
Orange = [240 100 10]/256;
plot(n_Values,log(1./sqrt(n_Values)),'LineWidth',2,'Color',Orange)

title('$\log$(MSE) Comparison','fontsize',27,'Interpreter','latex')
xlabel('Sample Size $n$','FontSize',21,'Interpreter','latex')
ylabel('$\log$(MSE)','FontSize',21,'Interpreter','latex')
legend({'MC Integration','IS $\sim$ Expo$(1)$', ...
    'IS $\sim$ $\mathcal{T}$Expo$(1,\pi)$','IS $\sim$ Gamma$(1,2)$', ...
    'IS $\sim$ Weibull$(2,\sqrt{2})$','$\mathcal{O}(n^{-1/2})$'}, ...
    'Location','northeast','FontSize',30,'Interpreter','latex');
hold off;