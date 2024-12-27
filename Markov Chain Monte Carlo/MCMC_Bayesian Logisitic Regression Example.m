clear

%% MCMC Example: Bayesian Logisitic Regression


%% Step 1: Generate Regressor Data
rng('default')
n = 10000;  % Length of dataset considered
x = random('Normal',5,7,[n 1]);
X = [ones(n,1) x];  % Design matrix



%% Step 2: Generate Response Observations
Beta_True = [4.2; 2.1];  % Column vector

Pi = zeros(n,1);  % Store succes probabilities
for i = 1:n
    Pi(i) = ( 1 + exp( -X(i,:)*Beta_True ) )^(-1);  % \pi_i
    % Note: operator '*' is used instead of '.*' to compute the 
    %   multiplication of a row vector times a column vector.
end

y = zeros(n,1);  % Store categorical response observations
for i = 1:n
    y(i) = random('Binomial',1,Pi(i));
    % Note: a binomial distribution with one trial equals a Bernoulli
    %   distribution.
end






%% Step 5: Auditing Generated Samples
Beta_Initials = [[0 0]; [1 1]; [2 2]; [4 2]];  % Matrix with 4 rows
Sigma2 = 3;

tic
Sample_1 = MH_Algorithm_BLR(Beta_Initials(1,:),Sigma2,y,X);
Sample_2 = MH_Algorithm_BLR(Beta_Initials(2,:),Sigma2,y,X);
Sample_3 = MH_Algorithm_BLR(Beta_Initials(3,:),Sigma2,y,X);
Sample_4 = MH_Algorithm_BLR(Beta_Initials(4,:),Sigma2,y,X);
toc  % Elapsed time is 322.9203 seconds.

Index = 2;  % Selects either the \beta_0s or \beta_1s from sample
% Index = 1 corresponds to \beta_0, and Index = 2 to \beta_1

%%% Plotting trace plots for the different samples
plot(1:n,Sample_1(:,Index),'Color','green','LineWidth',2)
hold on;
plot(1:n,Sample_2(:,Index),'Color','blue','LineWidth',2)
plot(1:n,Sample_3(:,Index),'Color','cyan','LineWidth',2)
plot(1:n,Sample_4(:,Index),'Color','red','LineWidth',2)

title(['Trace Plot of $\beta_1$-Markov Chains with $n = $' num2str(n)], ...
    'FontSize',27,'Interpreter','latex')
xlabel('$i$','FontSize',21,'Interpreter','latex')
ylabel('$\beta_{1,i}$','FontSize',21,'Interpreter','latex')
legend({'\mbox{\boldmath$\beta_1$} $= (0,0)^T, \ \sigma^2 = 3$', ...
    '\mbox{\boldmath$\beta_1$} $= (1,1)^T, \ \sigma^2 = 3$', ...
    '\mbox{\boldmath$\beta_1$} $= (2,2)^T, \ \sigma^2 = 3$', ...
    '\mbox{\boldmath$\beta_1$} $= (4,2)^T, \ \sigma^2 = 3$'}, ...
    'Location','northeast','FontSize',24,'Interpreter','latex');
% Note: Adjust the beta subscript when changing 'Index'
hold off;



%%% Gelman-Rubin Convergence Statistic
% i = 1,...n denotes the time index of a single M-H Markov chain
% j = 1,...,m labelling of different chains
n_GR = n;  % n corresponding to the Gelman-Rubin statistic
m = 4;  % Number of different Markov chains considerd

% Store the M-H Markov chains, where each column corresponds to one of 
% the m initial values & the rows correspond to the time index 
VBS_Chains = zeros(2,m);  % Variance between sequences
VWS_Chains = zeros(2,m);  % Variance within sequences

Merged_Samples = [Sample_1(1:n_GR, ...
    :) Sample_2(1:n_GR,:) Sample_3(1:n_GR,:) Sample_4(1:n_GR,:)];
B = [1 2]; R = [1 2];  % To store the relevant statistics

for l = 1:2  % For beta_0 and beta_1
    for j = 1:m 
        Sequence = Merged_Samples(:,j+(l-1));
        VBS_Chains(l,j) = mean(Sequence);
        VWS_Chains(l,j) = var(Sequence);  
    end
    
    % Variance between sequences
    B(l) = n_GR/(m-1) * sum( VBS_Chains(l,:) - mean(VBS_Chains(l,:)) );  
    W = mean(VWS_Chains(l,:));  % Variance within sequences
    R(l) = 1 - 1/n_GR + B(l)/(n_GR*W);  % Gelman-Rubin Convergence Statistic
end





%%% Autocorrelation Plot
% Let us use 'Sample_3' to approximate the integral of interest
Beta_Sample = Sample_3(500:end,:);  % Delte Burn-in period of 500 samples

MaxLag = 100;  % Maximum number of lags
AutoCorrelation = zeros(2,MaxLag);  % Store autocorrelations

for l = 1:2  % For beta_0 and beta_1
    for k = 0:MaxLag
        Cov_Xk = mean( (Beta_Sample(1:end-k,l) - mean( ...
            Beta_Sample(:,l))) .* (Beta_Sample(k+1:end,l) - ...
            mean(Beta_Sample(:,l))) ); 
        AutoCorrelation(l,k+1) = Cov_Xk/var(Beta_Sample(:,l));                   
    end
end

%%% Plotting the autocorrelation function
stem(0:MaxLag,AutoCorrelation(1,:), ...
    'Marker','none','LineWidth',2,'Color','green')
hold on;
stem(0:MaxLag,AutoCorrelation(2,:), ...
    'Marker','none','LineWidth',2,'Color','green')
title(['ACF for $n = $ ' num2str(n) '$, \ $ \mbox{\boldmath$\beta_0$}' ...
    ' $= (2,2)^T \ \& \ \ \sigma^2 =$ ' num2str(Sigma2)], ...
    'FontSize',27,'Interpreter','latex')
xlabel('Lag $k$','FontSize',21,'Interpreter','latex')
ylabel('$\rho_k$','FontSize',21,'Interpreter','latex')
hold off;








%% Step 6: Ergodic Theorem
n_AfterBurnIn = length(Beta_Sample(:,1));

x_new = [1 random('Normal',5,7)];  % x_{new}
Pi_new = zeros(n_AfterBurnIn-1,1);
for i = 1:(n_AfterBurnIn-1)  % \pi_{new,i}
    Pi_new(i) = ( 1 + exp(-x_new*transpose(Beta_Sample(i,:))) )^(-1);  
end

Predictive_Probability_Success = (1/n_AfterBurnIn)*sum(Pi_new);




%% Step 7: Visualizing Predictive Probability
x_new_Values = -20:0.1:20;
PPS_Image = zeros(length(x_new_Values),1);
for k = 1:length(x_new_Values)
    PPS_Image(k) = PPS(x_new_Values(k),Beta_Sample);
end

plot(x_new_Values,PPS_Image,'Color','magenta','LineWidth',2)
hold on;
plot(x_new_Values,( 1 + exp(-x_new_Values) ).^(-1), ...
    'Color','cyan','LineWidth',2)  % Plotting the Sigmoid function

title(['$p(y_{new} = 1 \mid$ \mbox{\boldmath$x$}$_{new},$ ' ...
    '\mbox{\boldmath$y$}, \mbox{\boldmath$X$}$)$'], ...
    'FontSize',27,'Interpreter','latex')
xlabel('$x_{new}$','FontSize',21,'Interpreter','latex')
legend({'Empirical Plot','$t \mapsto (1+e^{-t})^{-1}$'}, ...
    'Location','northeast','FontSize',30,'Interpreter','latex');
hold off;

