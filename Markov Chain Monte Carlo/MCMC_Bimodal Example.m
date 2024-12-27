clear

%% Bimodal PDF with Normal Proposal PDF
rng('default')
n = 100000;
x_Values = -10:0.1:20;
Tilde_f_X = 0.3*exp(-0.2*x_Values.^2) + 0.7*exp(-0.2*(x_Values-10).^2);
% Calculating Normalization Constant using Monte Carlo integration
Uniform_Realization = random('Uniform',-5,15, [1 n]);
Tilde_f_X_Uniform = 0.3*exp(-0.2*Uniform_Realization.^2) + 0.7*exp( ...
    -0.2*(Uniform_Realization-10).^2); 
NormalizationConstant = mean(Tilde_f_X_Uniform);  

%%% Plotting Tilde_f_X and f_X
plot(x_Values,Tilde_f_X,'Color','magenta','LineWidth',2);  
hold on;
plot(x_Values,NormalizationConstant*Tilde_f_X, ...
    'Color','green','LineWidth',2);  

title('Graphs of $f_X$ and $\tilde{f}_X$', ...
    'FontSize',34,'Interpreter','latex')
xlabel('x','FontSize',21,'Interpreter','latex')
legend({'$\tilde{f}_X$','$f_X$'}, ...
    'Location','northeast','FontSize',30,'Interpreter','latex');
hold off;


%% Visual Explanation of Metropolis-Hastings Algorithm Mechanism
% For example, suppose x_i = 14 and p(y | x_i) = N(x_i, 1)
x_i = 14;
PD_Normal = makedist('Normal','mu',x_t,'sigma',1);
PDF_Normal = pdf(PD_Normal,x_Values);

%%% Plotting Tilde_f_X and the proposal PDF p
plot(x_Values,Tilde_f_X,'Color','magenta','LineWidth',2);  
hold on;
plot(x_Values,PDF_Normal,'LineWidth',2,'Color','blue')
Grey = [0.7 0.7 0.7];
plot([x_i x_i],[0 max(PDF_Normal)],'--','LineWidth',2,'Color',Grey);
text(x_i-0.25,-0.02,'$x_t$','FontSize',21,'Interpreter','latex')
scatter(x_i,max(PDF_Normal),'LineWidth',2, ...
    'MarkerEdgeColor',Grey,'MarkerFaceColor',Grey)

% Proposal Y, 2 cases: 12 and 15.5
y = 15.5;
plot([y y],[0 PDF_Normal(256)],'--','LineWidth',2,'Color','green');
text(y-0.2,-0.02,'$y$','FontSize',21,'Interpreter','latex')
scatter(y,PDF_Normal(256),'LineWidth',2, ...
    'MarkerEdgeColor','green','MarkerFaceColor','green')
% y = 13 corresponds to x_Values' index (13 - -10)*10 + 1 = 231
% y = 15 corresponds to x_Values' index (15.5 - -10)*10 + 1 = 256
title(['Proposal Mechanism with $y= $ ' num2str(y)], ...
    'FontSize',34,'Interpreter','latex')
xlabel('x','FontSize',21,'Interpreter','latex')
legend({'$\tilde{f}_X$','$p(y \mid x_i) = \mathcal{N}(x_i,1)$'}, ...
    'Location','northwest','FontSize',30,'Interpreter','latex');
hold off;


%% Metropolis-Hastings Algorithm
% The MetropolisHastingsAlgorithm() function is defined in a separate
% function file.

%%% Plotting a Metropolis-Hastings Markov Chain
rng('default')
n = 100000; x_0 = 5; Sigma2 = 1;
MC = MH_Algorithm(n,x_0,Sigma2);

plot(1:length(MC),MC,'Color','blue')
hold on;
xlim tight  % Ensure the correct x-axis
title(['M-H Markov Chain with $n = $ ' num2str(n) '$, \ x_0 = $ ' ...
    num2str(x_0) '$\ \ \& \ \ p(y \mid x_i) = \mathcal{N}(x_i,$' ...
    num2str(Sigma2) '$)$'],'FontSize',27,'Interpreter','latex')
xlabel('Iteration $i$','FontSize',21,'Interpreter','latex')
ylabel('$X_i$','FontSize',21,'Interpreter','latex')
hold off;


%%% Plotting the histogram and theoretical PDF
histogram(MC,'Normalization','pdf','FaceColor','cyan')
hold on;
plot(x_Values,NormalizationConstant*Tilde_f_X, ...
    'LineWidth',2,'Color','magenta')

title(['M-H: $\tilde{f}_X$ with $n = $ ' num2str(n) '$, \ x_0 = $ ' ...
    num2str(x_0) '$\ \ \& \ \ p(y \mid x_i) = \mathcal{N}(x_i,$' ...
    num2str(Sigma2) '$)$'],'FontSize',27,'Interpreter','latex')
xlabel('$x$','FontSize',21,'Interpreter','latex')
legend({'MCMC $\tilde{f}_X$','True PDF'}, ...
    'Location','northeast','FontSize',30,'Interpreter','latex');
hold off;


%% Trace Plot of Multiple Metropolis-Hastings Markov Chains
TracePlot_Length = 1000;
rng('default')

plot(1:TracePlot_Length+1,MH_Algorithm(TracePlot_Length,5,10), ...
    'Color','green','LineWidth',1.5)
hold on;
plot(1:TracePlot_Length+1,MH_Algorithm(TracePlot_Length,42,10), ...
    'Color','blue','LineWidth',1.5)
plot(1:TracePlot_Length+1,MH_Algorithm(TracePlot_Length,100,10), ...
    'Color','cyan','LineWidth',1.5)
plot(1:TracePlot_Length+1,MH_Algorithm(TracePlot_Length,121,10), ...
    'Color','red','LineWidth',1.5)
Orange = [240 100 10]/256;
plot(1:TracePlot_Length+1,MH_Algorithm(TracePlot_Length,121,21), ...
    'Color',Orange,'LineWidth',1.5)
xlim tight  % Ensure the correct x-axis
title(['Trace Plot of Markov Chains with $n = $ ' ...
    num2str(TracePlot_Length)],'FontSize',27,'Interpreter','latex')
xlabel('Iteration $i$','FontSize',21,'Interpreter','latex')
ylabel('$X_i$','FontSize',21,'Interpreter','latex')
legend({'$x_0 = 5, \ \sigma^2 = 10$','$x_0 = 42, \ \sigma^2 = 10$', ...
    '$x_0 = 100, \ \sigma^2 = 10$','$x_0 = 121, \ \sigma^2 = 10$', ...
    '$x_0 = 121, \ \sigma^2 = 21$'}, ...
    'Location','northeast','FontSize',30,'Interpreter','latex');
hold off;


%% Gelman-Rubin Convergence Statistic
% i = 1,...n denotes the time index of a single M-H Markov chain
% j = 1,...,m labelling of different chains
rng('default')
n = 100;  % Number of iterations
x_0_mValues = [5 42 100 121];  % Initial state values
m = length(x_0_mValues);  % Number of different Markov chains considerd
Sigma2 = 10;

% Store the M-H Markov chains, where each column corresponds to one of 
% the m initial values & the rows correspond to the time index 
VBS_Chains = zeros(1,m);  % Variance between sequences
VWS_Chains = zeros(1,m);  % Variance within sequences

for j = 1:m 
    Exp_Sequence = mean( MH_Algorithm(n,x_0_mValues(j),Sigma2) );
    VBS_Chains(:,j) = Exp_Sequence;  % 1/n * sum_{i=1}^n X_{i,j}
    
    Var_Sequence = var( MH_Algorithm(n,x_0_mValues(j),Sigma2) );
    % 1/(n-1) * sum_{i=1}^n (X_{i,j} - mean(X_i))^2
    VWS_Chains(:,j) = Var_Sequence;  
end

% Variance between sequences
B = n/(m-1) * sum(VBS_Chains - mean(VBS_Chains));  
W = mean(VWS_Chains);  % Variance within sequences
R = 1 - 1/n + B/(n*W);  % Gelman-Rubin Convergence Statistic


%% Autocorrelation Plot
rng('default')
n = 100000; x_0 = 42; Sigma2 = 10;
MC = MH_Algorithm(n,x_0,Sigma2);

MaxLag = 100;  % Maximum number of lags
% Create vector to store autocorrelations
AutoCorrelation = zeros(1,MaxLag);  

for k = 0:MaxLag
    Cov_Xk = mean( (MC(1:end-k) - mean(MC)).*(MC(k+1:end) - mean(MC)) ); 
    AutoCorrelation(k+1) = Cov_Xk/var(MC);                   
end

%%% Plotting the autocorrelation function
stem(0:MaxLag,AutoCorrelation, ...
    'Marker','none','LineWidth',2,'Color','green')
hold on;
title(['ACF for $n = $ ' num2str(n) '$, \ x_0 = $ ' ...
    num2str(x_0) '$\ \ \& \ \ \sigma^2 =$ ' num2str(Sigma2)], ...
    'FontSize',27,'Interpreter','latex')
xlabel('Lag $k$','FontSize',21,'Interpreter','latex')
ylabel('$\rho_k$','FontSize',21,'Interpreter','latex')
hold off;