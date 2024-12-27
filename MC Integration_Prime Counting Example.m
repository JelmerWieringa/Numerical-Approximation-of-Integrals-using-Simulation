clear

%% Monte Carlo Integration Example: Prime Counting
rng('default')   % Control random number generator

tic
n_Values = 100:100:10000;   % Sequence from 100 to 10,000 by steps of 100
n = 1000;   % Sample size
% Next, create empty vectors to store values:
MC_Log_Integral = zeros(size(n_Values));
Prime_Counts = zeros(size(n_Values));
% Remark: Using length() instead of size() creates a matrix

for i = 1:length(n_Values)
    N = n_Values(i);
    % Generate n realizations from Unif(2,N)
    Uniform_N_Realization = random('Uniform',2,N,[n 1]);
    Log_Transformation = 1./log(Uniform_N_Realization);
    Mean_Value = mean(Log_Transformation);
    MC_Log_Integral(i) = (N-2)*Mean_Value;
    Prime_Counts(i) = length(primes(N));
end
toc


%% Plotting the Monte Carlo Integration estimates
plot(n_Values,MC_Log_Integral,'-og','MarkerFaceColor','green')   
hold on;
Light_Blue = [0 0.5 1];   % Plotting true prime counts
stem(n_Values,Prime_Counts,'Marker','none','LineWidth',1.5, ...
    'Color',Light_Blue)

title(['Monte Carlo Estimation of Logarithmic Integral ' ...
    'for Prime Counting'],'FontSize',34,'Interpreter','latex')
xlabel('$n$','FontSize',21,'Interpreter','latex')
ylabel('$\#$ primes up to integer n','FontSize',21,'Interpreter','latex')
legend({'MC Log Integral','Prime Counts $\pi$'}, ...
    'Location','northwest','FontSize',30,'Interpreter','latex');
hold off;

