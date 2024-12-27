clear

%% Monte Carlo Integration Example: Profit 
rng('default')   % Control random number generator

tic
n = 1000;   % Sample size
Volume_V = 365*100*100;

% Firstly, generate random samples:
X_Sample = round(100*rand(n,1));
Y_Sample = round(100*rand(n,1));
T_Sample = round(365*rand(n,1));

% Secondly, calulate profit pi for all samples:
Seasonal_Effect = (1/3)*cos((2*pi/365).*T_Sample - pi/6) + 1;
Growth = exp(T_Sample./1000);
Q_Function = (80 - 0.05.*X_Sample.*Seasonal_Effect - 0.08.* ...
    Y_Sample).*Growth;
Inflation = (1 + 0.02/365).^T_Sample;
P_Function = 5.*Inflation - (1/200).*X_Sample - (1/300)*Y_Sample;
C_Function = (2 + 0.015*X_Sample + 0.01*Y_Sample).*(1 + T_Sample./1000);

Pi_Function = (P_Function - C_Function).*Q_Function;

% Finally, calulate the total profit (big pi):
Pi = Volume_V*mean(Pi_Function);
toc