clear

%% Markov Chains Example: Jumping Frog
Lambda = [0.5 0.5];
P = [0.6 0.4; 0.8 0.2];

%%% Calculate the distrubtion mu_n for multiple n
tic
n = 20;  % Time steps
% Store the distribution vectors per column in a matrix
Mu = zeros(2,n);
for i = 1:n
    Mu(1:2,i) = Lambda*(P^i);  % Do not use '.*'
end
toc

%%% Plotting the distrubtion mu_n for multiple n
% Probability that the frog sits on plant 1
plot(0:n,[Lambda(1) Mu(1,:)],'-bo','MarkerFaceColor','magenta', ...
    'LineWidth',1,'MarkerSize',10)
hold on;
title('I\kern-0.15em P$(X_n = 1)$','FontSize',34,'Interpreter','latex')
xlabel('Time Period $n$','FontSize',21,'Interpreter','latex')
hold off;

% Probability that the frog sits on plant 1
plot(0:n,[Lambda(2) Mu(2,:)],'-bo','MarkerFaceColor','magenta', ...
    'LineWidth',1,'MarkerSize',10)
hold on;
title('I\kern-0.15em P$(X_n = 2)$','FontSize',34,'Interpreter','latex')
xlabel('Time Period $n$','FontSize',21,'Interpreter','latex')
hold off;




