clear

%% O(n^{-1/d})
n_Values = 0:1:10;   % Sequence from 0 to 10 by steps of 1
bigO_d1 = 1./n_Values;
bigO_d2 = 1./sqrt(n_Values); 
bigO_d3 = 1./(n_Values.^(1/3)); 

plot(n_Values,bigO_d1,'blue','LineWidth',2)
hold on;
plot(n_Values,bigO_d2,'green','LineWidth',2)
plot(n_Values,bigO_d3,'red','LineWidth',2)
%title('n^{-1/d}')
xlabel('n','FontSize',21,'Interpreter','latex')
ylabel('Value','FontSize',21,'Interpreter','latex')
legend({'$n^{-1}$','$n^{-1/2}$','$n^{-1/3}$'},'Interpreter','latex', ...
    'Location','northeast','FontSize', 30);
hold off;