%% IMPORTANCE SAMPLING SECTION
%https://nl.mathworks.com/help/stats/prob.normaldistribution.pdf.html
clear

%% Tails Comparison
x_Values = -5:0.1:5;  % Sequence from -5 to 5 by steps of 0.1

% Theoretical Standard Normal PDF
PD_Normal = makedist('Normal','mu',0,'sigma',1);
PDF_Normal = pdf(PD_Normal,x_Values);

% Theoretical Standard Cauchy PDF
PD_Cauchy = makedist('tLocationScale','mu',0,'sigma',1,'nu',1);
PDF_Cauchy = pdf(PD_Cauchy,x_Values);


plot(x_Values,PDF_Normal,'LineWidth',2)
hold on;
plot(x_Values,PDF_Cauchy,'LineWidth',2)

xlabel('x')
ylabel('PDF')
legend({'$\mathcal{N}(0,1)$','Cauchy$(0,1)$'},'Interpreter','latex', ...
    'Location','northeast','fontsize', 30);

hold off;