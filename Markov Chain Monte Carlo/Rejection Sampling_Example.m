clear

%% Rejection Sampling Example
x_Values = -4:0.001:4;

%%% Tilde_f
Trigonometric_Part = pi*cos(x_Values).^2.*sin(pi*x_Values).^2 + 1;
Tilde_f = zeros(1,length(x_Values));
for i=1:length(x_Values)
    if abs(x_Values(i)) <= 3
        Tilde_f(i) = exp(-0.5*x_Values(i).^2).*Trigonometric_Part(i);
    else
        Tilde_f(i) = 0;
    end
end

PD_StandardUniform = makedist('Uniform',-3,3);
PDF_StandardUniform = pdf(PD_StandardUniform,x_Values);
c = max(Tilde_f);
C = c*(3-(-3));

%%% Plotting the functions with dotted lines
plot(x_Values,Tilde_f,'LineWidth',2,'Color','green')
hold on;
plot(x_Values,PDF_StandardUniform,'LineWidth',2,'Color','blue')
plot(x_Values,C*PDF_StandardUniform,'LineWidth',2,'Color','magenta')

title(['$\tilde{f}_X(x) = e^{-\frac{1}{2}x^2} \cdot \left[ ' ...
    '\pi \cos^2(x) \cdot \sin^2(\pi x) + 1 \right]$'], ...
    'FontSize',27,'Interpreter','latex')
xlabel('x','FontSize',21,'Interpreter','latex')
Grey = [0.7 0.7 0.7];
plot(-4:-3,[max(Tilde_f) max(Tilde_f)],'--','LineWidth',2,'Color',Grey)
text(-5.1,max(Tilde_f),'$C \cdot p_X(-1.2)$', ...
    'FontSize',21,'Interpreter','latex')
xline(-1.2,'--','LineWidth',2,'Color',Grey)
text(-1.3,-0.08,'-1.2','FontSize',10)
f_Value = Tilde_f((-1.2+4)/0.01);  % f(-1.2)

plot(-1.2,f_Value,'Marker','.','MarkerSize',27, ...
    'LineWidth',2,'Color','red');  % Colour intersection
legend({'$\tilde{f}_X$','Unif$(-3,3)$','Scaled Unif$(-3,3)$'}, ...
    'Location','northeast','FontSize',30,'Interpreter','latex');
hold off;

%%% Selecting points below Tilde_f
n=1000;
Uniform_Points = random('Uniform',-3,3,[n 2]);
ScaledUniform_Points = [Uniform_Points(:,1) (c/3)*abs( ...
    Uniform_Points(:,2))];  % Merge first column and scaled second column

Red_Region = zeros(n,2);
Green_Region = zeros(n,2);
for j=1:n
    x_Coordinate = ScaledUniform_Points(j,1);
    y_Coordinate = ScaledUniform_Points(j,2);
    
    % Find index in x_Values for Tilde_f(x_Coordinate)
    Find_Index = (abs(min(x_Values)) + x_Coordinate) / 0.001; 
    % int16() formats it so that it can be used as an input
    Index = int16(Find_Index)+1;  % Small round-off error
    
    if y_Coordinate <= Tilde_f(Index)
        % Colour the point green
        Green_Region(j,1) = x_Coordinate;
        Green_Region(j,2) = y_Coordinate;
    else
        % Colour the point red
        Red_Region(j,1) = x_Coordinate;
        Red_Region(j,2) = y_Coordinate;
    end
end

%%% Plotting the points below the scaled uniform PDF
plot(x_Values,Tilde_f,'LineWidth',2,'Color','green')
hold on;
plot(x_Values,C*PDF_StandardUniform,'LineWidth',2,'Color','magenta')
scatter(Red_Region(:,1),Red_Region(:,2), ...  % Plot the red points
    'MarkerEdgeColor','red','MarkerFaceColor','red')
scatter(Green_Region(:,1),Green_Region(:,2), ...  % Plot the green points
    'MarkerEdgeColor','green','MarkerFaceColor','green')

title('$\tilde{f}_X \ \ \& \ $ Scaled Uniform Proposal PDF', ...
    'FontSize',27,'Interpreter','latex')
xlabel('x','FontSize',21,'Interpreter','latex')
legend({'$\tilde{f}_X$','$C \cdot$Unif$(-3,3)$'}, ...
    'Location','northeast','FontSize',30,'Interpreter','latex');
hold off;


%% New proposal PDF
n=1000;
PD_StandardNormal = makedist('Normal','mu',0,'sigma',1);
PDF_StandardNormal = pdf(PD_StandardNormal,x_Values);
C_N = 9;  % Constant that ensures that proposal PDF lies above target PDF

plot(x_Values,Tilde_f,'LineWidth',2,'Color','green')
hold on;
plot(x_Values,C_N*PDF_StandardNormal,'LineWidth',2,'Color','magenta')

title(['$\tilde{f}_X(x) = e^{-\frac{1}{2}x^2} \cdot \left[ ' ...
    '\pi \cos^2(x) \cdot \sin^2(\pi x) + 1 \right]$'], ...
    'FontSize',27,'Interpreter','latex')
xlabel('x','FontSize',21,'Interpreter','latex')
legend({'$\tilde{f}_X$','$9 \cdot \mathcal{N}(0,1)$'}, ...
    'Location','northeast','FontSize',30,'Interpreter','latex');
hold off;

%%% Generate points the fill the region below the standard normal curve
rng('default')
StandardNormal_Realizations = random('Normal',0,1,[n 1]);  % x-coordinates
StandardNormal_Points = (1/sqrt(2*pi))*exp( ...  % y-coordinates
    -0.5*StandardNormal_Realizations.^2);
Uniform_Scaling = random('Uniform',0,1,[n 1]);  % Spread points vertically
c = C_N;

ScaledStandardNormal_Points = [
    StandardNormal_Realizations c*Uniform_Scaling.*StandardNormal_Points];
scatter(ScaledStandardNormal_Points(:,1),ScaledStandardNormal_Points(:,2))

%%% Selecting points below Tilde_f
Red_Region = zeros(n,2);
Green_Region = zeros(n,2);
for k=1:n
    x_Coordinate = ScaledStandardNormal_Points(k,1);
    y_Coordinate = ScaledStandardNormal_Points(k,2);
    
    % Find index in x_Values for Tilde_f(x_Coordinate)
    Find_Index = (abs(min(x_Values)) + x_Coordinate) / 0.001; 
    % int16() formats it so that it can be used as an input
    Index = int16(Find_Index)+1;  % Small round-off error
    
    if y_Coordinate <= Tilde_f(Index)
        % Colour the point green
        Green_Region(k,1) = x_Coordinate;
        Green_Region(k,2) = y_Coordinate;
    else
        % Colour the point red
        Red_Region(k,1) = x_Coordinate;
        Red_Region(k,2) = y_Coordinate;
    end
end

%%% Plotting the points below the scaled standard normal PDF
plot(x_Values,Tilde_f,'LineWidth',2,'Color','green')
hold on;
plot(x_Values,C_N*PDF_StandardNormal,'LineWidth',2,'Color','magenta')
scatter(Red_Region(:,1),Red_Region(:,2), ...  % Plot the red points
    'MarkerEdgeColor','red','MarkerFaceColor','red')
scatter(Green_Region(:,1),Green_Region(:,2), ...  % Plot the green points
    'MarkerEdgeColor','green','MarkerFaceColor','green')

title('$\tilde{f}_X \ \ \& \ $ Scaled Normal Proposal PDF', ...
    'FontSize',27,'Interpreter','latex')
xlabel('x','FontSize',21,'Interpreter','latex')
legend({'$\tilde{f}(x)$','$9 \cdot \mathcal{N}(0,1)$'}, ...
    'Location','northeast','FontSize',30,'Interpreter','latex');
hold off;