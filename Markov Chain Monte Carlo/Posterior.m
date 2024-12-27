%% Step 3: Posterior Function
function [Output] = Posterior (Beta_0,Beta_1,y,X)
    % A function that accepts inputs and returns Output 
    % This is saved in a function file called 'Posterior.m'
    Prior_Distribution = exp( (-Beta_0^2 - Beta_1^2)/50 );  % p(Beta)
    
    n = length(y);
    Marginal_Likelihood = zeros(1,n);  % p(y | Beta, X)

    Pi = zeros(n,1);  % Store succes probabilities
    for i = 1:n
        Pi(i) = ( 1 + exp( -X(i,:)*[Beta_0; Beta_1] ) )^(-1);  % \pi_i
    end

    % Calculating the marginal likelihoods for all i
    for i = 1:n
        Marginal_Likelihood(i) = Pi(i)^(y(i)) * (1 - Pi(i))^(1-y(i));
    end
    
    Output = prod(Marginal_Likelihood) * Prior_Distribution;
end