%% Step 4: Metropolis-Hastings Algorithm for BLR
function [Beta] = MH_Algorithm_BLR (Beta_Initial,Sigma2,y,X)
    % A function that accepts inputs and returns output Beta
    % This is saved in a function file called 'Posterior.m'
    % Beta_Initial = Initial state of Markov chain (row vector)
    % Sigma2 = Deterministic parameter of proposal PDF q
    % y,X = Observed data needed to calculate the posterior
    
    n = length(y);
    Beta = zeros(n,2);  % Create vector to store the Markov chain
    Beta(1,:) = Beta_Initial;  % In MATLAB, the first index is defined by a 1
    
    %AcceptanceRate = zeros(n,1);
    %AR_Numerator = zeros(n,1);
    %AR_Denominator = zeros(n,1); 


    for i = 1:(n-1)  % For Beta_1,...,Beta_n. 
        U = random('Uniform',0,1);
        
        % Generate proposal Z = (z_0, z_1)
        z_0 = random('Normal',Beta(i,1),Sigma2);
        z_1 = random('Normal',Beta(i,2),Sigma2);
        % q(Beta_i | z) / q(z | Beta_i) = 1 by symmetry of normal PDF
        
        AR_Numerator = Posterior(z_0,z_1,y,X);
        AR_Denominator = Posterior(Beta(i,1),Beta(i,2),y,X);
        AcceptanceRate = AR_Numerator/AR_Denominator; 
        
        if U <= AcceptanceRate
            Beta(i+1,:) = [z_0 z_1];
        else
            Beta(i+1,:) = Beta(i,:);
        end
    end
end