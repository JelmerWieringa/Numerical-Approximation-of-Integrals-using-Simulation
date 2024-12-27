%% Function for Visualization of Probabilities
function[Predictive_Probability_Success] = PPS(x_new_i,Beta_Sample)
    % This is saved in a function file called 'PPS.m'
    x_new = [1 x_new_i];  % x_{new}
    
    N = length(Beta_Sample);
    Pi_new = zeros(N-1,1);
    for j = 1:(N-1)  % \pi_{new,i}
        Pi_new(j) = ( 1 + exp(-x_new*transpose(Beta_Sample(j,:))) )^(-1);  
    end
    
    Predictive_Probability_Success = (1/N)*sum(Pi_new);
end