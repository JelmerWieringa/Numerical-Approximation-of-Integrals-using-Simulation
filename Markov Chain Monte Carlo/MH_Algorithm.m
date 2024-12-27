function [X] = MH_Algorithm (n,x_0,Sigma2)
    % A function that accepts inputs and returns output X
    % x_0 = Initial state of Markov chain
    % Sigma2 = Deterministic parameter of proposal PDF p
    X = zeros(1,n);  % Create vector to store the Markov chain
    X(1) = x_0;  % In MATLAB, the first index is defined by a 1
    
    for i = 1:(n-1)  % For X_1,...,X_n. Note: x_0 corresponds to index 1 of X
        U = random('Uniform',0,1);
        % Proposal PDF p(y | x_i)
        p_y_xi = random('Normal',X(i),Sigma2);  
        p_xi_y = random('Normal',p_y_xi,Sigma2);  % p(x_i | y)
        
        Tilde_f_X_xi = 0.3*exp(-0.2*X(i)^2) + 0.7*exp( ...
            -0.2*(X(i)-10)^2);  % f_X(x_i)
        Tilde_f_X_y = 0.3*exp(-0.2*p_y_xi^2) + 0.7*exp( ...
            -0.2*(p_y_xi-10)^2);  % f_X(y)
        
        AcceptanceRate = (Tilde_f_X_y*p_xi_y)/(Tilde_f_X_xi*p_y_xi);  
        
        if U <= AcceptanceRate
            X(i+1) = p_y_xi;
        else
            X(i+1) = X(i);
        end
    end
end