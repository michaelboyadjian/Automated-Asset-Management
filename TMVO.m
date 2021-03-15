function  x = TMVO(mu, Q, x0)
    
    % Use this function to construct a TRADEOFF MVO portfolio.
    % *************** WRITE YOUR CODE HERE ***************
    %----------------------------------------------------------------------
    
    % Find the total number of assets
    n = size(Q,1); 
    
    % Set the target as the average expected return of all assets
    targetRet = mean(mu);
    
    % Disallow short sales
    lb = zeros(n,1);
    ub = ones(n,1)*0.2;

    % Add the expected return constraint
    A = -1 .* mu';
    b = -1 * targetRet;

    %constrain weights to sum to 1
    Aeq = ones(1,n);
    beq = 1;

    %LAMBA FOR TRADEOFF
    L = 0.5;
    
    % Set the quadprog options 
    options = optimoptions( 'quadprog', 'TolFun', 1e-9, 'Display','off');
    
    % Optimal asset weights
    x = quadprog( 2 * Q, -L*mu, [], [], Aeq, beq, lb, ub, [], options);
    
    %----------------------------------------------------------------------
    
end