function x = Project1_Function(periodReturns, periodFactRet, x0)

    % Use this function to implement your algorithmic asset management
    % strategy. You can modify this function, but you must keep the inputs
    % and outputs consistent.
    %
    % INPUTS: periodReturns, periodFactRet, x0 (current portfolio weights)
    % OUTPUTS: x (optimal portfolio)
    %
    % An example of an MVO implementation with OLS regression is given
    % below. Please be sure to include comments in your code.
    %
    % *************** WRITE YOUR CODE HERE ***************
    %----------------------------------------------------------------------

    % subset the data to consistently use the most recent 1 years
    % for parameter estimation (we experimented between 1 and 3 years)
    returns = periodReturns(end-11:end,:);
    factRet = periodFactRet(end-11:end,:);
    
    % We used a 5 factor FamaFrench model to estimate mu and Q
    [mu, Q] = FAMAOLS(returns, factRet);
    
    % Use Robust MVO that includes a turnover penalty (i.e RTURNMVO) to optimize our portfolio
    x = RTURNMVO(mu, Q,x0);

    %----------------------------------------------------------------------
end
