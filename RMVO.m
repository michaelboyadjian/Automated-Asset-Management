function  x = RMVO(mu, Q, x0)
    
    %----------------------------------------------------------------------
    
    % Find the total number of assets
    n = size(Q,1); 
    
    
    % Number of observations;
    T = size(mu, 1);
    
    % Risk aversion parameter
    lambda = 20;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ellipsoidal uncertainty set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Uncertainty set size
    Theta = diag (diag(Q) ) ./ T; 

    % Square root of Theta
    sqrtTh = sqrt(Theta); 

    % Confidence level
    alpha = 0.9;

    % Scaling parameter epsilon for uncertainty set
    ep = sqrt(chi2inv(alpha, n));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup our input parameters for fmincon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% min   lambda * (x' * Q x) - mu' x + epsilon * norm (sqrtTh * x)
% s.t.  sum(x) == 1
%       x >= 0

%-------------------------------------------------------------------------- 
% 3.1 Inequality constraints:
% fmincon accepts linear inequality constraints of the form "A x <= b".
%--------------------------------------------------------------------------

% The only inequality constraint is the one on short selling. However, this
% can be applied as a bound on our variable x.

% Linear inequality Constraint bounds
    b = [];
    A = []; 

% Lower and upper bounds on variables
    %lb = zeros(n,1); 
    lb = [];
    ub = [];

%--------------------------------------------------------------------------
% 3.2 Equality constraints: 
% We only have the budget constraint
%--------------------------------------------------------------------------

    beq = 1; 
    Aeq = ones(1,n);

%--------------------------------------------------------------------------
% 3.3 Initial solution: 
% fmincon requires an initial feasible solution 
%--------------------------------------------------------------------------
% Define an initial portfolio ("equally weighted" or "1/n portfolio")
    x0 = repmat(1.0/n,n,1);

    %x0 = zeros(n, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start 'fmincon' solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Solve using fmincon to find our optimal portfolio weights 
    x = fmincon(@(x)objFun(x, mu, Q, lambda, sqrtTh, ep), x0, A, b, Aeq, beq, lb, ub, @(x)nonlcon(x));
    
 
%--------------------------------------------------------------------------

    function f = objFun(x, mu, Q, lambda, sqrtTh, ep)

        f = (lambda * x' * Q * x) - mu' * x + ep * norm(sqrtTh * x); 

    end

%-------------------------------------------------------------------------- 
% 5.2 Define the equality and inequality nonlinear constraints:
% fmincon accepts nonlinear constraints, but these must be defined as
% separate functions. In our case, we do not have nonlinear constraints. 
%--------------------------------------------------------------------------
    function [c,ceq] = nonlcon(x)

        c = [];
        ceq = [];

    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program End
    
 
    
end