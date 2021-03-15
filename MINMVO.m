function  x = MINMVO(mu, Q,x0)
    
    %----------------------------------------------------------------------
    %IMPOSED INTEGER RESTRICTIONS SO THAT THERE WAS A MINIMUM ASSET
    %WEIGHTING BUY IN
    %This code (was based off of Matlabs opensource Integer Optimization)
    %We did not use the code too much in the lab to it is not worthwhile
    %documentiong
    % Find the total number of assets
    N = length(mu);
    
    xvars = optimvar('xvars',N,1,'LowerBound',0,'UpperBound',1);
    vvars = optimvar('vvars',N,1,'Type','integer','LowerBound',0,'UpperBound',1);
    zvar = optimvar('zvar',1,'LowerBound',0);
    
    M = 20;
    m = 5;
    qpprob = optimproblem('ObjectiveSense','maximize');
    qpprob.Constraints.mconstr = sum(vvars) <= M;
    qpprob.Constraints.mconstr2 = sum(vvars) >= m;
    
    fmin = 0.02;
    fmax = 0.2;
    
    qpprob.Constraints.fmaxconstr = xvars <= fmax*vvars;
    qpprob.Constraints.fminconstr = fmin*vvars <= xvars;
    
    qpprob.Constraints.allin = sum(xvars) == 1;
    
    lambda = 300;
    
    qpprob.Objective = mu'*xvars - lambda*zvar;
    
    options = optimoptions(@intlinprog,'Display','off'); % Suppress iterative display
    [xLinInt,fval,exitFlagInt,output] = solve(qpprob,'options',options);
    
    thediff = 1e-4;
    iter = 1; % iteration counter
    assets = xLinInt.xvars;
    truequadratic = assets'*Q*assets;
    zslack = xLinInt.zvar;
    
    
    history = [truequadratic,zslack];

    options = optimoptions(options,'LPOptimalityTolerance',1e-10,'RelativeGapTolerance',1e-8,...
                      'ConstraintTolerance',1e-9,'IntegerTolerance',1e-6);
                  
    while abs((zslack - truequadratic)/truequadratic) > thediff % relative error
        constr = 2*assets'*Q*xvars - zvar <= assets'*Q*assets;
        newname = ['iteration',num2str(iter)];
        qpprob.Constraints.(newname) = constr;
        % Solve the problem with the new constraints
        [xLinInt,fval,exitFlagInt,output] = solve(qpprob,'options',options);
        assets = (assets+xLinInt.xvars)/2; % Midway from the previous to the current
    %     assets = xLinInt(xvars); % Use the previous line or this one
        truequadratic = xLinInt.xvars'*Q*xLinInt.xvars;
        zslack = xLinInt.zvar;
        history = [history;truequadratic,zslack];
        iter = iter + 1;
    end
    

    
    % Optimal asset weights
    x = xLinInt.xvars

    
    %----------------------------------------------------------------------
    
end