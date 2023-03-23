function [x,f,gNorm] = lbfgs(x0, F, varargin)
    % limitet-memory bfgs optimization function
    %
    % Arguments:
    %   x0          initial guess
    %   F           Objective function
    %
    % varargin:
    %   storedVec   Number of vectors used to approximate hessian
    %   maxIt       Maximum number of iterations
    %   tol         Convergence tolerance. Convergence test is
    %               gradF(x) <= tol*gradF(x0)
    %   minStep     if steplength is less than minStep, the function
    %               returns
    %   inDomain    function with default value @(x) true(size(x,1),1).
    %               should return a array of logical where the i'th element
    %               is false if the i'th input point is outside your
    %               domain.
    %   
    % Returns:
    %   x           optimal point
    %   f           function value at each step
    %   gNorm       norm of the gradient at each step

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2016 Runar Lie Berge. See COPYRIGHT.TXT for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

    opt = struct('storedVec', 5,    ...
                 'maxIt',     1000, ...
                 'tol',       1e-6, ...
                 'minStep',   10*eps,...
                 'inDomain',  @(x) true(size(x,1),1));

    opt    = merge_options(opt,varargin{:});
    m      = opt.storedVec;
    maxIt  = opt.maxIt;
    tol    = opt.tol;
    
    x      = x0;
    [f,g]  = F(x0);
    gNorm  = norm(g,2);
    H      = speye(size(x,1));
    s      = zeros(size(x,1),m);
    y      = zeros(size(x,1),m);

    for k = 1:maxIt
        if gNorm(k)<=tol*gNorm(1)
            return
        end
        % find search direction
        p      = - twoLoopRec(g, H, s, y, m);
        % finde step length
        
        alpha  = backTracking(F, x, p, opt.inDomain);
        % take step
        x      = x + alpha*p;
        [f(k+1), gNew] = F(x);
        gNorm(k+1) = norm(g,2);
        % update hessian vectors
        s      = [s(:,2:m), alpha*p];
        y      = [y(:,2:m), gNew - g];
       
        fprintf('%3d: f=%10.3e, |df|=%10.3e, |xk+1-xk| = %10.3e\n', ...
                k, f(k+1), gNorm(k+1), norm(x-x0));
       
        % Test progress
        if norm(x-x0)<opt.minStep
            return
        else
            % Precondition hessian
            gamma  = s(:,m)'*y(:,m)/(y(:,m)'*y(:,m));
            H      = gamma*speye(size(x,1));
        end
        x0     = x;
        g      = gNew;
       
    end
    warning('Did not converge in maximum number of iterations');
end


function [r] = twoLoopRec(g, H, s, y, m)
    alpha = zeros(1,m);
    for i = m:-1:1
        rhoInv = y(:,i)'*s(:,i);
        if rhoInv>0
            alpha(i) = s(:,i)'*g/rhoInv;
            g = g - alpha(i)*y(:,i);
        end
    end
    r = H*g;
    for i = 1:m
        rhoInv = y(:,i)'*s(:,i);
        if rhoInv>0
           beta = y(:,i)'*r/rhoInv;
           r = r + s(:,i)*(alpha(i) - beta);
        end
    end
end


function [alpha] = backTracking(F, x, p, inDomain)
    alpha = 1; rho = 0.5; c = 1e-3;
    [f0, g0] = F(x);
    while any(~inDomain(x + p*alpha))&& alpha>1e-16
         alpha = rho*alpha;
    end
    beta = f0 + c*alpha*g0'*p;
    [f1,~] = F(x+alpha*p);
    while f1 > beta && alpha>1e-16
        alpha  = rho*alpha;
        [f1,~] = F(x+alpha*p);
    end
end
