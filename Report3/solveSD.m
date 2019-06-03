function [x, niter] = solveSD(A, b, maxiter)
     % This function solve A*x = b using steepest descent algorithm.
     % It returns the solution and the convergence information.
     % ARGS :
     %  - A : Positive definite NxN matrix
     %  - b : Right-hand side Nx1 column vector
     % RETURN :
     %  - x :  Nx1 solution vector
     %  - conv : Table showing the value of the residual at each step
    
     % Check if maxiter argument was given. If not, assign default value.
     if nargin < 3
       maxiter = 100000000;
    end
    
    % Initialisation of the parameters of the SD method
    N = size(A);
    x = eye(N(1),1); 
    niter = 1;

    % First initialisation of the varialbles (first iteration)
    r = b - A*x;
    delta = norm(r);%r'*r;
    conv = delta;

    % We continue to iterate until the convergence is reach
    while (delta > 1e-12)&&(maxiter > niter)
        % Here we just apply the well know algorithm
        q = A*r;
        alpha = (r'*r)/(r'*q);
        x = x + alpha*r;
        r = b - A*x;
        delta = norm(r);
        conv = [conv, delta];
        niter = niter + 1;
    end
end

