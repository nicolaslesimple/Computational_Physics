function [x, niter] = solveCG(A, b, maxiter)
    % This function solve A*x = b using Conjugate Gradients method.
    % It returns the solution and the convergence information.
    % ARGS :
    %  - A : Positive definite NxN matrix     
    %  - b : Right-hand side Nx1 column vector
    % RETURN :
    %  - x :  Nx1 solution vector
    %  - niter : Number of iterations performed

    % Initialisation of the parameters of the SD method
    tol = 1e-12;
    N = size(A);
    s = eye(N(1),1); 
    
    % Check if maxiter argument was given. If not, assign default value.   
    if nargin < 3
        maxiter = 100000000;
    end
   
    % First iteration
    x = s;
    r = b - A*s;
    d = r;         
    rho = r'*r;
    niter = 0;     % Initialize counter

    % We continue to iterate until the convergence is reach
    while (norm(r) > tol)&&(maxiter > niter)
        % Here we just apply the well know algorithm
        a = A*d;
        alpha = rho/(a'*d); % Calculate the step lenght
        x = x + alpha*d;
        r = r - alpha*a;
        rho_new = r'*r;
        d = r + rho_new/rho * d;
        rho = rho_new;
        niter = niter + 1;
    end
end


