function [vec, val] = eig_rq(input_matrix, target)
    % eig_rq computes the closest eigen vector and eigen value
    % of a given matrix with Rayleigh quotient iteration
    %
    % ARGS :
    %   - input matrix (2D complex Hermitian matrix) : matrix
    %   for the eigen value problem;
    %   - target (real scalar) : an estimation to the eigen value;
    % RETURN : 
    %   - a right eigen vector and the corresponding eigen value
    %   of a matrix.

    % Define main parameters 
    [m,n] = size(input_matrix);
    nmax = 10000;
    tol = 10^(-10);
    x0 = ones(n,1);
    sigma = target;
    
    if m~=n
          disp('matrix input_matrix  is not square')  ;
          return;
    end
    
    vec = x0;
    for k = 0 : nmax
         val = (vec'*input_matrix*vec)/(vec'*vec);
         xhat = (input_matrix-sigma * eye(n,n))\vec;
         vec = xhat/max(xhat);
         if  norm( (input_matrix-val * eye(n,n))*vec )  < tol
             return;
         end
    end
end