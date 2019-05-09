function [x] = solve_ls(A,b)
    % solve_ls solves the system of linear equations A*x=b, for a square matrix
    % A and a column vector b, using the LU decomposition and backward 
    % substitution
    % ARGS :
    %   - A : Square matrix representing the system.
    %   - b : Vector of size (n,1) representing the solution of the ls
    %   system.
    % RETURN :
    %   - x : Vector representing the solution of the ls system Ax = b.

    [L, U] = LU_square_matrix(A);
    y = Forward(L,b);
    x = Backward(U,y);

end

function x = Backward(U,y)
    % Solves the nonsingular upper triangular system  Ux = y.
    % ARGS :
    %   - U : Square matrix representing the system. U is (n,n)
    %   - y : Vector of size (n,1) representing the solution of the ls
    %   system Ly = b. y is (n,1).
    % RETURN :
    %   - x : Vector representing the solution of the ls system Ux = y and
    %   also Ax = b.
    %   x is (n,1).

    n = length(y);
    x = zeros(n,1);
    for j=n:-1:2
        x(j) = y(j)/U(j,j);
        y(1:j-1) = y(1:j-1) - x(j)*U(1:j-1,j);
    end
    x(1) = y(1)/U(1,1);
end

function y = Forward(L,b)
    % Solves the nonsingular lower triangular system  Ly = b 
    % ARGS :
    %   - L : Square matrix representing the system. L is (n,n)
    %   - b : Vector of size (n,1) representing the linear system.
    % RETURN :
    %   - y : Vector representing the solution of the ls system Ly = b.
    %   y is (n,1).

    n = length(b);
    y = zeros(n,1);
    for j=1:n-1
       y(j) = b(j)/L(j,j);
       b(j+1:n) = b(j+1:n) - L(j+1:n,j)*y(j);
    end
    y(n) = b(n)/L(n,n);
end
