function [L, U] = LU_square_matrix(A)
    % LU_square_matrix computes the LU decomposition of a square matrix
    % without pivoting.
    % ARGS :
    %   - Matrix on which we want to apply the decomposition.
    % RETURN :
    %   - lower matrix L, and upper matrix U such that A = L*U
   
    n = size(A, 1);
    L = eye(n); % We first fill up the lower triangular half
    for k = 1 : n
        % For each row k, access columns from k+1 to the end and divide by
        % the diagonal coefficient at A(k ,k)
        L(k+1:n,k) = A(k + 1 : n, k) / A(k, k);
        for l=k+1:n % For each row k+1 to the end, perform Gaussian elimination
            A(l, :) = A(l, :)-L(l, k)*A(k, :); % A is U at the end
        end
    end
    U = A;
end
