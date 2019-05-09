function [val,count] = eig_cj(input_matrix)
    % eig_j computes the eigenvalues of a matrix. 
    % The cyclic Jacobi's method of iteration is used
    % to compute the eigenpairs of a symmetric matrix.
    % ARGS :
    %   - input matrix (2D real symmetric matrix) : matrix
    %   for the eigen value problem.
    % RETURNS :
    %   -val : an array with eigenvalues.
    %   -count  Number of Jacobi rotation.
    
    % Initialization
    epsilon = 10^(-10);
    count = 0;
    D = input_matrix;
    [n,n] = size(input_matrix);
    V = eye(n);
    
    % Select element element of largest magnitude.
    % Use vector representation to be effective.
    [m1,p] = max(abs(D-diag(diag(D)))); 
    [m2,q] = max(m1);                   
    p = p(q);
    
    % Take care of a specific case
    if m2 == 0
        val = diag(D);
        return
    end
    
    % Main part of the module using epsilone in the while condition
    while(off(D) > epsilon)
        for p = 1:(n-1)
          for q = (p+1):n
              count = count+1;
              if D(p,q) ~= 0
                 t = D(p,q)/(D(q,q) - D(p,p));
                 c = 1/sqrt(t*t+1);
                 s = c*t;
              else
                 c=1;
                 s=0;
              end
              R = [c s; -s c];
              D([p q],:) = R'*D([p q],:);
              D(:,[p q]) = D(:,[p q])*R;
              V(:,[p q]) = V(:,[p q])*R;
          end
        end
    end
    D = diag(diag(D));
    val = diag(D);
end