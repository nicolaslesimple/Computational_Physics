function [vec,val_k]=eig_ipower(A,mu)
    % Inverse power method with shift.
    % This function computes the eigenvalue lambda of smallest module of 
    % the matrix A and the corresponding eigenvector.
    % A tolerance is set up to specifies the tolerance of the method.
    % A number max of itereations is also specify. 
    % x_0 is used as the initial guess.
    % ARGS :
    %   - A : matrix on which we want to find the eigenvalue lambda of 
    %   largest module
    %   - mu : It correspond to the shift.
    % RETURN :
    %   - vec : eigenvector corresponding to the smallest module with unit
    %   norm.
    %   - val : eigenvector corresponding to the smallest module.
    
    % Define the parameters and the initial values
    epsilon=1e-15;
    vec0=rand(length(A),1);
    
    % First definition of the return eigenvector and eigenvalue
    vec=vec0/(norm(vec0));
    val=vec'*A*vec;
    vec=A*vec;
    vec=vec/(norm(vec));
    val_k=vec'*A*vec;

    % Iterate and change the return values since the difference is higher
    % than epsilon
    while ((abs(val-val_k)) > epsilon)
        vec=solve_ls((A-mu*eye(length(A))),vec);
        vec=vec/(norm(vec));
        val=val_k;
        val_k=vec'*A*vec;
    end

end

