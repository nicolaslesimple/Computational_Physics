function [vec, val] = eig_power(A)
    % Power method for computing eigenvalues
    % This function computes the eigenvalue lambda of largest module of 
    % the matrix A and the corresponding eigenvector.
    % A tolerance is set up to specifies the tolerance of the method.
    % A number max of itereations is also specify. 
    % z_0 is used as the initial guess.
    % ARGS :
    %   - A : matrix on which we want to find the eigenvalue lambda of 
    %   largest module
    % RETURN :
    %   - vec : eigenvector corresponding to the largest module with unit
    %   norm.
    %   - val : eigenvector corresponding to the largest module
    
    % Define main parameters 
    [na,ma] = size(A);
    nmax = 10000;
    tol = 10^(-10);
    z0 = ones(na,1);
    
    % Check the size of the input matrix A :
    if na ~= ma
        disp('ERROR:Matrix A should be a square matrix')
        return
    end
    
    % Apply power method :
    q=z0/norm(z0);
    q2=q;
    relres=tol+1;
    iter=0;
    z=A*q;
    while ((relres(end)>=tol) && (iter<=nmax))
        q=z/norm(z); 
        z=A*q;
        val=q'*z;
        vec=q;
        z2=q2'*A; 
        q2=z2/norm(z2); 
        q2=q2';
        y1=q2; 
        costheta=abs(y1'*vec);
        iter=iter+1;
        temp=norm(z-val*q)/costheta;
        relres=[relres; temp];
    end
return


