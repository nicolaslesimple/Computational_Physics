function A=puzzleA(N)
    % This function allows the creation of matrix A defined in the 
    % question of the problem.
    % ARGS :
    %   - N : Interger corresponding to the size of the problem
    % RETURN :
    %   - A : Matrix A used to define the problem as linear equation
    T=full(gallery('tridiag',N,1,1,1));
    A=kron(T,eye(N))+kron(eye(N),T)-kron(eye(N),eye(N));
end