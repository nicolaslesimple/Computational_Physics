function A=puzzleA(N)
%%%%
%%%%
    T=full(gallery('tridiag',N,1,1,1));
    A=kron(T,eye(N))+kron(eye(N),T)-kron(eye(N),eye(N));
end