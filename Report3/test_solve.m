clear all;
close all;
clc;

%This script tests solve_SD.m and solve_CG.m

fprintf('Test 1: Identity matrix for SD ...');
for n=1:1:20
    A=eye(n);
    X=rand(n,1);
    B=A*X;
    sol = solveSD(A,B);
    assert(all(abs(sol - X)./X<1e-10));
end
fprintf('\tpassed\n');

fprintf('Test 2: Identity matrix for CG ...');
for n=1:1:20
    A=eye(n);
    X=rand(n,1);
    B=A*X;
    sol = solveCG(A,B);
    assert(all(abs(sol - X)./X<1e-10));
end
fprintf('\tpassed\n');

fprintf('Test 3: Random real matrices for SD ...')
for n=1:1:20
   A = rand(n,n); % generate a random n x n matrix
   A = A*A';% construct a symmetric matrix using either
   A = A + n*eye(n);
   X=rand(n,1);
   B=A*X;
   sol = solveSD(A,B);
   assert(all(abs(sol - X)./X<1e-10));
      
end
fprintf('\tpassed\n');
   
fprintf('Test 4: Random real matrices CG ...')
for n=1:1:10
   A = rand(n,n); % generate a random n x n matrix
   A = A*A';% construct a symmetric matrix using either
   A = A + n*eye(n);
   X=rand(n,1);
   B=A*X;
   sol = solveCG(A,B);
   assert(all(abs(sol - X)./X<1e-10));
      
end
fprintf('\tpassed\n');

fprintf('Test 5: Random imaginary matrices SD ...')
for n=1:1:10
    
    A=rand(n)+1i*rand(n);
    A = A*A';
    A = A + n*eye(n);
    X=rand(n,1)+1i*rand(n,1);
    B=A*X;
    sol = solveSD(A,B);
    assert(all(abs(sol - X)./X<1e-10));
end
fprintf('\tpassed\n');

fprintf('Test 6: Random imaginary matrices CG ...')
for n=1:1:10
    
    A=rand(n)+1i*rand(n);
    A = A*A';
    A = A + n*eye(n);
    X=rand(n,1)+1i*rand(n,1);
    B=A*X;
    sol = solveCG(A,B);
    assert(all(abs(sol - X)./X<1e-10));
end
fprintf('\tpassed\n');


fprintf('All tests passed!\n');