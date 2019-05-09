clear all;
close all;
clc;

% Question 1 :
N = 3;
A=puzzleA(N);
b=[1;2;0;2;1;1;0;1;0];
x=A\b;
solution = zeros(N,N);
for i = 1:N
    solution(i,:) = x((i-1)*N+1:i*N);
end
fprintf('The solution of the puzzle is :\n')
disp(solution)


