clear all;
close all;
clc;

A_1 = [1,0,1;1.001,1,0;0,-1,1];
b_1 = [2;1;1];
x_1 = A_1\b_1;
A_2 = [1,0,1;1.001,1,0;0,-1,1];
b_2 = [2.001;1;1];
x_2 = A_2\b_2;

fprintf('The solution for the first matrix is :\n')
disp(x_1)
fprintf('The solution for the second matrix is :\n')
disp(x_2)
fprintf('The condition number for the first matrix is :\n')
disp(cond(A_1))
fprintf('The condition number for the second matrix is :\n')
disp(cond(A_2))