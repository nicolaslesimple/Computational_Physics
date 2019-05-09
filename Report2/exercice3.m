clear all;
close all;
clc;

% Define the system
A = [0,10,50,0;-1,3,10,0;-2,10,0,0;0,0,0,1];
r = [5;0;0;1];

% Solve the system
x = A\r;

fprintf('To have a overall system at equilibrium, we need the following population 1 :')
disp(x(1))
fprintf('To have a overall system at equilibrium, we need the following population 2 :')
disp(x(2))
fprintf('To have a overall system at equilibrium, we need the following population 3 :')
disp(x(3))
fprintf('To have a overall system at equilibrium, we need the following population 4 :')
disp(x(4))