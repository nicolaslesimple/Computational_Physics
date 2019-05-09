clear all;
close all;
clc;

%Example of symbolic solution of a linear system

% declare the symbolic variables
syms a x y;

% solve the system with respect to (x,y). sol is an object contains the solution (x*,y*).
%You have to provide the equations (N.B.  replace '=' by '==') and the variables.   
sol=solve(2*a*x+y==0 , 3*x+a*y==1, x, y);

disp('x=');
disp(sol.x);
disp('y=')
disp(sol.y)