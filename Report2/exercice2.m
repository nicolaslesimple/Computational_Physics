clear all;
close all;
clc;

% Definition of the system
A = [0,9/12,0,0,0,0;
     1/2,0,0,0,0,0;
     0,0,3/14,10/14,0,0;
     0,0,5/14,12/14,0,0;
     0,0,0,0,3/32,7/32;
     6/32,0,0,0,4/32,0];
 b=[6;2;4;6;25/8;15/4];
 
 % We solve the system thus we find the mass per object
 x = A \ b;
 fprintf('The different mass are stores in the following vecotor :\n')
 disp(x);
 fprintf('Thus, we have :')
 fprintf('m_alpha = ')
 disp(x(1))
 fprintf('m_beta = ')
 disp(x(2))
 fprintf('m_gamma = ')
 disp(x(3))
 fprintf('m_epsilon = ')
 disp(x(4))
 fprintf('m_delta = ')
 disp(x(5))
 fprintf('m_zeta = ')
 disp(x(6))
 
 % We find the center of mass of the overall system using the formula
 R_x = (x(1)*0+x(2)*9+x(3)*3+x(4)*10+x(5)*3+x(6)*7)/(sum(x));
 R_y = (x(1)*6+x(2)*0+x(3)*5+x(4)*12+x(5)*4+x(6)*0)/(sum(x));
 fprintf('The center of mass of the overall system has the following coordinates :')
 disp([R_x,R_y])
 