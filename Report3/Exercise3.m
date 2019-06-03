%% Clean
clear all;
close all;
clc;

%% Question 1 : Solve the system
A = [1 1; 1 2; 1 3];
b=[1;1.2;1.1;];
fprintf('\nThe solution of the system is given by :\n')
% Solve the system with pinv
x = pinv(A)*b;
fprintf('x1 : ');
disp(x(1));
fprintf('\nx2 : ');
disp(x(2));

%% Question 2: We create a grid where each element will be the value if the residual
% norm corresponding to the coordinate of the grid.
x_grid = linspace(x(1)-1,x(1)+1);
y_grid = linspace(x(2)-1,x(2)+1);
mesh_value = zeros(100, 100);
for i=1:100
    for j=1:100
        % Give value to the grid
        mesh_value(i,j) = norm(A*[x_grid(i);y_grid(j)]-b);
    end
end

figure(1);
% We use contour to plot the value of the grid
contour(x_grid, y_grid, mesh_value, 'LevelList', -30:0.3:30, 'linewidth', 1);
hold on;
% We plot the minimum found in the previous step
plot(x(1),x(2),'r*','MarkerSize',8,'linewidth',1);
% We plot the color bar to clearly see the minimum
c = colorbar;
c.Label.String = ['r(x) =  ||Ax-b||'];
set(c,'TickLabelInterpreter','latex');
legend('Value of the residual r(x)','Solution of the system corresponding to the min of r(x)')
xlabel('x_{1}');
ylabel('x_{2}');
title(['Plot of magnitude of the norm of residual r(x) in the vicinity' ...
    ' of the solution x']);
filename='./plot/minimization_pinv.eps';
print(gcf,'-depsc',filename)