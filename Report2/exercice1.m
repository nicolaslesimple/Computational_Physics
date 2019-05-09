clear all;
close all;
clc;

% Question 1 : Verify the solution explicitly by
% comparing with the one obtained using Matlab?s matrix operations
A = [1,1,1,1; 2,3,0,-1; -3,4,1,2; 1,2,-1,1];
b = [13;-1;10;1];
x_matlab = A\b;
x = solve_ls(A,b);
fprintf('Here is the x results of the equation Ax = b, using matlab function :\n')
disp(x_matlab);
fprintf('Here is the x results of the equation Ax = b, using our implemented function :\n')
disp(x_matlab);
fprintf('Here is the mean  difference between the previous output :\n')
disp(mean(x-x_matlab));



% Question 2 : Implement the LU decomposition of a square matrix. 
% Test your implementation on random matrices of size up to 100.
diff_L = zeros(20,1);
diff_U = zeros(20,1);
diff_LU = zeros(20,1);
for i = 1:20
    m = randi(100,5*i,5*i);
    [L_matlab, U_matlab, P] = lu(m);
    LU_matlab = L_matlab*U_matlab;
    % In order to compare properly, run your LU function on A0 = P*A, 
    % where P is taken from Matlab's output.
    m_matlab = P*m;
    [L_me, U_me] = LU_square_matrix(m_matlab);
    LU_me = L_me*U_me;
    size_tmp = size(L_me);
    tmp_L = 0; 
    tmp_U = 0;
    tmp_LU= 0;
    for k = 1: size_tmp(1)
        for j = 1: size_tmp(2)
            tmp_L = tmp_L + (L_me(k,j)^2-L_matlab(k,j)^2);
            tmp_U = tmp_U + (U_me(k,j)^2-U_matlab(k,j)^2);
            tmp_LU = tmp_LU + (LU_me(k,j)^2-LU_matlab(k,j)^2);
        end
    end
    diff_L(i) = abs(tmp_L);
    diff_U(i) = abs(tmp_U);
    diff_LU(i) = abs(tmp_LU);
end

% Question 2 : Compare your implementation with Matlab's lu function 
% plotting the difference between the solutions provided by the two methods.
figure(1)
x = linspace(5,100,20);
plot(x, diff_L, '-r', 'LineWidth', 0.1)
leg1 = legend('L matrix Difference','Location','NorthEast');
set(leg1,'FontName','Arial','FontSize',10)
title1 = title({'Plot of the difference between L matrices obtained ',...
             'with our own implementation and with Matlab build in funciton :'});
set(title1,'FontName','Arial','FontSize',12)
xlabel('Size of the square random matrix (n*n)','FontName','Arial','FontSize',10);
ylabel('Difference [-]','FontName','Arial','FontSize',10);
grid on;
filename='./plot/LU_decomposition_square_matrix_comparison_L.eps';
print(gcf,'-depsc',filename)

figure(2)
plot(x, diff_U, '-b', 'LineWidth', 0.1)
leg1 = legend('U matrix Difference','Location','NorthEast');
set(leg1,'FontName','Arial','FontSize',10)
title1 = title({'Plot of the difference between U matrices obtained ',...
             'with our own implementation and with Matlab build in funciton :'});
set(title1,'FontName','Arial','FontSize',12)
xlabel('Size of the square random matrix (n*n)','FontName','Arial','FontSize',10);
ylabel('Difference [-]','FontName','Arial','FontSize',10);
grid on;
filename='./plot/LU_decomposition_square_matrix_comparison_U.eps';
print(gcf,'-depsc',filename)

figure(3)
plot(x, diff_LU, '-g', 'LineWidth', 0.1)
leg1 = legend('LU matrix Difference','Location','NorthEast');
set(leg1,'FontName','Arial','FontSize',10)
title1 = title({'Plot of the difference between LU matrices obtained ',...
             'with our own implementation and with Matlab build in funciton :'});
set(title1,'FontName','Arial','FontSize',12)
xlabel('Size of the square random matrix (n*n)','FontName','Arial','FontSize',10);
ylabel('Difference [-]','FontName','Arial','FontSize',10);
grid on;
filename='./plot/LU_decomposition_square_matrix_comparison_LU.eps';
print(gcf,'-depsc',filename)

% Question 3 : Test your implementation on the following matrix:
test_matrix = [0,-7,0;-3,2,6;5,-1,5];
[L_matlab, U_matlab, P] = lu(test_matrix);
LU_matlab = L_matlab * U_matlab;
[L_me, U_me] = LU_square_matrix(test_matrix);
LU_me = L_me * U_me;
size_tmp = size(L_me);
tmp_L = 0; 
tmp_U = 0;
tmp_LU= 0;
for k = 1: size_tmp(1)
    for j = 1: size_tmp(2)
        tmp_L = tmp_L + (L_me(k,j)^2-L_matlab(k,j)^2);
        tmp_U = tmp_U + (U_me(k,j)^2-U_matlab(k,j)^2);
        tmp_LU = tmp_LU + (LU_me(k,j)^2-LU_matlab(k,j)^2);
    end
end
diff_L = abs(tmp_L);
diff_U = abs(tmp_U);
diff_LU = abs(tmp_LU);
fprintf('Question 3 : \nThe diffenrence between the LU decomposition is :')
disp(diff_LU)
fprintf('The diffenrence between the U matrices is :')
disp(diff_U)
fprintf('The diffenrence between the L matrices is :')
disp(diff_L)

% Question 4 :Implement the LU decomposition with pivoting. Test your 
% implementation of LU decomposition on the above matrix.
[L_me,U_me,P_me] = lu_decomposition(test_matrix);
LU_me = inv(P_me)*L_me * U_me;
[L_matlab, U_matlab, P_matlab] = lu(test_matrix);
LU_matlab = inv(P_matlab)*L_matlab * U_matlab;
size_tmp = size(L_me);
tmp_L = 0; 
tmp_U = 0;
tmp_LU= 0;
tmp_P = 0;
for k = 1: size_tmp(1)
    for j = 1: size_tmp(2)
        tmp_L = tmp_L + (L_me(k,j)^2-L_matlab(k,j)^2);
        tmp_U = tmp_U + (U_me(k,j)^2-U_matlab(k,j)^2);
        tmp_P = tmp_P + (P_me(k,j)^2-P_matlab(k,j)^2);
        tmp_LU = tmp_LU + (LU_me(k,j)^2-LU_matlab(k,j)^2);
    end
end
diff_L = abs(tmp_L);
diff_U = abs(tmp_U);
diff_LU = abs(tmp_LU);
diff_P = abs(tmp_P);
fprintf('Question 4 : \nThe difference between the LU decomposition is :')
disp(diff_LU)
fprintf('The difference between the U matrices is :')
disp(diff_U)
fprintf('The difference between the L matrices is :')
disp(diff_L)
fprintf('The difference between the P matrices is :')
disp(diff_P)
