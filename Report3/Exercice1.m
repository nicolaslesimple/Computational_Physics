clear all;
close all;
clc;

%% Question 1 :
N = 3;
A=puzzleA(N);
b=[1;2;0;2;1;1;0;1;0];
x=A\b;
solution = zeros(N,N);
for i = 1:N
    solution(i,:) = x((i-1)*N+1:i*N);
end
fprintf('The solution of the puzzle is :\n')
disp(abs(solution))


%% Question 3
boards = importdata('./boards.mat');
board1 = boards.board1;
board2 = boards.board2;
boards_size_1 = size(board1);
boards_size_2 = size(board2);
N_1 = boards_size_1(1);
N_2 = boards_size_2(1);
A_1 = puzzleA(N_1);
A_2 = puzzleA(N_2);
b_1 = zeros(N_1^2,1);
b_2 = zeros(N_2^2,1);
% Put the input matrix in a vector format
for i = 1:N_1
    b_1((i-1)*N_1+1:i*N_1) = board1(i,:);
end
for i = 1:N_2
    b_2((i-1)*N_2+1:i*N_2) = board2(i,:);
end

% Solving system with different technics for board1
sol_1 = ((A_1*A_1)\(A_1*b_1));
x_1_SD = (solveSD(A_1*A_1,A_1*b_1));
x_1_CG = (solveCG(A_1*A_1,A_1*b_1));
% Solving system with different technics for board2
sol_2 = ((A_2*A_2)\(A_2*b_2));
x_2_SD = (solveSD(A_2*A_2,A_2*b_2));
x_2_CG = (solveCG(A_2*A_2,A_2*b_2));

% Initialisation of puzzle solution for board1
solution_1 = zeros(N_1,N_1);
sol_1_SD = zeros(N_1,N_1);
sol_1_CG = zeros(N_1,N_1);
% Initialisation of puzzle solution for board2
solution_2 = zeros(N_2,N_2);
sol_2_SD = zeros(N_2,N_2);
sol_2_CG = zeros(N_2,N_2);
% Reformating of the solution in a matrix format
for i = 1:N_1
    solution_1(i,:) = sol_1((i-1)*N_1+1:i*N_1);
    sol_1_SD(i,:) = x_1_SD((i-1)*N_1+1:i*N_1);
    sol_1_CG(i,:) = x_1_CG((i-1)*N_1+1:i*N_1);
end
for i = 1:N_2
    solution_2(i,:) = sol_2((i-1)*N_2+1:i*N_2);
    sol_2_SD(i,:) = x_2_SD((i-1)*N_2+1:i*N_2);
    sol_2_CG(i,:) = x_2_CG((i-1)*N_2+1:i*N_2);
end

figure(1)
imagesc(board1);            % Create a colored plot of the matrix values
colorbar;
textStrings = num2str(board1(:), '%0.2f');       % Create strings from the matrix values
textStrings = strtrim(cellstr(textStrings));  % Remove any space padding
[x, y] = meshgrid(1:10);  % Create x and y coordinates for the strings
hStrings = text(x(:), y(:), textStrings(:),'HorizontalAlignment', 'center');  % Plot the strings       
midValue = mean(get(gca, 'CLim'));  % Get the middle value of the color range
set(gca, 'XTick', 1:11,'YTick', 1:11, 'TickLength', [0 0]);
title1 = title({'Plot of the first puzzle called board1 :'});
set(title1,'FontName','Arial','FontSize',12)
xlabel('Index of the matrix','FontName','Arial','FontSize',10);
ylabel('Index of the matrix','FontName','Arial','FontSize',10);
grid on;
filename='./plot/board_1.eps';
print(gcf,'-depsc',filename)

figure(2)
imagesc(sol_1_SD);            % Create a colored plot of the matrix values
colorbar;
textStrings = num2str(sol_1_SD(:), '%0.2f');       % Create strings from the matrix values
textStrings = strtrim(cellstr(textStrings));  % Remove any space padding
[x, y] = meshgrid(1:10);  % Create x and y coordinates for the strings
hStrings = text(x(:), y(:), textStrings(:),'HorizontalAlignment', 'center');  % Plot the strings                
midValue = mean(get(gca, 'CLim'));  % Get the middle value of the color range
set(gca, 'XTick', 1:11,'YTick', 1:11, 'TickLength', [0 0]);
title1 = title({'Plot of the solution of the puzzle called board1 obtain with',...
             'Steepest Descent :'});
set(title1,'FontName','Arial','FontSize',12)
xlabel('Index of the solution matrix','FontName','Arial','FontSize',10);
ylabel('Index of the solution matrix','FontName','Arial','FontSize',10);
grid on;
filename='./plot/board_1_sg.eps';
print(gcf,'-depsc',filename)

figure(3)
imagesc(sol_1_CG);            % Create a colored plot of the matrix values
colorbar;
textStrings = num2str(sol_1_CG(:), '%0.2f');       % Create strings from the matrix values
textStrings = strtrim(cellstr(textStrings));  % Remove any space padding
[x, y] = meshgrid(1:10);  % Create x and y coordinates for the strings
hStrings = text(x(:), y(:), textStrings(:),'HorizontalAlignment', 'center');  % Plot the strings
midValue = mean(get(gca, 'CLim'));  % Get the middle value of the color range
set(gca, 'XTick', 1:11,'YTick', 1:11, 'TickLength', [0 0]);
title1 = title({'Plot of the solution of the puzzle called board1 obtain with',...
             ' Conjugate Gradient method:'});
set(title1,'FontName','Arial','FontSize',12)
xlabel('Index of the solution matrix','FontName','Arial','FontSize',10);
ylabel('Index of the solution matrix','FontName','Arial','FontSize',10);
grid on;
filename='./plot/board_1_cg.eps';
print(gcf,'-depsc',filename)

figure(4)
imagesc(board2);            % Create a colored plot of the matrix values
colorbar;
textStrings = num2str(board2(:), '%0.2f');       % Create strings from the matrix values
textStrings = strtrim(cellstr(textStrings));  % Remove any space padding
[x, y] = meshgrid(1:10);  % Create x and y coordinates for the strings
hStrings = text(x(:), y(:), textStrings(:),'HorizontalAlignment', 'center');  % Plot the strings             
midValue = mean(get(gca, 'CLim'));  % Get the middle value of the color range
set(gca, 'XTick', 1:11,'YTick', 1:11, 'TickLength', [0 0]);
title1 = title({'Plot of the second puzzle called board2 :'});
set(title1,'FontName','Arial','FontSize',12)
xlabel('Index of the matrix','FontName','Arial','FontSize',10);
ylabel('Index of the matrix','FontName','Arial','FontSize',10);
grid on;
filename='./plot/board_2.eps';
print(gcf,'-depsc',filename)

figure(5)
imagesc(sol_2_SD);            % Create a colored plot of the matrix values
colorbar;
textStrings = num2str(sol_2_SD(:), '%0.2f');       % Create strings from the matrix values
textStrings = strtrim(cellstr(textStrings));  % Remove any space padding
[x, y] = meshgrid(1:10);  % Create x and y coordinates for the strings
hStrings = text(x(:), y(:), textStrings(:),'HorizontalAlignment', 'center');  % Plot the strings                
midValue = mean(get(gca, 'CLim'));  % Get the middle value of the color range
set(gca, 'XTick', 1:11,'YTick', 1:11, 'TickLength', [0 0]);
title1 = title({'Plot of the solution of the puzzle called board2 obtain with',...
             'Steepest Descent :'});
set(title1,'FontName','Arial','FontSize',12)
xlabel('Index of the solution matrix','FontName','Arial','FontSize',10);
ylabel('Index of the solution matrix','FontName','Arial','FontSize',10);
grid on;
filename='./plot/board_2_sg.eps';
print(gcf,'-depsc',filename)

figure(6)
imagesc(sol_2_CG);            % Create a colored plot of the matrix values
colorbar;
textStrings = num2str(sol_2_CG(:), '%0.2f');       % Create strings from the matrix values
textStrings = strtrim(cellstr(textStrings));  % Remove any space padding
[x, y] = meshgrid(1:10);  % Create x and y coordinates for the strings
hStrings = text(x(:), y(:), textStrings(:),'HorizontalAlignment', 'center');  % Plot the strings
midValue = mean(get(gca, 'CLim'));  % Get the middle value of the color range
set(gca, 'XTick', 1:11,'YTick', 1:11, 'TickLength', [0 0]);
title1 = title({'Plot of the solution of the puzzle called board2 obtain with',...
             ' Conjugate Gradient method:'});
set(title1,'FontName','Arial','FontSize',12)
xlabel('Index of the solution matrix','FontName','Arial','FontSize',10);
ylabel('Index of the solution matrix','FontName','Arial','FontSize',10);
grid on;
filename='./plot/board_2_cg.eps';
print(gcf,'-depsc',filename)

fprintf('We check if our solution makes sense with solution found with the build in Matlab function: \n')
fprintf('Error for Steepest Descent method for board1: ')
disp(sum(sum(abs(solution_1-sol_1_SD))));
fprintf('Error for Conjugate Gradient method for board1: ')
disp(sum(sum(abs(solution_1-sol_1_CG))));
fprintf('Error for Steepest Descent for board2: ')
disp(sum(sum(abs(solution_2-sol_2_SD))));
fprintf('Error for Conjugate Gradient method for board2: ')
disp(sum(sum(abs(solution_2-sol_2_CG))));

%% Question 4 :
% Redefine the problem and the perfect solution
boards_size_1 = size(board1);
N_1 = boards_size_1(1);
A_1 = puzzleA(N_1);
b_1 = zeros(N_1^2,1);
for i = 1:N_1
    b_1((i-1)*N_1+1:i*N_1) = board1(i,:);
end
sol_1 = ((A_1*A_1)\(A_1*b_1));
solution_1 = zeros(N_1,N_1);
for i = 1:N_1
    solution_1(i,:) = sol_1((i-1)*N_1+1:i*N_1);
end

list_error_sd_1 = zeros(1,5e3);
list_error_cg_1 = zeros(1,5e3);
for j = 1:5e3
    x_SD = (solveSD(A_1*A_1,A_1*b_1,j));
    x_CG = (solveCG(A_1*A_1,A_1*b_1,j));
    sol_SD = zeros(N_1,N_1);
    sol_CG = zeros(N_1,N_1);
    for i = 1:N_1
        sol_SD(i,:) = x_SD((i-1)*N_1+1:i*N_1);
        sol_CG(i,:) = x_CG((i-1)*N_1+1:i*N_1);
    end
    list_error_sd_1(1,j) = mean(mean(abs(sol_SD - solution_1)));
    list_error_cg_1(1,j) = mean(mean(abs(sol_CG - solution_1)));
end 

figure(7)
semilogy( list_error_sd_1, '-b', 'LineWidth', 0.1)
hold on
semilogy( list_error_cg_1, '-r', 'LineWidth', 0.1)
leg1 = legend('Steepest Descent error','Conjugate Gradient error','Location','NorthEast');
set(leg1,'FontName','Arial','FontSize',10)
title1 = title({'Plot of the error for the board 1 problem obtained with Steepest Descent ',...
             'and Conugate Gradient method in function of the number of iteration:'});
set(title1,'FontName','Arial','FontSize',12)
xlabel('Number of iteration allowed','FontName','Arial','FontSize',10);
ylabel('Relative error [%]','FontName','Arial','FontSize',10);
grid on;
hold off;
filename='./plot/error_1.eps';
print(gcf,'-depsc',filename)

% Redefine the problem and the perfect solution for the second case
boards_size_2 = size(board2);
N_2 = boards_size_1(1);
A_2 = puzzleA(N_2);
b_2 = zeros(N_2^2,1);
for i = 1:N_2
    b_2((i-1)*N_2+1:i*N_2) = board2(i,:);
end
sol_2 = ((A_2*A_2)\(A_2*b_2));
solution_2 = zeros(N_2,N_2);
for i = 1:N_2
    solution_2(i,:) = sol_2((i-1)*N_2+1:i*N_2);
end

list_error_sd_2 = zeros(1,5e3);
list_error_cg_2 = zeros(1,5e3);
for j = 1:5e3
    x_SD = (solveSD(A_2*A_2,A_2*b_2,j));
    x_CG = (solveCG(A_2*A_2,A_2*b_2,j));
    sol_SD = zeros(N_2,N_2);
    sol_CG = zeros(N_2,N_2);
    for i = 1:N_2
        sol_SD(i,:) = x_SD((i-1)*N_2+1:i*N_2);
        sol_CG(i,:) = x_CG((i-1)*N_2+1:i*N_2);
    end
    list_error_sd_2(1,j) = mean(mean(abs(sol_SD - solution_2)));
    list_error_cg_2(1,j) = mean(mean(abs(sol_CG - solution_2)));
end 

figure(8)
semilogy( list_error_sd_2, '-b', 'LineWidth', 0.1)
hold on
semilogy( list_error_cg_2, '-r', 'LineWidth', 0.1)
leg1 = legend('Steepest Descent error','Conjugate Gradient error','Location','NorthEast');
set(leg1,'FontName','Arial','FontSize',10)
title1 = title({'Plot of the error for the board 2 problem obtained with Steepest Descent ',...
             'and Conugate Gradient method in function of the number of iteration:'});
set(title1,'FontName','Arial','FontSize',12)
xlabel('Number of iteration allowed','FontName','Arial','FontSize',10);
ylabel('Relative error [%]','FontName','Arial','FontSize',10);
grid on;
hold off;
filename='./plot/error_2.eps';
print(gcf,'-depsc',filename)

figure(9)
semilogy( list_error_sd_1, '-b', 'LineWidth', 0.1)
hold on
semilogy( list_error_sd_2, '-r', 'LineWidth', 0.1)
leg1 = legend('Steepest Descent error of board 1','Steepest Descent error of board 2','Location','NorthEast');
set(leg1,'FontName','Arial','FontSize',10)
title1 = title({'Plot of the error for the board 1 and 2 problem obtained with Steepest Descent ',...
             'in function of the number of iteration:'});
set(title1,'FontName','Arial','FontSize',12)
xlabel('Number of iteration allowed','FontName','Arial','FontSize',10);
ylabel('Relative error [%]','FontName','Arial','FontSize',10);
grid on;
hold off;
filename='./plot/error_sg.eps';
print(gcf,'-depsc',filename)

figure(10)
semilogy( list_error_cg_1, '-b', 'LineWidth', 0.1)
hold on
semilogy( list_error_cg_2, '-r', 'LineWidth', 0.1)
leg1 = legend('Conjugate Gradient error for board 1','Conjugate Gradient error for board 2','Location','NorthEast');
set(leg1,'FontName','Arial','FontSize',10)
title1 = title({'Plot of the error for the board 1 and 2 problem obtained with ',...
             'Conjugate Gradient method in function of the number of iteration:'});
set(title1,'FontName','Arial','FontSize',12)
xlim([0 100])
xlabel('Number of iteration allowed','FontName','Arial','FontSize',10);
ylabel('Relative error [%]','FontName','Arial','FontSize',10);
grid on;
hold off;
filename='./plot/error_cg.eps';
print(gcf,'-depsc',filename)

