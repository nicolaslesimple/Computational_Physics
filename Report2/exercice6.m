clear all;
close all;
clc;

%%%%%% First Part : N = 20 %%%%%%%

% Hamiltonian Matrix creation
N = 20;
H_20 = zeros(N,N);
H_20(1,N) = -1;
H_20(N,N-1) = -1;
H_20(N,1) = -1;
for i = 1:N-1
    H_20(i,i+1)=-1;
    H_20(i+1,i)=-1;
end

% Calculus of the largest eigenvalue thanks to power method.
[vec_max, val_max] = eig_ipower(H_20,10);
% Calculus of the smallest eigenvalue thanks to power method.
[vec_min, val_min] = eig_ipower(H_20,-10);

% Display the eigenvalues :
fprintf('\nThe smallest eigenvalue corresponding to the Hamiltonian of our problem is :')
disp(val_min)
fprintf('The largest eigenvalue corresponding to the Hamiltonian of our problem is :')
disp(val_max)

% Plot the two corresponding eigenstates
figure(1)
plot(vec_min, '-r', 'LineWidth', 0.1)
leg1 = legend('\psi_{1}','Location','NorthEast');
set(leg1,'FontName','Arial','FontSize',10)
title1 = title({'Plot of the eigenstates of the smallest ',...
             'eigenvalues of our Hamiltoninan with N = 20:'});
set(title1,'FontName','Arial','FontSize',12)
xlabel('Position of the atom','FontName','Arial','FontSize',10);
ylabel('\psi_{1} Eigenstate value','FontName','Arial','FontSize',10);
grid on;
filename='./plot/eigenstates_1.eps';
print(gcf,'-depsc',filename)

figure(2)
plot(vec_max, '-b', 'LineWidth', 0.1)
leg1 = legend('\psi_{N}','Location','NorthEast');
set(leg1,'FontName','Arial','FontSize',10)
title1 = title({'Plot of the eigenstates of the largest ',...
             'eigenvalues of our Hamiltoninan with N = 20:'});
set(title1,'FontName','Arial','FontSize',12)
xlabel('Position of the atom','FontName','Arial','FontSize',10);
ylabel('\psi_{N} Eigenstate value','FontName','Arial','FontSize',10);
grid on;
filename='./plot/eigenstates_N.eps';
print(gcf,'-depsc',filename)

% Plot the probility density
figure(3)
plot((vec_max.*vec_max)/norm(vec_max), '-r', 'LineWidth', 0.1)
hold on 
plot((vec_min.*vec_min)/norm(vec_min), '-b', 'LineWidth', 0.1)
leg1 = legend('Probability density of \psi_{N}','Probability density of \psi_{1}','Location','NorthEast');
set(leg1,'FontName','Arial','FontSize',10)
title1 = title({'Probability density of the Hamiltonian corresponding to the',... 
    'largest and smallest eigenvalue with N = 20'});
set(title1,'FontName','Arial','FontSize',12)
xlabel('Position of the atom','FontName','Arial','FontSize',10);
ylabel('Probability density Amplitude','FontName','Arial','FontSize',10);
grid on;
filename='./plot/proba_density.eps';
print(gcf,'-depsc',filename)


%%%%%% Second Part : N = 21 %%%%%%%

% Hamiltonian Matrix creation
N = 21;
H = zeros(N,N);
H(1,N) = -1;
H(N,N-1) = -1;
H(N,1) = -1;
for i = 1:N-1
    H(i,i+1)=-1;
    H(i+1,i)=-1;
end

% Calculus of the largest eigenvalue thanks to power method.
[vec_max, val_max] = eig_ipower(H,10);
% Calculus of the smallest eigenvalue thanks to power method.
[vec_min, val_min] = eig_ipower(H,-10);

% Display the eigenvalues :
fprintf('\nThe smallest eigenvalue corresponding to the Hamiltonian of our problem is :')
disp(val_min)
fprintf('The largest eigenvalue corresponding to the Hamiltonian of our problem is :')
disp(val_max)

% Plot the two corresponding eigenstates
figure(4)
plot(vec_min, '-r', 'LineWidth', 0.1)
leg1 = legend('\psi_{1}','Location','NorthEast');
set(leg1,'FontName','Arial','FontSize',10)
title1 = title({'Plot of the eigenstates of the smallest ',...
             'eigenvalues of our Hamiltoninan with N = 21:'});
set(title1,'FontName','Arial','FontSize',12)
xlabel('Position of the atom','FontName','Arial','FontSize',10);
ylabel('\psi_{1} Eigenstate value','FontName','Arial','FontSize',10);
grid on;
filename='./plot/eigenstates_1_21.eps';
print(gcf,'-depsc',filename)

figure(5)
plot(vec_max, '-b', 'LineWidth', 0.1)
leg1 = legend('\psi_{N}','Location','NorthEast');
set(leg1,'FontName','Arial','FontSize',10)
title1 = title({'Plot of the eigenstates of the largest ',...
             'eigenvalues of our Hamiltoninan with N = 21 :'});
set(title1,'FontName','Arial','FontSize',12)
xlabel('Position of the atom','FontName','Arial','FontSize',10);
ylabel('\psi_{N} Eigenstate value','FontName','Arial','FontSize',10);
grid on;
filename='./plot/eigenstates_N_21.eps';
print(gcf,'-depsc',filename)


% Plot the probility density
figure(6)
plot((vec_max.*vec_max)/norm(vec_max), '-r', 'LineWidth', 0.1)
hold on 
plot((vec_min.*vec_min)/norm(vec_min), '-b', 'LineWidth', 0.1)
leg1 = legend('Probability density of \psi_{N}','Probability density of the \psi_{1}','Location','NorthEast');
set(leg1,'FontName','Arial','FontSize',10)
title1 = title({'Probability density of the Hamiltonian corresponding',...
    'to the largest and smallest eigenvalue with N = 21:'});
set(title1,'FontName','Arial','FontSize',12)
xlabel('Position of the atom','FontName','Arial','FontSize',10);
ylabel('Probability density Amplitude','FontName','Arial','FontSize',10);
grid on;
filename='./plot/proba_density_21.eps';
print(gcf,'-depsc',filename)