clear all;
close all;
clc;

% Define the parameters that characterized the problem
c=1; 
N=100; 
R=1e-9;  
a=5*R; 
dx=a/(N+1); 
e=1.60217656535e-19;
mel=9.1093829e-31; 
hbar= 1.05457172647e-34; 
V0=-1.5*e; 

H= gallery('tridiag',N^2,1,-4,1);
psi=zeros(size(H));
for i=1+N:N^2
    H(i,i-N)=1;
    H(i-N,i)=1;
end

for i=2:N
    H(1+N*(i-1), N*(i-1)) = 0;
    H( N*(i-1), 1+N*(i-1)) = 0;
end

H=-hbar^2/(2*mel*dx^2)*H;

for i=1:N
    for j=1:N
        if (i*dx-a/2)^2+((j*dx-a/2)/c)^2<R^2
            H(N*(i-1)+j,N*(i-1)+j)=H(N*(i-1)+j,N*(i-1)+j)+V0;
        end 
    end
end

% Apply eigs function  to find eigenvalues and eigenvectors
[psi, E]=eigs(H,11,'SA'); 
% Conversion of eigenvalues in eV and take the diagonal
E=diag(E/e);

% Find the number of bound solution, ie the negative eigenvalues
number = size(find(E<0));%-1*10^(-30)));
number = number(1);
fprintf('The number of bound solutions is :')
disp(number)

% Find the three lowest energy states, ie the three lowest eigenvalues 
min_eig = sort(E);
min_eig = min_eig(1:3);
fprintf('The three lowest energy states, i.e the three lowest eigenvalues are :')
disp(min_eig)
min_index_1 = find(E==min_eig(1));
min_index_2 = find(E==min_eig(2));
min_index_3 = find(E==min_eig(3));

% Find the corresponding eigenvector
vector_min_1 = psi(:,min_index_1);
vector_min_2 = psi(:,min_index_2);
vector_min_3 = psi(:,min_index_3);


% Reconstruct the matrix corresponding to each vector
matrix_min_1 = eye(N,N);
matrix_min_2 = eye(N,N);
matrix_min_3 = eye(N,N);
matrix_near_1 = eye(N,N);
x = linspace(1,N^2,N^2/N);
for i=1:N
    for j=1:N
        matrix_min_1(i,j) = vector_min_1((i-1)*N+j);
        matrix_min_2(i,j) = vector_min_2((i-1)*N+j);
        matrix_min_3(i,j) = vector_min_3((i-1)*N+j);
    end
end

% Plot the results
plot_pcolor(matrix_min_1,'Plot of \psi_{1} wavefunctions (eigenvectors)','Plot of \psi_{1} probability densities','./plot/pcolor_1.eps','./plot/pcolor_proba_1.eps')
plot_pcolor(matrix_min_2,'Plot of \psi_{2} wavefunctions (eigenvectors)','Plot of \psi_{2} probability densities','./plot/pcolor_2.eps','./plot/pcolor_proba_2.eps')
plot_pcolor(matrix_min_3,'Plot of \psi_{3} wavefunctions (eigenvectors)','Plot of \psi_{3} probability densities','./plot/pcolor_3.eps','./plot/pcolor_proba_3.eps')

% Investigate one of the scattering solutions at a positive energy
[psi, E]=eigs(H,1000,'SA');
E=diag(E/e);
vector_near_1 = psi(:,54);
for i=1:N
    for j=1:N
        matrix_near_1(i,j) = vector_near_1((i-1)*N+j);
    end
end
plot_pcolor(matrix_near_1,'Plot of \psi near 1eV wavefunctions (eigenvectors)','Plot of \psi near 1eV probability densities','./plot/pcolor_1_1.eps','./plot/pcolor_proba_1_1.eps')

% Calculation of the probability y of finding particle inside the
% quantum well
proba_1 = 0;
proba_min_1 = 0;
proba_min_2 = 0;
proba_min_3 = 0;
for i=1:N
    for j=1:N
        if (i*dx-a/2)^2+((j*dx-a/2)/c)^2<R^2
            proba_1 = proba_1+ + (matrix_near_1(i,j)*matrix_near_1(i,j))/norm(matrix_near_1)^2;
            proba_min_1 = proba_min_1+ + (matrix_min_1(i,j)*matrix_min_1(i,j))/norm(matrix_min_1)^2;
            proba_min_2 = proba_min_2+ + (matrix_min_2(i,j)*matrix_min_2(i,j))/norm(matrix_min_2)^2;
            proba_min_3 = proba_min_3+ + (matrix_min_3(i,j)*matrix_min_3(i,j))/norm(matrix_min_3)^2;
        end 
    end
end
fprintf('Probability of finding particle inside the quantum well for energy of 1ev is :')
disp(proba_1*100)
fprintf('Probability of finding particle inside the quantum well for the first bound solution is :')
disp(proba_min_1*100)
fprintf('Probability of finding particle inside the quantum well for the second bound solution is :')
disp(proba_min_2*100)
fprintf('Probability of finding particle inside the quantum well for the third bound solution is :')
disp(proba_min_3*100)

% Question 3 : number of bound states as a function of c
bound = eye(1,16);
count = 0;
for c=linspace(2,10,16)
    count = count + 1;
    H= gallery('tridiag',N^2,1,-4,1);
    psi=zeros(size(H));
    for i=1+N:N^2
        H(i,i-N)=1;
        H(i-N,i)=1;
    end

    for i=2:N
        H(1+N*(i-1), N*(i-1)) = 0;
        H( N*(i-1), 1+N*(i-1)) = 0;
    end

    H=-hbar^2/(2*mel*dx^2)*H;

    for i=1:N
        for j=1:N
            if (i*dx-a/2)^2+((j*dx-a/2)/c)^2<R^2
                H(N*(i-1)+j,N*(i-1)+j)=H(N*(i-1)+j,N*(i-1)+j)+V0;
            end 
        end
    end
    % Apply eigs function  to find eigenvalues and eigenvectors
    [psi, E]=eigs(H,500,'SA'); 
    % Conversion of eigenvalues in eV and take the diagonal
    E=diag(E/e);
    % Find the number of bound solution, ie the negative eigenvalues
    number = size(find(E<0));
    bound(count) = number(1);
    disp(c)
    
    if (c == 2 || c ==3.6000 || c ==10)
        % Find the three lowest energy states, ie the three lowest eigenvalues 
        min_eig = sort(E);
        min_eig = min_eig(1:3);
        min_index_1 = find(E==min_eig(1));
        min_index_2 = find(E==min_eig(2));
        min_index_3 = find(E==min_eig(3));

        % Find the corresponding eigenvector
        vector_min_1 = psi(:,min_index_1);
        vector_min_2 = psi(:,min_index_2);
        vector_min_3 = psi(:,min_index_3);


        % Reconstruct the matrix corresponding to each vector
        matrix_min_1 = eye(N,N);
        matrix_min_2 = eye(N,N);
        matrix_min_3 = eye(N,N);
        matrix_near_1 = eye(N,N);
        x = linspace(1,N^2,N^2/N);
        for i=1:N
            for j=1:N
                matrix_min_1(i,j) = vector_min_1((i-1)*N+j);
                matrix_min_2(i,j) = vector_min_2((i-1)*N+j);
                matrix_min_3(i,j) = vector_min_3((i-1)*N+j);
            end
        end

        % Plot the results
        plot_pcolor(matrix_min_1,['Plot of \psi_{1} wavefunctions (eigenvectors) for c = ', num2str(c)],['Plot of \psi_{1} probability densities for c = ',num2str(c)],['./plot/pcolor_1' num2str(c) '.eps'],['./plot/pcolor_proba_1_',num2str(c),'.eps'])
        plot_pcolor(matrix_min_2,['Plot of \psi_{2} wavefunctions (eigenvectors) for c = ', num2str(c)],['Plot of \psi_{2} probability densities for c = ',num2str(c)],['./plot/pcolor_2' num2str(c) '.eps'],['./plot/pcolor_proba_2_',num2str(c),'.eps'])
        plot_pcolor(matrix_min_3,['Plot of \psi_{3} wavefunctions (eigenvectors) for c = ', num2str(c)],['Plot of \psi_{3} probability densities for c = ',num2str(c)],['./plot/pcolor_3' num2str(c) '.eps'],['./plot/pcolor_proba_3_',num2str(c),'.eps'])
        
    end
end
    

figure
plot(linspace(2,10,16), bound, '-b', 'LineWidth', 0.3)
title1 = title({'Number of bound states as a function of c'});
set(title1,'FontName','Arial','FontSize',12)
xlabel('c [-]','FontName','Arial','FontSize',10);
ylabel('Number of bound states [-]','FontName','Arial','FontSize',10);
grid on;
filename='./plot/c.eps';
print(gcf,'-depsc',filename)