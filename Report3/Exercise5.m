clear all;
close all;
clc;

% Initialization of the variable to store the results
result = zeros(36,1);
% We iterate on each file given on moodle
for i=1:36
    % Loading of each file
    Q = csvread(char(['./wavefunctions/wavefunctions_' num2str(i) '.csv']));
    [n,~]=size(Q);
    
    % We compute the equation given in the exercise paper
    mean_k = zeros(n,1);
    for k=1:n
        O = [Q(:,k) Q(:,k+n)];
        [S,T,V] = svd(O);
        matrix = conj(S')*Q;
        mean_k(k) = (matrix(1,k)*matrix(2,k+n)-matrix(1,k+n)*matrix(2,k))^2;
    end
    result(i) = sum(mean_k)/n;
end

% Plot of the expectation value of double occupancy as a function of U/t
data = csvread('./wavefunctions/values_u.csv');
figure(1);
plot(data, result, 'LineWidth',1,'MarkerSize',3)
set(gca,'FontName','Arial','FontSize',15)
title1=title('Plot of the expectation value of double occupancy as a function of U/t');
set(title1,'FontName','Arial','FontSize',12)
xlabel('U / t','FontName','Arial','FontSize',16);
ylabel('Expectation value of double occupancy','FontName','Arial','FontSize',16);
grid on;
filename='./plot/ex5.eps';
print(gcf,'-depsc',filename)